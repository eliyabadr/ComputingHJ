// cell.hpp - hyperrectangular cells and the adaptively refined cell tree.
//
// Cells live in one arena (`cells_`) and are addressed by id; `leaves_` holds
// the ids of the current leaves in the same order the Python list keeps them.
#pragma once

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <ctime>
#include <limits>
#include <vector>

#include "env.hpp"
#include "rtree.hpp"
#include "timing.hpp"

namespace dhj {

struct Cell {
    Bounds3 bounds{};
    State3 center{};
    std::uint32_t cell_id = 0;
    std::int64_t parent = -1;
    std::int64_t child0 = -1, child1 = -1;
    bool is_leaf = true;

    double V_upper = 0.0, V_lower = 0.0;
    double l_upper = 0.0, l_lower = 0.0;
    double r_upper = 0.0, r_lower = 0.0;

    double range(int d) const { return bounds[d][1] - bounds[d][0]; }

    // Dimension of maximum range; ties (within 1e-10) go to the lowest index,
    // matching Cell.get_max_range_dim.
    int max_range_dim() const {
        double mx = range(0);
        for (std::size_t d = 1; d < kDim; ++d) mx = std::max(mx, range(static_cast<int>(d)));
        constexpr double tol = 1e-10;
        for (std::size_t d = 0; d + 1 < kDim; ++d)
            if (std::fabs(range(static_cast<int>(d)) - mx) < tol) return static_cast<int>(d);
        return static_cast<int>(kDim - 1);
    }

    double max_range() const { return range(max_range_dim()); }

    bool contains_point(const State3& p) const {
        for (std::size_t d = 0; d < kDim; ++d)
            if (p[d] < bounds[d][0] || p[d] > bounds[d][1]) return false;
        return true;
    }
};

inline State center_of(const Bounds& b) {
    State c{};
    for (std::size_t d = 0; d < kDim; ++d) c[d] = 0.5 * (b[d][0] + b[d][1]);
    return c;
}

class CellTree {
public:
    CellTree(const Bounds3& root_bounds, int initial_resolution) : root_bounds_(root_bounds) {
        create_initial_grid(initial_resolution);
        build_spatial_index();
    }

    // Restore constructor, used when resuming from a checkpoint. The caller
    // supplies an arena and leaf list that were already consistent when saved;
    // only the spatial index has to be rebuilt.
    CellTree(const Bounds3& root_bounds, std::vector<Cell>&& cells,
             std::vector<std::uint32_t>&& leaves)
        : root_bounds_(root_bounds), cells_(std::move(cells)), leaves_(std::move(leaves)) {
        build_spatial_index();
    }

    // `index_` stores a pointer into `leaf_boxes_`, a sibling member, so a
    // copy or move would leave the new object's index aimed at the old
    // object's vector -- a dangling read as soon as the source dies. The tree
    // is only ever owned through unique_ptr, so forbid both outright rather
    // than silently rebuilding the index.
    CellTree(const CellTree&) = delete;
    CellTree& operator=(const CellTree&) = delete;
    CellTree(CellTree&&) = delete;
    CellTree& operator=(CellTree&&) = delete;

    const Bounds3& root_bounds() const { return root_bounds_; }
    const std::vector<std::uint32_t>& leaves() const { return leaves_; }
    std::size_t num_leaves() const { return leaves_.size(); }
    Cell& cell(std::uint32_t id) { return cells_[id]; }
    const Cell& cell(std::uint32_t id) const { return cells_[id]; }
    std::size_t num_cells() const { return cells_.size(); }

    // Splits every cell in `to_refine` (given in leaf order) at the midpoint of
    // its widest dimension. Leaf ordering matches CellTree.refine_cell called in
    // a loop: surviving leaves keep their order, then the children are appended
    // in refinement order.
    std::vector<std::uint32_t> refine_cells(const std::vector<std::uint32_t>& to_refine) {
        std::vector<bool> refining(cells_.size(), false);
        for (std::uint32_t id : to_refine) refining[id] = true;

        std::vector<std::uint32_t> new_cells;
        new_cells.reserve(2 * to_refine.size());
        for (std::uint32_t id : to_refine) {
            Cell& parent = cells_[id];
            if (!parent.is_leaf) continue;
            const int d = parent.max_range_dim();
            const double mid = 0.5 * (parent.bounds[d][0] + parent.bounds[d][1]);

            Cell c0{}, c1{};
            c0.bounds = parent.bounds; c0.bounds[d][1] = mid;
            c1.bounds = parent.bounds; c1.bounds[d][0] = mid;
            c0.center = center_of(c0.bounds);
            c1.center = center_of(c1.bounds);
            c0.cell_id = static_cast<std::uint32_t>(cells_.size());
            c1.cell_id = c0.cell_id + 1;
            c0.parent = c1.parent = static_cast<std::int64_t>(id);

            cells_[id].is_leaf = false;
            cells_[id].child0 = c0.cell_id;
            cells_[id].child1 = c1.cell_id;
            cells_.push_back(c0);
            cells_.push_back(c1);
            new_cells.push_back(c0.cell_id);
            new_cells.push_back(c1.cell_id);
        }

        std::vector<std::uint32_t> kept;
        kept.reserve(leaves_.size() + new_cells.size());
        for (std::uint32_t id : leaves_)
            if (id >= refining.size() || !refining[id]) kept.push_back(id);
        kept.insert(kept.end(), new_cells.begin(), new_cells.end());
        leaves_.swap(kept);
        return new_cells;
    }

    void rebuild_spatial_index() { build_spatial_index(); }

    // Leaves intersecting `b`, with theta compared modulo 2*pi as in
    // CellTree.get_intersecting_cells. `fn` is called once per distinct leaf.
    //
    // `stamp` (one slot per leaf) and `epoch` do the de-duplication the Python
    // gets from a Python set. The caller owns `epoch` and bumps it whenever a
    // new de-duplication group starts, so several reach boxes belonging to the
    // same (cell, action) pair can share one group.
    template <typename Fn>
    void for_each_intersecting(const Bounds3& b, std::vector<std::uint32_t>& stamp,
                               std::uint32_t epoch, Fn&& fn) const {
        const double shifts[3] = {0.0, -2.0 * kPi, 2.0 * kPi};
        for (double sh : shifts) {
            // only the periodic dimension wraps; every other axis is taken as-is
            double lo[kDim], hi[kDim];
            for (std::size_t d = 0; d < kDim; ++d) {
                const double off = (d == kPeriodicDim) ? sh : 0.0;
                lo[d] = b[d][0] + off;
                hi[d] = b[d][1] + off;
            }
            index_.query(lo, hi, [&](std::size_t entry) {
                const std::uint32_t leaf_pos = static_cast<std::uint32_t>(entry);
                if (stamp[leaf_pos] == epoch) return;
                stamp[leaf_pos] = epoch;
                fn(leaves_[leaf_pos]);
            });
        }
    }

    // Leaf containing `p`, or -1. Used to rasterize theta slices for plotting.
    // A point on a shared face lies in both neighbours; the lowest cell id wins
    // so the result does not depend on index traversal order (and so the C++
    // and Rust rasterizers agree exactly).
    std::int64_t find_leaf_containing(const State3& p) const {
        double lo[kDim];
        for (std::size_t d = 0; d < kDim; ++d) lo[d] = p[d];
        std::int64_t found = -1;
        index_.query(lo, lo, [&](std::size_t entry) {
            const std::int64_t id = static_cast<std::int64_t>(leaves_[entry]);
            if (found < 0 || id < found) found = id;
        });
        return found;
    }

    double index_build_seconds() const { return index_build_seconds_; }

private:
    void create_initial_grid(int resolution) {
        std::vector<std::vector<double>> edges(kDim);
        for (std::size_t d = 0; d < kDim; ++d) {
            edges[d].resize(static_cast<std::size_t>(resolution) + 1);
            const double a = root_bounds_[d][0], b = root_bounds_[d][1];
            for (int i = 0; i <= resolution; ++i)
                edges[d][static_cast<std::size_t>(i)] =
                    a + (b - a) * static_cast<double>(i) / static_cast<double>(resolution);
            edges[d][static_cast<std::size_t>(resolution)] = b;  // exact endpoint, as np.linspace
        }

        std::size_t total = 1;
        for (std::size_t d = 0; d < kDim; ++d) total *= static_cast<std::size_t>(resolution);
        cells_.reserve(total);
        leaves_.reserve(total);

        // itertools.product ordering: last dimension varies fastest.
        std::array<std::size_t, kDim> idx{};
        for (std::size_t n = 0; n < total; ++n) {
            Cell c{};
            for (std::size_t d = 0; d < kDim; ++d)
                c.bounds[d] = {{edges[d][idx[d]], edges[d][idx[d] + 1]}};
            c.center = center_of(c.bounds);
            c.cell_id = static_cast<std::uint32_t>(cells_.size());
            leaves_.push_back(c.cell_id);
            cells_.push_back(c);

            for (std::size_t d = kDim; d-- > 0;) {
                if (++idx[d] < static_cast<std::size_t>(resolution)) break;
                idx[d] = 0;
            }
        }
    }

    void build_spatial_index() {
        std::printf("  Building spatial index for %zu cells...\n", leaves_.size());
        const double wall0 = now_seconds();
        leaf_boxes_.resize(leaves_.size());
        for (std::size_t i = 0; i < leaves_.size(); ++i) leaf_boxes_[i] = cells_[leaves_[i]].bounds;
        index_.build(leaf_boxes_);
        index_build_seconds_ = now_seconds() - wall0;
        std::printf("   Spatial index built in %.2fs\n", index_build_seconds_);
    }

    Bounds3 root_bounds_;
    std::vector<Cell> cells_;
    std::vector<std::uint32_t> leaves_;
    std::vector<Bounds3> leaf_boxes_;
    BoxIndex index_;
    double index_build_seconds_ = 0.0;
};

}  // namespace dhj
