// checkpoint.hpp - binary checkpoint/resume for the adaptive refinement loop.
//
// A checkpoint holds exactly what is needed to restart Algorithm 2 at the top of
// a refinement phase: the cell arena (geometry, tree links and the six interval
// values per cell), the leaf list in order, and the two loop counters. Every
// other piece of solver state is either derivable or cheap to rebuild:
//
//   * the BVH             rebuilt by CellTree's restore constructor
//   * cell centers        recomputed from bounds
//   * the successor cache NOT stored - see below
//
// The successor cache is deliberately left out. It is keyed by (cell_id, action)
// and at the cell counts that make checkpointing worth doing it runs to several
// GB, which would dominate both the file size and the time to write it. Leaving
// it empty on resume is *exactly* correct rather than approximate: the cache is
// a pure function of the leaf set and the dynamics, and AdaptiveRefinement's
// first call to local_value_iteration re-derives it in full (with an empty cache
// every leaf reports a missing action, so `affected` covers the whole leaf set).
// Resuming therefore costs one full successor recompute, which is far less than
// replaying the phases that led here.
//
// ---------------------------------------------------------------------------
// On-disk format, little-endian, shared byte-for-byte with the Rust port so a
// checkpoint written by either implementation can be resumed by the other.
//
//   magic              char[8]  "DHJCKPT1"
//   version            u32      = 1
//   dim                u32      must equal kDim of the running build
//   mode               u32      must equal the running Mode
//   phase              u32      number of refinement phases already completed
//   total_refined      u64      cumulative parent cells split so far
//   initial_resolution u32
//   _pad               u32      = 0
//   epsilon            f64
//   root_bounds        f64[2*dim]
//   num_cells          u64
//   num_leaves         u64
//   cells              record[num_cells]
//   leaves             u32[num_leaves]
//
// One cell record (16*dim + 56 bytes; 120 bytes at dim 4):
//   bounds   f64[2*dim]   lo,hi interleaved per dimension
//   parent   i32
//   child0   i32
//   values   f64[6]       V_upper, V_lower, l_upper, l_lower, r_upper, r_lower
//
// cell_id, child1 and is_leaf are derived on load: cell_id is the array index,
// child1 == child0 + 1 because refine_cells always appends the pair together,
// and is_leaf == (child0 < 0).
#pragma once

#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include "args.hpp"
#include "cell.hpp"
#include "dim.hpp"

namespace dhj {

constexpr char kCkptMagic[8] = {'D', 'H', 'J', 'C', 'K', 'P', 'T', '1'};
constexpr std::uint32_t kCkptVersion = 1;

// Bytes per serialized cell record.
constexpr std::size_t kCellRecordBytes = 16 * kDim + 56;

// Loop state that lives in AdaptiveRefinement rather than in the tree.
struct CheckpointMeta {
    std::uint32_t phase = 0;
    std::uint64_t total_refined = 0;
    std::uint32_t initial_resolution = 0;
    double epsilon = 0.0;
    Bounds3 root_bounds{};
};

// --------------------------------------------------------------------- write --

namespace detail {

inline void put_u32(unsigned char*& p, std::uint32_t v) { std::memcpy(p, &v, 4); p += 4; }
inline void put_i32(unsigned char*& p, std::int32_t v) { std::memcpy(p, &v, 4); p += 4; }
inline void put_u64(unsigned char*& p, std::uint64_t v) { std::memcpy(p, &v, 8); p += 8; }
inline void put_f64(unsigned char*& p, double v) { std::memcpy(p, &v, 8); p += 8; }

inline std::uint32_t get_u32(const unsigned char*& p) { std::uint32_t v; std::memcpy(&v, p, 4); p += 4; return v; }
inline std::int32_t get_i32(const unsigned char*& p) { std::int32_t v; std::memcpy(&v, p, 4); p += 4; return v; }
inline std::uint64_t get_u64(const unsigned char*& p) { std::uint64_t v; std::memcpy(&v, p, 8); p += 8; return v; }
inline double get_f64(const unsigned char*& p) { double v; std::memcpy(&v, p, 8); p += 8; return v; }

inline void must_write(std::FILE* f, const void* p, std::size_t n, const std::string& path) {
    if (std::fwrite(p, 1, n, f) != n)
        throw std::runtime_error("checkpoint: short write to " + path);
}

inline void must_read(std::FILE* f, void* p, std::size_t n, const std::string& path) {
    if (std::fread(p, 1, n, f) != n)
        throw std::runtime_error("checkpoint: truncated file " + path);
}

}  // namespace detail

// Writes a checkpoint for `tree` to `path`. The file is built under a .tmp
// sibling and renamed into place, so a kill part-way through leaves the previous
// checkpoint intact rather than a half-written one.
inline void save_checkpoint(const std::string& path, Mode mode, const CellTree& tree,
                            const CheckpointMeta& meta) {
    const std::string tmp = path + ".tmp";
    std::error_code ec;
    std::filesystem::create_directories(std::filesystem::path(path).parent_path(), ec);

    std::FILE* f = std::fopen(tmp.c_str(), "wb");
    if (!f) throw std::runtime_error("checkpoint: cannot open " + tmp + " for writing");

    try {
        unsigned char head[8 + 4 * 4 + 8 + 4 + 4 + 8 + 16 * kDim + 16];
        unsigned char* p = head;
        std::memcpy(p, kCkptMagic, 8); p += 8;
        detail::put_u32(p, kCkptVersion);
        detail::put_u32(p, static_cast<std::uint32_t>(kDim));
        detail::put_u32(p, static_cast<std::uint32_t>(mode));
        detail::put_u32(p, meta.phase);
        detail::put_u64(p, meta.total_refined);
        detail::put_u32(p, meta.initial_resolution);
        detail::put_u32(p, 0);
        detail::put_f64(p, meta.epsilon);
        for (std::size_t d = 0; d < kDim; ++d) {
            detail::put_f64(p, meta.root_bounds[d][0]);
            detail::put_f64(p, meta.root_bounds[d][1]);
        }
        detail::put_u64(p, static_cast<std::uint64_t>(tree.num_cells()));
        detail::put_u64(p, static_cast<std::uint64_t>(tree.num_leaves()));
        detail::must_write(f, head, static_cast<std::size_t>(p - head), tmp);

        // Cells, in blocks so the whole arena is never duplicated in memory.
        constexpr std::size_t kBlock = 65536;
        std::vector<unsigned char> buf(kBlock * kCellRecordBytes);
        std::size_t written = 0;
        while (written < tree.num_cells()) {
            const std::size_t n = std::min(kBlock, tree.num_cells() - written);
            unsigned char* q = buf.data();
            for (std::size_t i = 0; i < n; ++i) {
                const Cell& c = tree.cell(static_cast<std::uint32_t>(written + i));
                for (std::size_t d = 0; d < kDim; ++d) {
                    detail::put_f64(q, c.bounds[d][0]);
                    detail::put_f64(q, c.bounds[d][1]);
                }
                detail::put_i32(q, static_cast<std::int32_t>(c.parent));
                detail::put_i32(q, static_cast<std::int32_t>(c.child0));
                detail::put_f64(q, c.V_upper);
                detail::put_f64(q, c.V_lower);
                detail::put_f64(q, c.l_upper);
                detail::put_f64(q, c.l_lower);
                detail::put_f64(q, c.r_upper);
                detail::put_f64(q, c.r_lower);
            }
            detail::must_write(f, buf.data(), n * kCellRecordBytes, tmp);
            written += n;
        }

        detail::must_write(f, tree.leaves().data(), tree.num_leaves() * sizeof(std::uint32_t), tmp);
    } catch (...) {
        std::fclose(f);
        std::filesystem::remove(tmp, ec);
        throw;
    }

    if (std::fclose(f) != 0) {
        std::filesystem::remove(tmp, ec);
        throw std::runtime_error("checkpoint: failed to close " + tmp);
    }
    std::filesystem::rename(tmp, path, ec);
    if (ec) throw std::runtime_error("checkpoint: cannot rename " + tmp + " -> " + path);
}

// ---------------------------------------------------------------------- read --

// Loads `path`, returning the restored tree by out-parameter and the loop state
// as the result. Throws if the file was written by a build with a different
// dimension or mode, because silently continuing would produce wrong answers.
inline CheckpointMeta load_checkpoint(const std::string& path, Mode mode,
                                      std::unique_ptr<CellTree>& tree_out) {
    std::FILE* f = std::fopen(path.c_str(), "rb");
    if (!f) throw std::runtime_error("checkpoint: cannot open " + path + " for reading");

    CheckpointMeta meta;
    std::size_t num_cells = 0, num_leaves = 0;
    std::vector<Cell> cells;
    std::vector<std::uint32_t> leaves;

    try {
        unsigned char head[8 + 4 * 4 + 8 + 4 + 4 + 8 + 16 * kDim + 16];
        const std::size_t head_bytes = 8 + 4 + 4 + 4 + 4 + 8 + 4 + 4 + 8 + 16 * kDim + 8 + 8;
        detail::must_read(f, head, head_bytes, path);
        const unsigned char* p = head;

        if (std::memcmp(p, kCkptMagic, 8) != 0)
            throw std::runtime_error("checkpoint: " + path + " is not a DHJ checkpoint");
        p += 8;

        const std::uint32_t version = detail::get_u32(p);
        if (version != kCkptVersion)
            throw std::runtime_error("checkpoint: " + path + " has version " +
                                     std::to_string(version) + ", this build reads version " +
                                     std::to_string(kCkptVersion));

        const std::uint32_t dim = detail::get_u32(p);
        if (dim != kDim)
            throw std::runtime_error("checkpoint: " + path + " is " + std::to_string(dim) +
                                     "-D but this build is " + std::to_string(kDim) +
                                     "-D (rebuild with the matching DHJ_DIM)");

        const std::uint32_t file_mode = detail::get_u32(p);
        if (file_mode != static_cast<std::uint32_t>(mode))
            throw std::runtime_error("checkpoint: " + path + " was written by a different "
                                     "solver mode; resume it with the binary that produced it");

        meta.phase = detail::get_u32(p);
        meta.total_refined = detail::get_u64(p);
        meta.initial_resolution = detail::get_u32(p);
        (void)detail::get_u32(p);  // pad
        meta.epsilon = detail::get_f64(p);
        for (std::size_t d = 0; d < kDim; ++d) {
            meta.root_bounds[d][0] = detail::get_f64(p);
            meta.root_bounds[d][1] = detail::get_f64(p);
        }
        num_cells = static_cast<std::size_t>(detail::get_u64(p));
        num_leaves = static_cast<std::size_t>(detail::get_u64(p));

        cells.resize(num_cells);
        constexpr std::size_t kBlock = 65536;
        std::vector<unsigned char> buf(kBlock * kCellRecordBytes);
        std::size_t done = 0;
        while (done < num_cells) {
            const std::size_t n = std::min(kBlock, num_cells - done);
            detail::must_read(f, buf.data(), n * kCellRecordBytes, path);
            const unsigned char* q = buf.data();
            for (std::size_t i = 0; i < n; ++i) {
                Cell& c = cells[done + i];
                for (std::size_t d = 0; d < kDim; ++d) {
                    c.bounds[d][0] = detail::get_f64(q);
                    c.bounds[d][1] = detail::get_f64(q);
                }
                const std::int32_t parent = detail::get_i32(q);
                const std::int32_t child0 = detail::get_i32(q);
                c.V_upper = detail::get_f64(q);
                c.V_lower = detail::get_f64(q);
                c.l_upper = detail::get_f64(q);
                c.l_lower = detail::get_f64(q);
                c.r_upper = detail::get_f64(q);
                c.r_lower = detail::get_f64(q);

                c.center = center_of(c.bounds);
                c.cell_id = static_cast<std::uint32_t>(done + i);
                c.parent = parent;
                c.child0 = child0;
                c.child1 = child0 < 0 ? -1 : child0 + 1;
                c.is_leaf = child0 < 0;
            }
            done += n;
        }

        leaves.resize(num_leaves);
        detail::must_read(f, leaves.data(), num_leaves * sizeof(std::uint32_t), path);
    } catch (...) {
        std::fclose(f);
        throw;
    }
    std::fclose(f);

    for (std::uint32_t id : leaves)
        if (id >= num_cells || !cells[id].is_leaf)
            throw std::runtime_error("checkpoint: " + path + " has an inconsistent leaf list");

    // parent/child are raw file data and are used as arena subscripts, so range
    // check them here rather than trusting the file. -1 means "none".
    const std::int64_t n = static_cast<std::int64_t>(num_cells);
    for (std::size_t i = 0; i < num_cells; ++i) {
        const Cell& c = cells[i];
        const bool bad_parent = c.parent < -1 || c.parent >= n;
        const bool bad_child  = c.child0 < -1 || c.child0 >= n ||
                                c.child1 < -1 || c.child1 >= n;
        if (bad_parent || bad_child)
            throw std::runtime_error("checkpoint: " + path + " has an out-of-range "
                                     "parent/child index at cell " + std::to_string(i));
    }

    tree_out.reset(new CellTree(meta.root_bounds, std::move(cells), std::move(leaves)));
    return meta;
}

// ------------------------------------------------------------------- helpers --

inline std::string checkpoint_path(const std::string& dir, int phase) {
    char buf[64];
    std::snprintf(buf, sizeof buf, "/checkpoint_phase_%04d.bin", phase);
    return dir + buf;
}

// Highest-numbered checkpoint in `dir`, or an empty string if there is none.
// Lets `--resume auto` pick up wherever the last run stopped.
inline std::string latest_checkpoint(const std::string& dir) {
    std::error_code ec;
    if (!std::filesystem::is_directory(dir, ec)) return {};
    std::string best;
    int best_phase = -1;
    constexpr const char* kPrefix = "checkpoint_phase_";
    constexpr const char* kSuffix = ".bin";
    const std::size_t plen = std::strlen(kPrefix), slen = std::strlen(kSuffix);

    for (const auto& e : std::filesystem::directory_iterator(dir, ec)) {
        const std::string name = e.path().filename().string();
        // Anchored at both ends on purpose: sscanf would happily accept the
        // "...bin.tmp" file a concurrent save is still writing, and resuming
        // from a half-written checkpoint is exactly the failure this is meant
        // to prevent.
        if (name.size() <= plen + slen) continue;
        if (name.compare(0, plen, kPrefix) != 0) continue;
        if (name.compare(name.size() - slen, slen, kSuffix) != 0) continue;

        const std::string digits = name.substr(plen, name.size() - plen - slen);
        if (digits.empty() || digits.find_first_not_of("0123456789") != std::string::npos) continue;

        const int phase = std::atoi(digits.c_str());
        if (phase > best_phase) {
            best_phase = phase;
            best = e.path().string();
        }
    }
    return best;
}

}  // namespace dhj
