// plot.hpp - SVG rendering of the value-function slices.
//
// The Python draws with matplotlib; there is no equivalent to link against
// here, so the same figure (3 rows x 6 theta slices of coloured leaf cells,
// plus obstacle/target circles and a colourbar) is emitted as a standalone SVG.
// Filenames keep the base name from the Python and swap the .png extension.
#pragma once

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <string>
#include <vector>

#include "args.hpp"
#include "cell.hpp"
#include "env.hpp"
#include "image.hpp"

namespace dhj {

#if DHJ_DIM == 3


struct PlotStyle {
    bool draw_value_cells = true;  // avoid_nodiscount.py leaves these panels blank
    bool cell_labels = false;      // RA_discount.py prints the value inside big cells
    bool draw_target = true;
};

// matplotlib's RdYlGn, an 11-anchor ColorBrewer map with linear interpolation.
inline void rd_yl_gn(double t, int rgb[3]) {
    static const int anchors[11][3] = {
        {0xa5, 0x00, 0x26}, {0xd7, 0x30, 0x27}, {0xf4, 0x6d, 0x43}, {0xfd, 0xae, 0x61},
        {0xfe, 0xe0, 0x8b}, {0xff, 0xff, 0xbf}, {0xd9, 0xef, 0x8b}, {0xa6, 0xd9, 0x6a},
        {0x66, 0xbd, 0x63}, {0x1a, 0x98, 0x50}, {0x00, 0x68, 0x37}};
    t = std::min(1.0, std::max(0.0, t));
    const double x = t * 10.0;
    int i = static_cast<int>(x);
    if (i > 9) i = 9;
    const double f = x - i;
    for (int c = 0; c < 3; ++c)
        rgb[c] = static_cast<int>(std::lround(anchors[i][c] + f * (anchors[i + 1][c] - anchors[i][c])));
}

inline std::string hex_color(int rgb[3]) {
    char buf[16];
    std::snprintf(buf, sizeof buf, "#%02x%02x%02x", rgb[0], rgb[1], rgb[2]);
    return buf;
}

inline std::string with_svg_extension(const std::string& path) {
    const std::size_t dot = path.rfind(".png");
    if (dot != std::string::npos && dot + 4 == path.size()) return path.substr(0, dot) + ".svg";
    return path + ".svg";
}

class ValueFunctionPlotter {
public:
    ValueFunctionPlotter(const Environment& env, PlotStyle style) : env_(env), style_(style) {}

    void plot(const CellTree& tree, const std::string& png_filename, int iteration) const {
        static const double thetas[6] = {0.0, kPi, kPi / 4, -kPi / 4, kPi / 2, -kPi / 2};
        constexpr int kPanel = 460;   // drawing area, px
        constexpr int kPadL = 60, kPadR = 96, kPadT = 44, kPadB = 52;
        const int cell_w = kPanel + kPadL + kPadR;
        const int cell_h = kPanel + kPadT + kPadB;
        const int width = cell_w * 6;
        const int height = cell_h * 3 + 40;

        std::ofstream out(with_svg_extension(png_filename));
        if (!out) {
            std::fprintf(stderr, "  [plot] could not open %s\n", with_svg_extension(png_filename).c_str());
            return;
        }
        out << "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"" << width << "\" height=\"" << height
            << "\" viewBox=\"0 0 " << width << " " << height << "\">\n"
            << "<rect width=\"100%\" height=\"100%\" fill=\"white\"/>\n"
            << "<text x=\"" << width / 2 << "\" y=\"26\" font-family=\"sans-serif\" font-size=\"22\""
               " text-anchor=\"middle\">Safety Value Function - Iteration " << iteration << "</text>\n";

        for (int row = 0; row < 3; ++row)
            for (int col = 0; col < 6; ++col) {
                const int ox = col * cell_w + kPadL;
                const int oy = 40 + row * cell_h + kPadT;
                out << "<g>\n";
                panel(out, tree, thetas[col], row, ox, oy, kPanel);
                out << "</g>\n";
            }
        out << "</svg>\n";
    }

    // Rasterizes the six theta slices onto an n x n grid and writes them as raw
    // f64 (V_lower then V_upper per slice), with a sidecar text file describing
    // the geometry. At the cell counts adaptive refinement reaches, per-cell
    // vector output is impractical; this stays small and renders directly.
    void dump_slices(const CellTree& tree, const std::string& stem, int n) const {
        static const double thetas[6] = {0.0, kPi, kPi / 4, -kPi / 4, kPi / 2, -kPi / 2};
        const Bounds3& gb = env_.state_bounds();
        const double x_min = gb[0][0], x_max = gb[0][1];
        const double y_min = gb[1][0], y_max = gb[1][1];

        std::ofstream bin(stem + ".bin", std::ios::binary);
        if (!bin) {
            std::fprintf(stderr, "  [slices] could not open %s.bin\n", stem.c_str());
            return;
        }

        std::vector<double> lower(static_cast<std::size_t>(n) * n);
        std::vector<double> upper(static_cast<std::size_t>(n) * n);
        std::vector<std::vector<double>> lower_all, upper_all;
        for (double theta : thetas) {
            for (int iy = 0; iy < n; ++iy) {
                // Pixel centres, so a query never lands exactly on a grid face.
                const double y = y_min + (y_max - y_min) * (iy + 0.5) / n;
                for (int ix = 0; ix < n; ++ix) {
                    const double x = x_min + (x_max - x_min) * (ix + 0.5) / n;
                    const std::size_t k = static_cast<std::size_t>(iy) * n + ix;
                    const std::int64_t id = tree.find_leaf_containing(State3{x, y, theta});
                    if (id < 0) {
                        lower[k] = upper[k] = std::nan("");
                    } else {
                        const Cell& c = tree.cell(static_cast<std::uint32_t>(id));
                        lower[k] = c.V_lower;
                        upper[k] = c.V_upper;
                    }
                }
            }
            bin.write(reinterpret_cast<const char*>(lower.data()),
                      static_cast<std::streamsize>(lower.size() * sizeof(double)));
            bin.write(reinterpret_cast<const char*>(upper.data()),
                      static_cast<std::streamsize>(upper.size() * sizeof(double)));
            lower_all.push_back(lower);
            upper_all.push_back(upper);
        }

        render_slice_png(stem + ".png", n, thetas, lower_all, upper_all);

        std::ofstream meta(stem + ".meta");
        meta.precision(17);
        meta << "n " << n << "\n"
             << "x_min " << x_min << "\nx_max " << x_max << "\n"
             << "y_min " << y_min << "\ny_max " << y_max << "\n"
             << "n_slices 6\nthetas";
        for (double t : thetas) meta << ' ' << t;
        meta << "\nleaves " << tree.num_leaves() << "\n";

        std::array<double, 2> obs{}, tgt{};
        double obs_r = 0.0, tgt_r = 0.0;
        if (const auto* d = dynamic_cast<const DubinsCarEnvironment*>(&env_)) {
            obs = d->obstacle_position(); obs_r = d->obstacle_radius();
            tgt = d->target_position();   tgt_r = d->target_radius();
        } else if (const auto* e = dynamic_cast<const EvasionEnvironment*>(&env_)) {
            obs = e->obstacle_position(); obs_r = e->obstacle_radius();
            tgt = e->target_position();   tgt_r = e->target_radius();
        }
        meta << "obstacle " << obs[0] << ' ' << obs[1] << ' ' << obs_r << "\n"
             << "target " << tgt[0] << ' ' << tgt[1] << ' ' << tgt_r << "\n"
             << "draw_target " << (style_.draw_target ? 1 : 0) << "\n";
    }

    // Composites the six slices into one figure: row 0 upper bound, row 1 lower
    // bound, row 2 classification (green safe / red unsafe / grey undetermined).
    // Same layout and colours as the matplotlib figure in the Python.
    void render_slice_png(const std::string& path, int n, const double* thetas,
                          const std::vector<std::vector<double>>& lower_all,
                          const std::vector<std::vector<double>>& upper_all) const {
        const int cols = static_cast<int>(lower_all.size());
        const int gap = 8;
        Canvas cv(cols * n + (cols + 1) * gap, 3 * n + 4 * gap);
        cv.rect(0, 0, cv.w, cv.h, 245, 245, 245);

        const Bounds3& gb = env_.state_bounds();
        const double x_min = gb[0][0], x_max = gb[0][1];
        const double y_min = gb[1][0], y_max = gb[1][1];

        std::array<double, 2> obs{}, tgt{};
        double obs_r = 0.0, tgt_r = 0.0;
        if (const auto* d = dynamic_cast<const DubinsCarEnvironment*>(&env_)) {
            obs = d->obstacle_position(); obs_r = d->obstacle_radius();
            tgt = d->target_position();   tgt_r = d->target_radius();
        } else if (const auto* e = dynamic_cast<const EvasionEnvironment*>(&env_)) {
            obs = e->obstacle_position(); obs_r = e->obstacle_radius();
            tgt = e->target_position();   tgt_r = e->target_radius();
        }

        for (int col = 0; col < cols; ++col) {
            const std::vector<double>& lo = lower_all[static_cast<std::size_t>(col)];
            const std::vector<double>& up = upper_all[static_cast<std::size_t>(col)];
            (void)thetas;

            for (int row = 0; row < 3; ++row) {
                const int ox = gap + col * (n + gap);
                const int oy = gap + row * (n + gap);

                double vmin = 0.0, vmax = 0.0;
                if (row < 2) {
                    const std::vector<double>& v = (row == 0) ? up : lo;
                    vmin = 1e300; vmax = -1e300;
                    for (double t : v) {
                        if (!std::isfinite(t)) continue;
                        vmin = std::min(vmin, t);
                        vmax = std::max(vmax, t);
                    }
                    if (vmax <= 0.0) vmax = 0.0; else if (vmin >= 0.0) vmin = 0.0;
                    if (!(vmax > vmin)) { vmin -= 1e-12; vmax += 1e-12; }
                }

                for (int iy = 0; iy < n; ++iy) {
                    // Slice rows run bottom-up; image rows run top-down.
                    const int src_y = n - 1 - iy;
                    for (int ix = 0; ix < n; ++ix) {
                        const std::size_t k = static_cast<std::size_t>(src_y) * n + ix;
                        int rgb[3];
                        if (row < 2) {
                            const double value = (row == 0) ? up[k] : lo[k];
                            if (!std::isfinite(value)) { rgb[0] = rgb[1] = rgb[2] = 255; }
                            else rd_yl_gn((value - vmin) / (vmax - vmin), rgb);
                        } else if (lo[k] > 0.0) {
                            rgb[0] = 0x2c; rgb[1] = 0xa0; rgb[2] = 0x2c;   // safe
                        } else if (up[k] <= 0.0) {
                            rgb[0] = 0xd6; rgb[1] = 0x27; rgb[2] = 0x28;   // unsafe
                        } else {
                            rgb[0] = rgb[1] = rgb[2] = 0x7f;               // undetermined
                        }
                        cv.set(ox + ix, oy + iy, rgb[0], rgb[1], rgb[2]);
                    }
                }

                const double sx = n / (x_max - x_min);
                const double sy = n / (y_max - y_min);
                const double cx = ox + (obs[0] - x_min) * sx;
                const double cy = oy + n - (obs[1] - y_min) * sy;
                cv.circle(cx, cy, obs_r * sx, obs_r * sy, 0x00, 0x00, 0x8b, 1);
                if (style_.draw_target) {
                    const double tx = ox + (tgt[0] - x_min) * sx;
                    const double ty = oy + n - (tgt[1] - y_min) * sy;
                    cv.dashed_circle(tx, ty, tgt_r * sx, tgt_r * sy, 0xff, 0xa5, 0x00, 1);
                }
                cv.rect(ox - 1, oy - 1, ox + n + 1, oy, 0, 0, 0);
                cv.rect(ox - 1, oy + n, ox + n + 1, oy + n + 1, 0, 0, 0);
                cv.rect(ox - 1, oy, ox, oy + n, 0, 0, 0);
                cv.rect(ox + n, oy, ox + n + 1, oy + n, 0, 0, 0);
            }
        }

        if (!write_png(path, cv.w, cv.h, cv.px))
            std::fprintf(stderr, "  [slices] could not write %s\n", path.c_str());
    }

    // Plain CSV of every leaf: bounds and both value bounds.
    void dump_csv(const CellTree& tree, const std::string& path) const {
        std::ofstream out(path);
        if (!out) return;
        out << "cell_id,x_lo,x_hi,y_lo,y_hi,theta_lo,theta_hi,V_lower,V_upper,l_lower,l_upper,r_lower,r_upper\n";
        out.precision(17);
        for (std::uint32_t id : tree.leaves()) {
            const Cell& c = tree.cell(id);
            out << c.cell_id << ',' << c.bounds[0][0] << ',' << c.bounds[0][1] << ','
                << c.bounds[1][0] << ',' << c.bounds[1][1] << ',' << c.bounds[2][0] << ','
                << c.bounds[2][1] << ',' << c.V_lower << ',' << c.V_upper << ',' << c.l_lower
                << ',' << c.l_upper << ',' << c.r_lower << ',' << c.r_upper << '\n';
        }
    }

private:
    // row 0: upper bound, row 1: lower bound, row 2: classification.
    void panel(std::ofstream& out, const CellTree& tree, double theta, int row, int ox, int oy,
               int size) const {
        const Bounds3& gb = env_.state_bounds();
        const double x_min = gb[0][0], x_max = gb[0][1];
        const double y_min = gb[1][0], y_max = gb[1][1];
        const double sx = size / (x_max - x_min);
        const double sy = size / (y_max - y_min);
        auto px = [&](double x) { return ox + (x - x_min) * sx; };
        auto py = [&](double y) { return oy + size - (y - y_min) * sy; };

        std::vector<std::uint32_t> hits;
        std::vector<double> values;
        for (std::uint32_t id : tree.leaves()) {
            const Cell& c = tree.cell(id);
            if (!(c.bounds[2][0] <= theta && theta <= c.bounds[2][1])) continue;
            hits.push_back(id);
            if (row < 2) values.push_back(row == 0 ? c.V_upper : c.V_lower);
        }

        const char* titles[3] = {"Upper Bound V̄_g", "Lower Bound V_g", "Cell Classification"};
        char header[128];
        std::snprintf(header, sizeof header, "%s (theta=%.2f rad)", titles[row], theta);
        out << "<text x=\"" << ox + size / 2 << "\" y=\"" << oy - 12
            << "\" font-family=\"sans-serif\" font-size=\"13\" text-anchor=\"middle\">" << header
            << "</text>\n";

        if (row < 2 && values.empty()) {
            out << "<text x=\"" << ox + size / 2 << "\" y=\"" << oy + size / 2
                << "\" font-family=\"sans-serif\" font-size=\"13\" text-anchor=\"middle\">"
                   "No data for this slice</text>\n";
            frame(out, ox, oy, size, x_min, x_max, y_min, y_max);
            return;
        }

        double vmin = 0.0, vmax = 0.0;
        if (row < 2) {
            vmin = *std::min_element(values.begin(), values.end());
            vmax = *std::max_element(values.begin(), values.end());
            if (std::fabs(vmax - vmin) <= 1e-8 * std::max(std::fabs(vmin), std::fabs(vmax))) {
                vmin -= 1e-12;
                vmax += 1e-12;
            }
            if (vmax <= 0.0) vmax = 0.0;
            else if (vmin >= 0.0) vmin = 0.0;
        }

        const bool draw_cells = (row == 2) || style_.draw_value_cells;
        if (draw_cells) {
            std::size_t vi = 0;
            for (std::uint32_t id : hits) {
                const Cell& c = tree.cell(id);
                std::string fill;
                double value = 0.0, tnorm = 0.0;
                if (row < 2) {
                    value = values[vi++];
                    tnorm = (vmax > vmin) ? (value - vmin) / (vmax - vmin) : 0.5;
                    int rgb[3];
                    rd_yl_gn(tnorm, rgb);
                    fill = hex_color(rgb);
                } else {
                    fill = c.V_lower > 0.0 ? "#2ca02c" : (c.V_upper <= 0.0 ? "#d62728" : "#7f7f7f");
                }
                const double x0 = px(c.bounds[0][0]), x1 = px(c.bounds[0][1]);
                const double y0 = py(c.bounds[1][1]), y1 = py(c.bounds[1][0]);
                out << "<rect x=\"" << x0 << "\" y=\"" << y0 << "\" width=\"" << (x1 - x0)
                    << "\" height=\"" << (y1 - y0) << "\" fill=\"" << fill
                    << "\" stroke=\"black\" stroke-width=\"0.05\" shape-rendering=\"crispEdges\"/>\n";

                if (row < 2 && style_.cell_labels) {
                    const double w = c.bounds[0][1] - c.bounds[0][0];
                    const double h = c.bounds[1][1] - c.bounds[1][0];
                    if (w > 0.05 && h > 0.05) {
                        const double fs = std::min(std::min(w * 15.0, h * 15.0), 6.0);
                        char lbl[32];
                        std::snprintf(lbl, sizeof lbl, "%.1e", value);
                        out << "<text x=\"" << 0.5 * (x0 + x1) << "\" y=\"" << 0.5 * (y0 + y1)
                            << "\" font-family=\"sans-serif\" font-size=\"" << fs
                            << "\" text-anchor=\"middle\" dominant-baseline=\"middle\" fill=\""
                            << (tnorm > 0.5 ? "black" : "white") << "\">" << lbl << "</text>\n";
                    }
                }
            }
        }

        circles(out, px, py, sx, sy);
        frame(out, ox, oy, size, x_min, x_max, y_min, y_max);
        if (row < 2) colorbar(out, ox + size + 14, oy, size, vmin, vmax);
    }

    template <typename PxFn, typename PyFn>
    void circles(std::ofstream& out, PxFn px, PyFn py, double sx, double sy) const {
        std::array<double, 2> obs{}, tgt{};
        double obs_r = 0.0, tgt_r = 0.0;
        if (const auto* d = dynamic_cast<const DubinsCarEnvironment*>(&env_)) {
            obs = d->obstacle_position(); obs_r = d->obstacle_radius();
            tgt = d->target_position();   tgt_r = d->target_radius();
        } else if (const auto* e = dynamic_cast<const EvasionEnvironment*>(&env_)) {
            obs = e->obstacle_position(); obs_r = e->obstacle_radius();
            tgt = e->target_position();   tgt_r = e->target_radius();
        } else {
            return;
        }
        out << "<ellipse cx=\"" << px(obs[0]) << "\" cy=\"" << py(obs[1]) << "\" rx=\"" << obs_r * sx
            << "\" ry=\"" << obs_r * sy << "\" fill=\"none\" stroke=\"darkblue\" stroke-width=\"2\"/>\n";
        if (style_.draw_target)
            out << "<ellipse cx=\"" << px(tgt[0]) << "\" cy=\"" << py(tgt[1]) << "\" rx=\"" << tgt_r * sx
                << "\" ry=\"" << tgt_r * sy
                << "\" fill=\"none\" stroke=\"orange\" stroke-width=\"2\" stroke-dasharray=\"6 4\"/>\n";
    }

    void frame(std::ofstream& out, int ox, int oy, int size, double x_min, double x_max, double y_min,
               double y_max) const {
        out << "<rect x=\"" << ox << "\" y=\"" << oy << "\" width=\"" << size << "\" height=\"" << size
            << "\" fill=\"none\" stroke=\"black\" stroke-width=\"1\"/>\n";
        char b[64];
        std::snprintf(b, sizeof b, "%.1f", x_min);
        out << "<text x=\"" << ox << "\" y=\"" << oy + size + 16
            << "\" font-family=\"sans-serif\" font-size=\"11\" text-anchor=\"middle\">" << b << "</text>\n";
        std::snprintf(b, sizeof b, "%.1f", x_max);
        out << "<text x=\"" << ox + size << "\" y=\"" << oy + size + 16
            << "\" font-family=\"sans-serif\" font-size=\"11\" text-anchor=\"middle\">" << b << "</text>\n";
        out << "<text x=\"" << ox + size / 2 << "\" y=\"" << oy + size + 34
            << "\" font-family=\"sans-serif\" font-size=\"13\" text-anchor=\"middle\">x</text>\n";
        std::snprintf(b, sizeof b, "%.1f", y_min);
        out << "<text x=\"" << ox - 6 << "\" y=\"" << oy + size
            << "\" font-family=\"sans-serif\" font-size=\"11\" text-anchor=\"end\">" << b << "</text>\n";
        std::snprintf(b, sizeof b, "%.1f", y_max);
        out << "<text x=\"" << ox - 6 << "\" y=\"" << oy + 10
            << "\" font-family=\"sans-serif\" font-size=\"11\" text-anchor=\"end\">" << b << "</text>\n";
        out << "<text x=\"" << ox - 40 << "\" y=\"" << oy + size / 2
            << "\" font-family=\"sans-serif\" font-size=\"13\" text-anchor=\"middle\">y</text>\n";
    }

    void colorbar(std::ofstream& out, int x, int y, int size, double vmin, double vmax) const {
        constexpr int kSteps = 64;
        const int w = 14;
        for (int i = 0; i < kSteps; ++i) {
            const double t = 1.0 - static_cast<double>(i) / (kSteps - 1);
            int rgb[3];
            rd_yl_gn(t, rgb);
            out << "<rect x=\"" << x << "\" y=\"" << y + i * size / kSteps << "\" width=\"" << w
                << "\" height=\"" << (size / kSteps + 1) << "\" fill=\"" << hex_color(rgb)
                << "\" shape-rendering=\"crispEdges\"/>\n";
        }
        out << "<rect x=\"" << x << "\" y=\"" << y << "\" width=\"" << w << "\" height=\"" << size
            << "\" fill=\"none\" stroke=\"black\" stroke-width=\"0.7\"/>\n";
        char b[32];
        std::snprintf(b, sizeof b, "%.1e", vmax);
        out << "<text x=\"" << x + w + 4 << "\" y=\"" << y + 9
            << "\" font-family=\"sans-serif\" font-size=\"10\">" << b << "</text>\n";
        std::snprintf(b, sizeof b, "%.1e", vmin);
        out << "<text x=\"" << x + w + 4 << "\" y=\"" << y + size
            << "\" font-family=\"sans-serif\" font-size=\"10\">" << b << "</text>\n";
        if (std::fabs(vmin) >= 0.1 && std::fabs(vmax) >= 0.1 && vmin < 0.0 && vmax > 0.0) {
            const double t = (0.0 - vmin) / (vmax - vmin);
            const double yy = y + size - t * size;
            out << "<text x=\"" << x + w + 4 << "\" y=\"" << yy
                << "\" font-family=\"sans-serif\" font-size=\"10\">0.0e+00</text>\n";
        }
    }

    const Environment& env_;
    PlotStyle style_;
};

#else  // DHJ_DIM != 3

// The plotter rasterizes (x, y) images at fixed theta, which only determines a
// slice when theta is the LAST dimension. At kDim==4 every v-slab would land on
// the same pixel and overdraw, so plotting is disabled here rather than emitting
// a silently misleading picture. Use the offline slicer instead:
//   Scratchpad/4DNewDynamics/plot_slices_4d.py   (pins theta AND v)
struct PlotStyle {
    bool draw_value_cells = true;
    bool cell_labels = false;
    bool draw_target = true;
};

class ValueFunctionPlotter {
public:
    ValueFunctionPlotter(const Environment&, PlotStyle) {}

    // Rasterizing is still disabled at 4D (see the note above), but the cell
    // data is dumped so the offline slicer can pin theta AND v and draw a
    // genuine slice. Columns match the 3D dump plus v_lo/v_hi.
    void plot(const CellTree&, const std::string&, int) const {}

    void dump_csv(const CellTree& tree, const std::string& path) const {
        std::ofstream out(path);
        if (!out) return;
        out << "cell_id,x_lo,x_hi,y_lo,y_hi,theta_lo,theta_hi,v_lo,v_hi,"
               "V_lower,V_upper,l_lower,l_upper,r_lower,r_upper\n";
        out.precision(17);
        for (std::uint32_t id : tree.leaves()) {
            const Cell& c = tree.cell(id);
            out << c.cell_id << ',' << c.bounds[0][0] << ',' << c.bounds[0][1] << ','
                << c.bounds[1][0] << ',' << c.bounds[1][1] << ',' << c.bounds[2][0] << ','
                << c.bounds[2][1] << ',' << c.bounds[3][0] << ',' << c.bounds[3][1] << ','
                << c.V_lower << ',' << c.V_upper << ',' << c.l_lower << ',' << c.l_upper
                << ',' << c.r_lower << ',' << c.r_upper << '\n';
        }
    }

    void dump_slices(const CellTree&, const std::string&, int) const {}
};

#endif  // DHJ_DIM == 3

}  // namespace dhj
