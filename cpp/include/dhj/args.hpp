// args.hpp - command-line interface matching the argparse setup in the Python.
#pragma once

#include <cstdio>
#include <cstdlib>
#include <stdexcept>
#include <string>
#include <thread>

#include "env.hpp"
#include "pyfmt.hpp"

namespace dhj {

// Which of the three Python scripts this binary reproduces.
enum class Mode {
    RANoDiscount,      // RA_nodiscount.py
    RADiscount,        // RA_discount.py
    AvoidNoDiscount,   // avoid_nodiscount.py
};

inline const char* mode_script(Mode m) {
    switch (m) {
        case Mode::RANoDiscount: return "RA_nodiscount";
        case Mode::RADiscount: return "RA_discount";
        default: return "avoid_nodiscount";
    }
}

inline const char* mode_results_root(Mode m) {
    switch (m) {
        case Mode::RANoDiscount: return "./results/RA_no_discount";
        case Mode::RADiscount: return "./results/RAdiscount";
        default: return "./results/avoid_nodiscount";
    }
}

struct Args {
    int algorithm = 0;                 // required: 1 or 2
    Dynamics dynamics = Dynamics::Dubins;
    std::string dynamics_name = "dubins";
    double velocity = 1.0;
    double dt = 0.1;
    double tau = 1.0;
    double gamma = 0.96;
    int resolution = 10;
    int iterations = 200;              // default overridden per mode
    double tolerance = 1e-13;
    int plot_freq = 1000;
    double epsilon = 0.1;
    int initial_resolution = 15;
    int refinements = 100;
    int vi_iterations = 20000;
    bool conservative = false;
    double delta_max = 1e-6;
    int workers = 0;                   // 0 => cpu_count - 1
    bool precompute = false;

    // Additions over the Python CLI, for benchmarking and for running without a
    // plotting dependency.
    std::string out_dir;               // override the derived results directory
    bool no_plot = false;              // skip SVG output entirely
    bool dump_csv = false;             // write the value function as CSV per phase
    bool seq_queries = false;          // run successor queries single-threaded, as the Python does
    bool draw_value_cells = false;     // see mode_default_draw_value_cells()
    int dump_slices = 0;               // rasterize theta slices at this resolution

    // Checkpoint / resume for the Algorithm 2 refinement loop.
    int checkpoint_every = 0;          // 0 disables; N writes one every N phases
    std::string checkpoint_dir;        // default: <out_dir>/checkpoints
    std::string resume;                // path to a .bin, or "auto" for the latest
    int keep_checkpoints = 2;          // prune all but the newest K (0 keeps all)
};

inline int default_iterations(Mode m) {
    switch (m) {
        case Mode::RANoDiscount: return 200;
        case Mode::RADiscount: return 2000;
        default: return 20000;
    }
}

// avoid_nodiscount.py has the value-rectangle loop of _plot_slice commented out,
// so its upper/lower panels come out empty. Reproduced here by default; pass
// --draw-value-cells to fill them in.
inline bool mode_default_draw_value_cells(Mode m) { return m != Mode::AvoidNoDiscount; }

// RA_discount.py additionally prints the numeric value inside each large cell.
inline bool mode_default_cell_labels(Mode m) { return m == Mode::RADiscount; }

// avoid_nodiscount.py draws no target circle (there is no target in an
// avoid-only problem).
inline bool mode_draws_target(Mode m) { return m != Mode::AvoidNoDiscount; }

[[noreturn]] inline void usage_and_exit(const char* prog, Mode m, int code) {
    std::fprintf(stderr,
        "usage: %s --algorithm {1,2} [options]   (port of %s.py)\n\n"
        "  --dynamics {dubins,evasion,dubins_speed}  dynamics model (default: dubins)\n""                                (dubins_speed is 4D; needs the DHJ_DIM=4 build)\n"
        "  --velocity FLOAT              constant velocity (default: 1.0)\n"
        "  --dt FLOAT                    checkpoint interval (default: 0.1)\n"
        "  --tau FLOAT                   control duration (default: 1.0)\n"
        "  --gamma FLOAT                 discount factor (default: 0.96)\n"
        "  --resolution INT              grid resolution, algorithm 1 (default: 10)\n"
        "  --iterations INT              max value iterations, algorithm 1 (default: %d)\n"
        "  --tolerance/--delta-min FLOAT convergence tolerance (default: 1e-13)\n"
        "  --plot-freq INT               plot every N iterations (default: 1000)\n"
        "  --epsilon FLOAT               refinement error tolerance (default: 0.1)\n"
        "  --initial-resolution INT      initial coarse grid (default: 15)\n"
        "  --refinements INT             max refinement phases (default: 100)\n"
        "  --vi-iterations INT           VI iterations per refinement (default: 20000)\n"
        "  --conservative                conservative stopping (Algorithm 3)\n"
        "  --delta-max FLOAT             delta_max for conservative stopping (default: 1e-6)\n"
        "  --workers INT                 worker threads (default: cores - 1)\n"
        "  --precompute                  precompute successor sets up front\n"
        "\n  extras not present in the Python CLI:\n"
        "  --out-dir PATH                override the results directory\n"
        "  --no-plot                     skip SVG rendering (timing runs)\n"
        "  --dump-csv                    write value_function_phase_N.csv per phase\n"
        "  --seq-queries                 single-threaded successor queries\n"
        "  --draw-value-cells            force-fill the value panels\n"
        "  --dump-slices N               rasterize the theta slices onto an NxN grid\n"
        "\n  checkpoint / resume (algorithm 2 only):\n"
        "  --checkpoint-every N          write a checkpoint every N refinement phases\n"
        "  --checkpoint-dir PATH         where to write them (default: <out-dir>/checkpoints)\n"
        "  --resume PATH|auto            restart from a checkpoint ('auto' takes the latest\n"
        "                                in --checkpoint-dir)\n"
        "  --keep-checkpoints K          keep only the newest K (default: 2; 0 keeps all)\n",
        prog, mode_script(m), default_iterations(m));
    std::exit(code);
}

inline double parse_double(const char* s, const char* flag) {
    char* end = nullptr;
    const double v = std::strtod(s, &end);
    if (end == s || (end && *end != '\0')) throw std::invalid_argument(std::string("bad float for ") + flag);
    return v;
}

inline long parse_long(const char* s, const char* flag) {
    char* end = nullptr;
    const long v = std::strtol(s, &end, 10);
    if (end == s || (end && *end != '\0')) throw std::invalid_argument(std::string("bad int for ") + flag);
    return v;
}

inline Args parse_args(int argc, char** argv, Mode mode) {
    Args a;
    a.iterations = default_iterations(mode);
    a.draw_value_cells = mode_default_draw_value_cells(mode);

    bool draw_value_cells_set = false;
    auto need = [&](int i) -> const char* {
        if (i + 1 >= argc) { std::fprintf(stderr, "error: %s needs a value\n", argv[i]); std::exit(2); }
        return argv[i + 1];
    };

    for (int i = 1; i < argc; ++i) {
        const std::string f = argv[i];
        if (f == "-h" || f == "--help") usage_and_exit(argv[0], mode, 0);
        else if (f == "--algorithm") { a.algorithm = static_cast<int>(parse_long(need(i), "--algorithm")); ++i; }
        else if (f == "--dynamics") {
            a.dynamics_name = need(i); ++i;
            if (a.dynamics_name == "dubins") a.dynamics = Dynamics::Dubins;
            else if (a.dynamics_name == "evasion") a.dynamics = Dynamics::Evasion;
            else if (a.dynamics_name == "dubins_speed") a.dynamics = Dynamics::DubinsSpeed;
            else { std::fprintf(stderr, "error: --dynamics must be dubins, evasion or dubins_speed\n"); std::exit(2); }
        }
        else if (f == "--velocity") { a.velocity = parse_double(need(i), "--velocity"); ++i; }
        else if (f == "--dt") { a.dt = parse_double(need(i), "--dt"); ++i; }
        else if (f == "--tau") { a.tau = parse_double(need(i), "--tau"); ++i; }
        else if (f == "--gamma") { a.gamma = parse_double(need(i), "--gamma"); ++i; }
        else if (f == "--resolution") { a.resolution = static_cast<int>(parse_long(need(i), "--resolution")); ++i; }
        else if (f == "--iterations") { a.iterations = static_cast<int>(parse_long(need(i), "--iterations")); ++i; }
        else if (f == "--tolerance" || f == "--delta-min") { a.tolerance = parse_double(need(i), "--tolerance"); ++i; }
        else if (f == "--plot-freq") { a.plot_freq = static_cast<int>(parse_long(need(i), "--plot-freq")); ++i; }
        else if (f == "--epsilon") { a.epsilon = parse_double(need(i), "--epsilon"); ++i; }
        else if (f == "--initial-resolution") { a.initial_resolution = static_cast<int>(parse_long(need(i), "--initial-resolution")); ++i; }
        else if (f == "--refinements") { a.refinements = static_cast<int>(parse_long(need(i), "--refinements")); ++i; }
        else if (f == "--vi-iterations") { a.vi_iterations = static_cast<int>(parse_long(need(i), "--vi-iterations")); ++i; }
        else if (f == "--conservative") { a.conservative = true; }
        else if (f == "--delta-max") { a.delta_max = parse_double(need(i), "--delta-max"); ++i; }
        else if (f == "--workers") { a.workers = static_cast<int>(parse_long(need(i), "--workers")); ++i; }
        else if (f == "--precompute") { a.precompute = true; }
        else if (f == "--out-dir") { a.out_dir = need(i); ++i; }
        else if (f == "--no-plot") { a.no_plot = true; }
        else if (f == "--dump-csv") { a.dump_csv = true; }
        else if (f == "--seq-queries") { a.seq_queries = true; }
        else if (f == "--draw-value-cells") { a.draw_value_cells = true; draw_value_cells_set = true; }
        else if (f == "--dump-slices") { a.dump_slices = static_cast<int>(parse_long(need(i), "--dump-slices")); ++i; }
        else if (f == "--checkpoint-every") { a.checkpoint_every = static_cast<int>(parse_long(need(i), "--checkpoint-every")); ++i; }
        else if (f == "--checkpoint-dir") { a.checkpoint_dir = need(i); ++i; }
        else if (f == "--resume") { a.resume = need(i); ++i; }
        else if (f == "--keep-checkpoints") { a.keep_checkpoints = static_cast<int>(parse_long(need(i), "--keep-checkpoints")); ++i; }
        else { std::fprintf(stderr, "error: unknown argument %s\n", argv[i]); usage_and_exit(argv[0], mode, 2); }
    }
    (void)draw_value_cells_set;

    if (a.algorithm != 1 && a.algorithm != 2) {
        std::fprintf(stderr, "error: --algorithm is required and must be 1 or 2\n");
        usage_and_exit(argv[0], mode, 2);
    }
    if (a.workers <= 0) {
        const unsigned hc = std::thread::hardware_concurrency();
        a.workers = static_cast<int>(hc > 1 ? hc - 1 : 1);
    }
    return a;
}

// ---------------------------------------------------------------------------
// Results directory names, reproduced verbatim from the f-strings in the Python
// so that a C++ run lands beside the corresponding Python run.
// ---------------------------------------------------------------------------

inline std::string vi_output_dir(Mode mode, const Args& a, const Environment& env) {
    std::string suffix;
    if (mode == Mode::RADiscount) {
        suffix = "dynamics_" + a.dynamics_name + "_" +
                 "gamma_" + py_fixed(a.gamma, 3) + "_" +
                 "dt_" + py_fixed(env.dt(), 3) + "_" +
                 "tau_" + py_fixed(env.tau(), 3) + "_" +
                 "tol_" + py_exp(a.tolerance, 1) + "_" +
                 "eps_" + py_fixed(a.epsilon, 3) +
                 "vi-iterations_" + std::to_string(a.vi_iterations) +
                 "conservative_" + py_bool(a.conservative) +
                 "delta-max_" + py_repr(a.delta_max) +
                 "init_resol_" + std::to_string(a.initial_resolution);
    } else {
        suffix = "dynamics_" + a.dynamics_name + "_" +
                 "delta_min_" + py_exp(a.tolerance, 1) + "_" +
                 "delta-max_" + py_repr(a.delta_max) +
                 "init_resol_" + std::to_string(a.initial_resolution) +
                 "vi-iterations_" + std::to_string(a.vi_iterations) +
                 "dt_" + py_fixed(env.dt(), 3) + "_" +
                 "tau_" + py_fixed(env.tau(), 3) + "_";
    }
    return std::string(mode_results_root(mode)) + "/" + suffix;
}

inline std::string refine_output_dir(Mode mode, const Args& a, const Environment& env) {
    std::string suffix;
    if (mode == Mode::RADiscount) {
        suffix = "dynamics_" + a.dynamics_name + "_" +
                 "gamma_" + py_fixed(a.gamma, 3) + "_" +
                 "dt_" + py_fixed(env.dt(), 3) + "_" +
                 "tau_" + py_fixed(env.tau(), 3) + "_" +
                 "delta_min" + py_exp(a.tolerance, 1) + "_" +
                 "vi-iterations_" + std::to_string(a.vi_iterations) +
                 "conservative_" + py_bool(a.conservative) +
                 "delta-max_" + py_repr(a.delta_max) +
                 "init_resol_" + std::to_string(a.initial_resolution);
    } else {
        suffix = "dynamics_" + a.dynamics_name + "_" +
                 "delta_min_" + py_exp(a.tolerance, 1) + "_" +
                 "delta-max_" + py_repr(a.delta_max) +
                 "init_resol_" + std::to_string(a.initial_resolution) +
                 "vi-iterations_" + std::to_string(a.vi_iterations) +
                 "dt_" + py_fixed(env.dt(), 3) + "_" +
                 "tau_" + py_fixed(env.tau(), 3) + "_";
    }
    return std::string(mode_results_root(mode)) + "/" + suffix;
}

inline std::string algorithm1_output_dir(const Args& a, const Environment& env) {
    return "./results/algorithm1_dynamics_" + a.dynamics_name +
           "_resol_" + std::to_string(a.resolution) +
           ",tol_" + py_exp(a.tolerance, 1) +
           "_tau_" + py_fixed(env.tau(), 3) +
           "_dt_" + py_fixed(env.dt(), 3) + "_";
}

}  // namespace dhj
