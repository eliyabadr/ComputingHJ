// dim.hpp - compile-time state-space dimension.
//
// The port was written for the paper's two 3D systems and baked `3` into the
// state/bounds types, the BVH and the reach-box construction. Rather than
// template every core type on N (which would touch every header and every call
// site), the dimension is a single compile-time constant: build the normal
// binaries at DHJ_DIM=3 and a separate binary at DHJ_DIM=4. The 3D build is
// bit-identical to before - nothing about it is conditional.
#pragma once

#include <cstddef>

#ifndef DHJ_DIM
#define DHJ_DIM 3
#endif

namespace dhj {

constexpr std::size_t kDim = DHJ_DIM;

// Index of the 2*pi-periodic (wrapping) state dimension. Both 3D systems and
// the 4D Dubins-with-speed put theta at index 2; 4D adds v at index 3, which is
// NOT periodic. Everything that wraps, and everything that decides whether a
// reach box left the domain, keys off this.
constexpr std::size_t kPeriodicDim = 2;

static_assert(kDim == 3 || kDim == 4, "DHJ_DIM must be 3 or 4");
static_assert(kPeriodicDim < kDim, "periodic dimension must be in range");

}  // namespace dhj
