# Plain-make build, for when cmake is not on hand.
CXX      ?= c++
CXXFLAGS ?= -std=c++17 -O3 -DNDEBUG -march=native -fno-math-errno -Iinclude -Wall -Wextra
LDFLAGS  ?= -pthread

BINS := build/ra_nodiscount build/ra_discount build/avoid_nodiscount
DIM4 := build/ra_nodiscount_4d
HDRS := $(wildcard include/dhj/*.hpp)

.PHONY: all clean dim4
all: $(BINS)

build/%: src/bin/%.cpp $(HDRS)
	@mkdir -p build
	$(CXX) $(CXXFLAGS) $< -o $@ $(LDFLAGS)

# 4D build (state = x, y, theta, v). Same sources, DHJ_DIM=4.
dim4: $(DIM4)

build/ra_nodiscount_4d: src/bin/ra_nodiscount.cpp $(HDRS)
	@mkdir -p build
	$(CXX) $(CXXFLAGS) -DDHJ_DIM=4 $< -o $@ $(LDFLAGS)

# Sanitizer build: AddressSanitizer + UndefinedBehaviorSanitizer.
# LeakSanitizer is not available on macOS/arm64 -- use /usr/bin/leaks instead.
.PHONY: san
san: build_san/ra_nodiscount_asan
build_san/ra_nodiscount_asan: src/bin/ra_nodiscount.cpp $(HDRS)
	@mkdir -p build_san
	$(CXX) -std=c++17 -O1 -g -fsanitize=address,undefined \
	    -fno-omit-frame-pointer -fno-sanitize-recover=undefined \
	    -Iinclude -Wall -Wextra $< -o $@ $(LDFLAGS)

# Hardened release: libc++ bounds-checks operator[] on the standard containers
# and the stack gets a guard. Slower than `all`; use it to validate a change.
.PHONY: hardened
hardened: build_hardened/ra_nodiscount
build_hardened/ra_nodiscount: src/bin/ra_nodiscount.cpp $(HDRS)
	@mkdir -p build_hardened
	$(CXX) -std=c++17 -O2 -g -march=native -fno-math-errno \
	    -D_LIBCPP_HARDENING_MODE=_LIBCPP_HARDENING_MODE_EXTENSIVE \
	    -fstack-protector-strong -Iinclude -Wall -Wextra $< -o $@ $(LDFLAGS)

clean:
	rm -rf build build_san build_hardened
