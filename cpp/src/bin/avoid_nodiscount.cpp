// avoid_nodiscount - C++ port of avoid_nodiscount.py
#include "dhj/run.hpp"

int main(int argc, char** argv) {
    return dhj::run_main(dhj::Mode::AvoidNoDiscount, argc, argv);
}
