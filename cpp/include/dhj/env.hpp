// env.hpp - dynamics environments (Dubins car and 3D evasion).
#pragma once

#include <array>
#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>

#include "dim.hpp"
#include "ode.hpp"

namespace dhj {

constexpr double kPi = 3.14159265358979323846;

using Bounds = std::array<std::array<double, 2>, kDim>;  // [dim][lo/hi]
using Bounds3 = Bounds;  // legacy spelling

enum class Dynamics { Dubins, Evasion, DubinsSpeed };

class Environment {
public:
    Environment(double dt, double tau) : dt_(dt), tau_(tau) {
        const double ratio = tau / dt;
        const long n = std::lround(ratio);
        if (std::fabs(ratio - static_cast<double>(n)) > 1e-9 * (1.0 + std::fabs(static_cast<double>(n))))
            throw std::invalid_argument("tau must be evenly divisible by dt");
        n_steps_ = static_cast<int>(n);
    }
    virtual ~Environment() = default;

    virtual const Bounds3& state_bounds() const = 0;
    virtual const std::vector<double>& actions() const = 0;
    virtual State3 ode(const State3& s, double u) const = 0;
    virtual double failure_function(const State3& s) const = 0;
    virtual double reward_function(const State3& s) const = 0;
    virtual const char* name() const = 0;

    double L_f() const { return L_f_; }
    double L_l() const { return L_l_; }
    double L_r() const { return L_r_; }
    double dt() const { return dt_; }
    double tau() const { return tau_; }
    int n_steps() const { return n_steps_; }
    static int state_dim() { return static_cast<int>(kDim); }

    // Applied to each reported state (never fed back into the integrator), exactly
    // as the Python mutates copies of the odeint output rather than its state.
    virtual void normalize(State& s) const { s[kPeriodicDim] = wrap_angle(s[kPeriodicDim]); }

    // Mirrors Environment.dynamics_multi_step: integrate with a constant action
    // and report the state at t = dt, 2*dt, ..., duration (theta wrapped).
    void dynamics_multi_step(const State3& s0, double u, double duration, double dt,
                             std::vector<State3>& out) const {
        const int n = static_cast<int>(std::ceil(duration / dt));
        out.clear();
        out.reserve(static_cast<std::size_t>(n));
        State3 y = s0;
        double t_prev = 0.0;
        for (int i = 1; i <= n; ++i) {
            const double t = duration * static_cast<double>(i) / static_cast<double>(n);
            y = integrate([this, u](const State3& s) { return this->ode(s, u); }, y, t - t_prev);
            t_prev = t;
            State w = y;
            normalize(w);
            out.push_back(w);
        }
    }

    // Mirrors Environment.dynamics: a single step of length tau.
    State3 dynamics(const State3& s0, double u) const {
        State y = integrate([this, u](const State& s) { return this->ode(s, u); }, s0, tau_);
        normalize(y);
        return y;
    }

protected:
    double dt_, tau_;
    int n_steps_ = 1;
    double L_f_ = 1.0, L_l_ = 1.0, L_r_ = 1.0;
};

#if DHJ_DIM == 3

class DubinsCarEnvironment : public Environment {
public:
    DubinsCarEnvironment(double v_const, double dt, double tau, double obstacle_radius = 0.5,
                         double target_radius = 0.5)
        : Environment(dt, tau), v_const_(v_const), obstacle_radius_(obstacle_radius),
          target_radius_(target_radius) {
        bounds_ = {{{{-3.0, 3.0}}, {{-3.0, 3.0}}, {{-kPi, kPi}}}};
        obstacle_position_ = {0.0, 0.0};
        target_position_ = {2.5, 0.0};
        L_f_ = v_const;
        L_l_ = std::sqrt(2.0);
        L_r_ = std::sqrt(2.0);
        actions_ = {-1.0, 0.0, 1.0};
    }

    const Bounds3& state_bounds() const override { return bounds_; }
    const std::vector<double>& actions() const override { return actions_; }
    const char* name() const override { return "DubinsCarEnvironment"; }

    State3 ode(const State3& s, double u) const override {
        return {v_const_ * std::cos(s[2]), v_const_ * std::sin(s[2]), u};
    }

    double failure_function(const State3& s) const override {
        const double dx = s[0] - obstacle_position_[0], dy = s[1] - obstacle_position_[1];
        return std::sqrt(dx * dx + dy * dy) - obstacle_radius_;
    }

    double reward_function(const State3& s) const override {
        const double dx = s[0] - target_position_[0], dy = s[1] - target_position_[1];
        return -(std::sqrt(dx * dx + dy * dy) - target_radius_);
    }

    std::array<double, 2> obstacle_position() const { return obstacle_position_; }
    double obstacle_radius() const { return obstacle_radius_; }
    std::array<double, 2> target_position() const { return target_position_; }
    double target_radius() const { return target_radius_; }

private:
    double v_const_;
    double obstacle_radius_, target_radius_;
    std::array<double, 2> obstacle_position_{}, target_position_{};
    Bounds3 bounds_{};
    std::vector<double> actions_;
};

class EvasionEnvironment : public Environment {
public:
    EvasionEnvironment(double v_const, double dt, double tau, double obstacle_radius = 1.0,
                       double target_radius = 0.5)
        : Environment(dt, tau), v_const_(v_const), obstacle_radius_(obstacle_radius),
          target_radius_(target_radius) {
        bounds_ = {{{{-3.0, 3.0}}, {{-3.0, 3.0}}, {{-kPi, kPi}}}};
        obstacle_position_ = {0.0, 0.0};
        target_position_ = {2.5, 0.0};
        L_f_ = 1.0 + v_const;
        L_l_ = std::sqrt(2.0);
        L_r_ = std::sqrt(2.0);
        actions_.resize(5);
        for (int i = 0; i < 5; ++i) actions_[i] = -1.0 + 2.0 * i / 4.0;  // np.linspace(-1, 1, 5)
    }

    const Bounds3& state_bounds() const override { return bounds_; }
    const std::vector<double>& actions() const override { return actions_; }
    const char* name() const override { return "EvasionEnvironment"; }

    State3 ode(const State3& s, double u) const override {
        const double v = v_const_;
        return {-v + v * std::cos(s[2]) + u * s[1],
                v * std::sin(s[2]) - u * s[0],
                -u};
    }

    double failure_function(const State3& s) const override {
        const double dx = s[0] - obstacle_position_[0], dy = s[1] - obstacle_position_[1];
        return std::sqrt(dx * dx + dy * dy) - obstacle_radius_;
    }

    double reward_function(const State3& s) const override {
        const double dx = s[0] - target_position_[0], dy = s[1] - target_position_[1];
        return -(std::sqrt(dx * dx + dy * dy) - target_radius_);
    }

    std::array<double, 2> obstacle_position() const { return obstacle_position_; }
    double obstacle_radius() const { return obstacle_radius_; }
    std::array<double, 2> target_position() const { return target_position_; }
    double target_radius() const { return target_radius_; }

private:
    double v_const_;
    double obstacle_radius_, target_radius_;
    std::array<double, 2> obstacle_position_{}, target_position_{};
    Bounds3 bounds_{};
    std::vector<double> actions_;
};

#endif  // DHJ_DIM == 3

#if DHJ_DIM == 4

// Dubins car with speed as a fourth state: (x, y, theta, v).
//
//   dx/dt = v cos(theta),  dy/dt = v sin(theta),  dtheta/dt = omega,  dv/dt = a
//
// Mirrors DubinsSpeedEnvironment in Scratchpad/4DNewDynamics/dynamics_4d.py.
// The Jacobian is nonzero only in rows 0-1 (columns theta and v), bounded by
// |v| <= v_max and 1, so the max-row-sum Lipschitz constant is L_f = v_max + 1.
//
// The action is a PAIR (omega, a). The engine carries actions as opaque doubles,
// so the pair table is indexed and the index travels as the action value - the
// same trick the Python gets for free by passing tuples around.
class DubinsSpeedEnvironment : public Environment {
public:
    DubinsSpeedEnvironment(double dt, double tau, double obstacle_radius = 1.3,
                           double target_radius = 0.5)
        : Environment(dt, tau), obstacle_radius_(obstacle_radius), target_radius_(target_radius) {
        bounds_ = {{{{-3.0, 3.0}}, {{-3.0, 3.0}}, {{-kPi, kPi}}, {{0.2, 1.5}}}};
        obstacle_position_ = {0.0, 0.0};
        target_position_ = {2.5, 0.0};

        const double v_max = bounds_[3][1];
        L_f_ = v_max + 1.0;
        L_l_ = std::sqrt(2.0);
        L_r_ = std::sqrt(2.0);

        // itertools.product((-1,0,1), (-0.5,0,0.5)) - omega outer, accel inner.
        const double omegas[3] = {-1.0, 0.0, 1.0};
        const double accels[3] = {-0.5, 0.0, 0.5};
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j) pairs_.push_back({omegas[i], accels[j]});
        actions_.resize(pairs_.size());
        for (std::size_t i = 0; i < pairs_.size(); ++i) actions_[i] = static_cast<double>(i);
    }

    const Bounds& state_bounds() const override { return bounds_; }
    const std::vector<double>& actions() const override { return actions_; }
    const char* name() const override { return "DubinsSpeedEnvironment"; }

    State ode(const State& s, double u) const override {
        const std::size_t k = static_cast<std::size_t>(u + 0.5);
        const double omega = pairs_[k][0], accel = pairs_[k][1];
        // v is clamped where it is USED, matching the Python vector field; v itself
        // still integrates freely (dv/dt = a) and is clamped only on output.
        const double v = std::min(std::max(s[3], bounds_[3][0]), bounds_[3][1]);
        return {v * std::cos(s[2]), v * std::sin(s[2]), omega, accel};
    }

    // theta wraps, v saturates at its bounds
    void normalize(State& s) const override {
        s[2] = wrap_angle(s[2]);
        s[3] = std::min(std::max(s[3], bounds_[3][0]), bounds_[3][1]);
    }

    double failure_function(const State& s) const override {
        const double dx = s[0] - obstacle_position_[0], dy = s[1] - obstacle_position_[1];
        return std::sqrt(dx * dx + dy * dy) - obstacle_radius_;
    }

    double reward_function(const State& s) const override {
        const double dx = s[0] - target_position_[0], dy = s[1] - target_position_[1];
        return -(std::sqrt(dx * dx + dy * dy) - target_radius_);
    }

    std::array<double, 2> obstacle_position() const { return obstacle_position_; }
    double obstacle_radius() const { return obstacle_radius_; }
    std::array<double, 2> target_position() const { return target_position_; }
    double target_radius() const { return target_radius_; }

private:
    double obstacle_radius_, target_radius_;
    std::array<double, 2> obstacle_position_{}, target_position_{};
    Bounds bounds_{};
    std::vector<std::array<double, 2>> pairs_;
    std::vector<double> actions_;
};

#endif  // DHJ_DIM == 4

}  // namespace dhj
