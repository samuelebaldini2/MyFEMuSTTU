#pragma once

#include <vector>
#include <stdexcept>
#include <cstddef>

#if defined(__GNUC__) || defined(__clang__)
  #define FEMUS_RESTRICT __restrict__
#else
  #define FEMUS_RESTRICT
#endif

class AdvectMarkers {
public:
  // Avanza i marker: X^{n} -> X^{n+1}
  static inline void forward(std::vector<std::vector<double>>& Xp,
                             const std::vector<std::vector<double>>& markerVelocity,
                             const double dt)
  {
    if (dt == 0.0) return;
    advect_inplace(Xp, markerVelocity, dt);
  }

  // Torna indietro: X^{n+1} -> X^{n}
  static inline void backward(std::vector<std::vector<double>>& Xp,
                              const std::vector<std::vector<double>>& markerVelocity,
                              const double dt)
  {
    if (dt == 0.0) return;
    advect_inplace(Xp, markerVelocity, -dt);
  }

private:
  static inline void advect_inplace(std::vector<std::vector<double>>& Xp,
                                    const std::vector<std::vector<double>>& V,
                                    const double dt)
  {
    const std::size_t dim = Xp.size();

    if (dim != 2 && dim != 3) {
      throw std::runtime_error("AdvectMarkers::advect_inplace: Xp.size() must be 2 or 3");
    }

    if (V.size() != dim) {
      throw std::runtime_error("AdvectMarkers::advect_inplace: velocity dimension mismatch");
    }

    if (Xp[0].empty()) return;

    const std::size_t nPts = Xp[0].size();

    for (std::size_t d = 1; d < dim; ++d) {
      if (Xp[d].size() != nPts) {
        throw std::runtime_error("AdvectMarkers::advect_inplace: inconsistent Xp[d].size()");
      }
      if (V[d].size() != nPts) {
        throw std::runtime_error("AdvectMarkers::advect_inplace: inconsistent V[d].size()");
      }
    }

    if (V[0].size() != nPts) {
      throw std::runtime_error("AdvectMarkers::advect_inplace: inconsistent V[0].size()");
    }

    double* FEMUS_RESTRICT x  = Xp[0].data();
    double* FEMUS_RESTRICT y  = Xp[1].data();
    double* FEMUS_RESTRICT z  = (dim == 3) ? Xp[2].data() : nullptr;

    const double* FEMUS_RESTRICT u = V[0].data();
    const double* FEMUS_RESTRICT v = V[1].data();
    const double* FEMUS_RESTRICT w = (dim == 3) ? V[2].data() : nullptr;

    if (dim == 2) {
      for (std::size_t i = 0; i < nPts; ++i) {
        x[i] += dt * u[i];
        y[i] += dt * v[i];
      }
      return;
    }

    for (std::size_t i = 0; i < nPts; ++i) {
      x[i] += dt * u[i];
      y[i] += dt * v[i];
      z[i] += dt * w[i];
    }
  }
};

#undef FEMUS_RESTRICT
