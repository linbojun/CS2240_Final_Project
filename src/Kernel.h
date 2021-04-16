#pragma once

//#include "core/Common.h"
//#include "core/Vector.h"
#include <math.h>
#include <cmath>
#include <Eigen/Dense>
#include <Eigen/StdVector>



#define cube(a) (a)*(a)*(a)
#define sqr(a) (a)*(a)

using namespace Eigen;
// SPH Kernels
// Kernels are split into constant and variable part.
// Note: r is supposed to be within the filter support (e.g. |r| <= h)
// Arguments are as follows:
// r  = displacement vector
// r2 = |r|^2 (squared norm of r)
// rn = |r|   (norm of r)
struct Kernel {
    double h;
    double h2;
    double halfh;

    void init(double h_) {
        h = h_;
        h2 = h*h;
        halfh = 0.5f * h;
        poly6Constant = 365.f / (64.f * M_PI * std::pow(h, 9.f));
        poly6GradConstant = -945.f / (32.f * M_PI * std::pow(h, 9.f));
        poly6LaplaceConstant = -945.f / (32.f * M_PI * std::pow(h, 9.f));
        spikyConstant = 15.f / (M_PI * std::pow(h, 6.f));
        spikyGradConstant = -45.f / (M_PI * std::pow(h, 6.f));
        spikyLaplaceConstant = -90.f / (M_PI * std::pow(h, 6.f));
        viscosityLaplaceConstant = 45.f / (M_PI * std::pow(h, 6.f));

        surfaceTensionConstant = 32.f / (M_PI * std::pow(h, 9.f));
        surfaceTensionOffset = -std::pow(h, 6.f) / 64.f;
    }

    double poly6Constant;
    inline double poly6(double r2) const {
        return cube(h2 - r2);
    }

    double poly6GradConstant;
    inline Vector3d poly6Grad(const Vector3d &r, double r2) const {
        return sqr(h2 - r2) * r;
    }

    double poly6LaplaceConstant;
    inline double poly6Laplace(double r2) {
        return (h2 - r2) * (3.f * h2 - 7.f * r2);
    }

    double spikyConstant;
    inline double spiky(double rn) const {
        return cube(h - rn);
    }

    double spikyGradConstant;
    inline Vector3d spikyGrad(const Vector3d &r, double rn) const {
        return sqr(h - rn) * r * (1.f / rn);
    }

    double spikyLaplaceConstant;
    inline double spikyLaplace(double rn) const {
        return (h - rn) * (h - 2.f * rn) / rn;
    }

    double viscosityLaplaceConstant;
    inline double viscosityLaplace(double rn) const {
        return (h - rn);
    }

    double surfaceTensionConstant;
    double surfaceTensionOffset;
    inline double surfaceTension(double rn) const {
        if (rn < halfh) {
            return 2.f * cube(h - rn) * cube(rn) + surfaceTensionOffset;
        } else {
            return cube(h - rn) * cube(rn);
        }
    }
};

