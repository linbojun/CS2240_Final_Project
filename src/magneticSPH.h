
#ifndef MAGNETICSPH_H
#define MAGNETICSPH_H
#include "Eigen/Dense"
#include "Eigen/Sparse"
#include "Eigen/SparseQR"
#include "magneticinit.h"
#include "SPH.h"
#define _USE_MATH_DEFINES
#include <cmath>
#include <math.h>



const float _chi = 0.1f;
const float _mu_0 = 4 * M_PI * 1e-7;

class MagneticSPH : public SPH
{

public:
    
    MagneticSPH(int n, float radius, double h);
    void update(float time_step) override;

private:
    int m_t = 0;
    int m_subupdate = 1;


    Eigen::VectorXd m_magneticForce;
    Eigen::VectorXd m_Bext;
    double m_h;

    Eigen::VectorXd m_guess;
    bool m_isFirst = true;

    Eigen::ConjugateGradient<Eigen::MatrixXd, Eigen::Lower| Eigen::Upper> m_cg;

    double w (const double q) const;
    double w_avr (const double q) const;
    
    double W (const double q) const;
    double W_avr (const double q) const;

    double Wprime (const double q) const;
    double wprime (const double q) const;

    double A (const double q) const;
    double Aprime (const double q) const;

    double C1 (const double q) const;
    double C2 (const double q) const;

    Eigen::Matrix3d objectToWorld (const Eigen::Vector3d r) const;

    Eigen::Matrix3d delH (const Eigen::Vector3d r, const Eigen::Vector3d m) const;


    double delta (const int i,const int j) const;

    Eigen::VectorXd calculateMagneticField(const Eigen::MatrixXd& A);

    Eigen::VectorXd calculateMagneticForce(Eigen::VectorXd &F);

    void buildProblem(Eigen::MatrixXd &mat);
    
    
    float getG (const SPH &particles, int i, int j, int k, int l);

    double getMagneticSusceptibility() const;
    Eigen::VectorXd getExternalB() const;
    double getPermeability() const;
    double getGamma() const;

};

#endif
