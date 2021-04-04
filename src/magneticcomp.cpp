#include "magneticcomp.h"
#include "SPH.h"
#define _USE_MATH_DEFINES
#include <math.h>

using namespace std;
using namespace Eigen;

MagneticComp::MagneticComp(double h):
    m_h(h)
{

}

 VectorXd MagneticComp::calculateMagneticField(const SPH &particles){
     MatrixXd A;
     buildProblem(A, particles);
     VectorXd h_ext = particles.getExternalB();
     m_cg.compute(A);
     if (m_isFirst){
         m_isFirst = false;
         m_guess = m_cg.solve(h_ext);
     } else {
         m_guess = m_cg.solveWithGuess(h_ext, m_guess);
     }
     return m_guess;
 }


void MagneticComp::buildProblem(MatrixXd& mat, const SPH &particles){
    double Gamma = particles.getGamma();
    // room for paralellization
    for (int particle = 0; particle <  particles.getNumParticle(); ++particle) {
        for (int neighbor = 0; neighbor <  particles.getNumParticle(); ++neighbor) {
            
            Vector3d r_ik = particles.getPos(particle) -  particles.getPos(neighbor);
            double l_ik = r_ik.norm();

            double W_avr_ik = W_avr(l_ik);
            double W_ik = W(l_ik);

            for (const int& j : {0, 1, 2}){
                for (const auto& l : {0, 1, 2}){
                    if (particle != neighbor || j == l){
                        double G  = 0.0;
                        if (particle != neighbor) G += r_ik(l) * r_ik(j) / pow(l_ik, 2.0) * (W_avr_ik - W_ik);

                        if (j == l) {
                            if (particle != neighbor) {
                                G += W_ik - W_avr_ik /3.0;
                            } else {
                                G += 2.0 * W_ik /3.0;
                            }

                        }
                        
                        int idx1 = particle * 3 + j;
                        int idx2 = neighbor * 3 + l;

                        mat(idx1, idx2) = delta(idx1, idx2) - Gamma * G;
                    }
                }
            }
        }
    }
}

VectorXd MagneticComp::calculateMagneticForce(VectorXd &F, const SPH &particles) {
    VectorXd m = particles.getGamma() * calculateMagneticField(particles);
    double mu_0 = particles.getPermeability();
    // room for paralellization
    for (int target = 0; target < particles.getNumParticle(); ++target) {
        Matrix3d U = Matrix3d::Zero();
        Vector3d m_target = Vector3d(m(3*target), m(3*target + 1), m(3*target + 2));
        for (int source = 0; source < particles.getNumParticle(); ++source) {
            Vector3d r = particles.getPos(target) - particles.getPos(source);
            double l = r.norm();
            Vector3d m_source = Vector3d(m(3*source), m(3*source + 1), m(3*source + 2));
            if (l >= 4 * m_h) {
                U += mu_0 * delH(r, m_source);
            } else {
                double q = l/m_h;
                Matrix3d R = objectToWorld(r);
                Vector3d local_m_source = R.transpose() * m_source;
                double c1 = C1(q)/pow(m_h, 4.0);
                double c2 = C2(q)/pow(m_h, 4.0);
                Matrix3d local_tensor;
                local_tensor.col(0) = Vector3d(c1*local_m_source(0), 0.0, c1*local_m_source(0));
                local_tensor.col(1) = Vector3d(0.0, c1*local_m_source(2), c1*local_m_source(1));
                local_tensor.col(2) = Vector3d(c1*local_m_source(0), c1*local_m_source(1), c2*local_m_source(2));
                U += R * local_tensor * R.transpose();
            }
        }
        Vector3d targetForce = U * m_target;
        F(3*target) = targetForce(0);
        F(3*target+1) = targetForce(1);
        F(3*target+2) = targetForce(2);
    }
    return F;
}

Matrix3d MagneticComp::objectToWorld (const Vector3d r) const {
    Vector3d v3 = r.normalized();
    Vector3d v2 = Vector3d(0.0, 1.0, 0.0);
    if (fabs(v3.dot(v2) == 1.0)) {
        v2 = Vector3d(1.0, 0.0, 0.0);
    }
    v2 = (v2 - v2.dot(v3)*v2).normalized();
    Vector3d v1 = v2.cross(v3);
    Matrix3d R;
    R << v1, v2, v3;
    return R;
}

Matrix3d MagneticComp::delH (const Vector3d r, const Vector3d m) const {
    double l = r.norm();
    return (r.transpose() * m * Matrix3d::Identity() + r * m.transpose() + m * r.transpose()) * A(l) + (r * (r.transpose() * m * r.transpose())/r.norm()) * Aprime(l);
}

double MagneticComp::w (const double q) const {
    if (0 <= q && q < 1){
        return 0.25 * pow(2.0 - q, 3.0) - pow(1.0 - q, 3.0) / M_PI;
    } else if (1 <= q && q < 2){
        return 0.25 * pow(2.0 - q, 3.0) / M_PI;
    } else {
        return 0.0;
    }
}

double MagneticComp::w_avr (const double q) const {
    if (0 <= q && q < 1){
        return (15.0 * pow(q, 3) - 36.0 * pow(q, 1.0) + 40.0)/ (40.0 * M_PI) ;
    } else if (1 <= q && q < 2){
        return -3.0/(4.0 * M_PI * pow(q, 3.0)) * (pow(q, 6.0)/6.0 - 6.0 * pow(q, 5.0)/5.0 
        + 3.0 * pow(q, 4.0) - 8.0 * pow(q, 3.0) / 3.0+ 1.0/15.0);
    } else {
        return 3.0/(4.0 * pow(q, 3.0)* M_PI);
    } 
}

double MagneticComp::W (const double q) const {
    return w(q/ m_h) / pow(m_h, 3.0);
}

double MagneticComp::delta(const int i,const int j) const {
    return i == j ? 1.0 : 0.0;
}

double MagneticComp::W_avr (const double q) const{
    return w_avr(q/ m_h) / pow(m_h, 3.0);
}

double MagneticComp::wprime (const double q) const {
    if (0 <= q && q < 1) {
        return (-0.75 * pow(2.0-q, 2.0) + 3.0 * pow(1.0-q, 2.0))/M_PI;
    } else if (1 <= q && q < 2) {
        return (-0.75 * pow(2.0-q, 2.0))/M_PI;
    } else {
        return 0.0;
    }
}

double MagneticComp::Wprime (const double q) const  {
    return wprime(q/m_h)/pow(m_h, 4.0);
}

double MagneticComp::A (const double q) const  {
    return (W(q) - W_avr(q))/pow(q, 2.0);
}

double MagneticComp::Aprime (const double q) const  {
    return 5.0 * (W(q) - W_avr(q))/pow(q, 3.0) - Wprime(q)/pow(q, 2.0);
}

double MagneticComp::C1 (const double q) const  {
    VectorXd qexp(5);
    qexp << pow(q, 4.0), pow(q, 3.0), pow(q, 2.0), pow(q, 1.0), 1.0;
    VectorXd coeff(5);
    if (q >= 0 && q < 1) {
        coeff << 9.978e-9, -2.979e-8, 2.389e-9, 4.531e8, 2.446e-11;
    } else if (q >=1 && q < 2) {
        coeff << -2.764e-9, 2.869e-8, -9.945e-8, 1.251e-7, -2.370e-8;
    } else if (q >=2 && q < 3) {
        coeff << -1.096e-9, 9.770e-9, -2.547e-8, 2.650e-9, 5.007e-8;
    } else if (q >=3 && q < 4) {
        coeff << 3.799e-10, -6.263e-9, 3.947e-8, -1.135e-7, 1.274e-7;
    }
    return qexp.dot(coeff);
}

double MagneticComp::C2 (const double q) const  {
    VectorXd qexp(5);
    qexp << pow(q, 4.0), pow(q, 3.0), pow(q, 2.0), pow(q, 1.0), 1.0;
    VectorXd coeff(5);
    if (q >= 0 && q < 1) {
        coeff << 6.695e-8, -1.617e-7, 1.682e-8, 1.345e-7, 1.109e-10;
    } else if (q >=1 && q < 2) {
        coeff << -3.084e-8, 2.291e-7, -5.883e-7, 5.611e-7, -1.144e-10;
    } else if (q >=2 && q < 3) {
        coeff << 3.504e-9, -5.259e-8, 2.788e-7, -6.241e-7, 4.918e-7;
    } else if (q >=3 && q < 4) {
        coeff << 7.334e-10, -9.588e-9, 4.379e-8, -7.480e-8, 2.341e-8;
    }
    return qexp.dot(coeff);
}
