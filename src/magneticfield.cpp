#include "magneticfield.h"
#include "SPH.h"
#define _USE_MATH_DEFINES
#include <math.h>

using namespace std;
using namespace Eigen;

MagneticField::MagneticField(double h):
    m_h(h)
{

}

 VectorXd MagneticField::calculateMagneticField(const _SPHApi &particles){
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


void MagneticField::buildProblem(MatrixXd& mat, const _SPHApi &particles){

    double V = particles.getVolume();
    double chi = particles.getMagneticSusceptibility();
    double Gamma = V * chi / (1 + chi);
    // room for paralellization
    for (int particle = 0; particle <  particles.getParticles(); ++particle) {
        for (int neighbor = 0; neighbor <  particles.getParticles(); ++neighbor) {
            
            Vector3d r_ik = particles.getPos(particle) -  particles.getPos(neighbor);
            double l_ik = r_ik.norm();

            double W_avr_ik = W_avr(l_ik);
            double W_ik = W(l_ik);

            for (const int& j : {0, 1, 2}){
                for (const auto& l : {0, 1, 2}){
                    if (particle != neighbor || j == l){
                        double G  = 0.0;
                        if (particle != neighbor) G += r_ik(l) * r_ik(j) / pow(l_ik, 2.0) * (W_avr_ik - W_ik);
                        if (j == l) G += W_ik - W_avr_ik /3.0;

                        int idx1 = particle * 3 + j;
                        int idx2 = neighbor * 3 + l;

                        mat(idx1, idx2) = delta(idx1, idx2) - Gamma * G;
                    }
                }
            }
        }
    }
}





double MagneticField::w (const double q) const {
    if (0 <= q && q < 1){
        return 0.25 * pow(2.0 - q, 3.0) - pow(1.0 - q, 3.0);
    } else if (1 <= q && q < 2){
        return 0.25 * pow(2.0 - q, 3.0);
    } else {
        return 0.0;
    }
}

double MagneticField::w_avr (const double q) const {
    if (0 <= q && q < 1){
        return 1/ (40.0 * M_PI) * (15.0 * pow(q, 3) - 36.0 * pow(q, 1.0) + 40.0);
    } else if (1 <= q && q < 2){
        return -3.0/(4.0 * M_PI * pow(q, 3.0)) * (pow(q, 6.0)/6.0 - 6.0 * pow(q, 5.0)/5.0 
        + 3.0 * pow(q, 4.0) - 8.0 * pow(q, 3.0) / 3.0+ 1.0/15.0);
    } else {
        return 3.0/(4.0 * pow(q, 3.0));
    } 
}

double MagneticField::W (const double q) const {
    return w(q/ m_h) / pow(m_h, 3.0);
}

double MagneticField::delta(const int i,const int j) const {
    return i == j ? 1.0 : 0.0;
}

double MagneticField::W_avr (const double q) const{
    return w_avr(q/ m_h) / pow(m_h, 3.0);
}

