#include "magneticfield.h"
#include <math.h>

using namespace std;
using namespace Eigen;

MagneticField::MagneticField(double h):
    m_h(h)
{

}


void MagneticField::buildProblem(SpMat& mat,const SPH &particles){
    
    vector<T> coefficients; 

    for (const int &particle: particles.getParticles()){
        for (const int& neighbor: particles.getNeighbors(particle, 2.0 * m_h)){
            
            Vector3f r_ik = particles.getPos(particle) -  particles.getPos(neighbor);
            float l_ik = r_ik.norm();

            float W_avr_ik = W_avr(l_ik);
            float W_ik = W_avr(l_ik);

            for (const int& j : {0, 1, 2}){
                for (const auto& l : {0, 1, 2}){
                    if (particle != neighbor || j == l){
                        float G  = 0.0;
                        if (particle != neighbor) G += r_ik(l) * r_ik(j) / l_ik * (W_avr_ik - W_ik);
                        if (j == l) G += W_ik - W_avr_ik /3.0;
                        
                        int idx1 = particle * 3 + j;
                        int idx2 = neighbor * 3 + l;

                        coefficients.push_back(T(idx1, idx2, G));
                    }
                }
            }
        }
    }

    mat.setFromTriplets(coefficients.begin(), coefficients.end());
}





const double MagneticField::w (const double q) const {
    if (0 <= q && q < 1){
        return 0.25 * pow(2.0 - q, 3.0) - pow(1.0 - q, 3.0);
    } else if (1 <= q && q < 2){
        return 0.25 * pow(2.0 - q, 3.0);
    } else {
        return 0.0;
    }
}

const double MagneticField::w_avr (const double q) const {
    if (0 <= q && q < 1){
        return 1/ (40.0 * M_PI) * (15.0 * pow(q, 3) - 36.0 * pow(q, 1.0) + 40.0);
    } else if (1 <= q && q < 2){
        return -3.0/(4.0 * M_PI * pow(q, 3.0)) * (pow(q, 6.0)/6.0 - 6.0 * pow(q, 5.0)/5.0 
        + 3.0 * pow(q, 4.0) - 8.0 * pow(q, 3.0) / 3.0+ 1.0/15.0);
    } else {
        return 0.0;
    } 
}

const double MagneticField::W (const double q) const {
    return w(q/ m_h) / pow(m_h, 3.0);
}

const double MagneticField::W_avr (const double q) const{
    return w_avr(q/ m_h) / pow(m_h, 3.0);
}


