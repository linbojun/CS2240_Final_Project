#include "magneticwcsph.h"
#define _USE_MATH_DEFINES
#include <math.h>
#include <iostream>
#include <Eigen/Core>
#include "shapes.h"

using namespace std;
using namespace Eigen;
const float _chi = 3.f;
const float _mu_0 = 4 * M_PI * 1e-7;

// TODO: n, radius are fully ignored
MagneticWCSPH::MagneticWCSPH(int n, float radius, double h):
    WCSPH(radius, 1), m_h(h), m_subupdate(5), m_Binit()
{
    m_isFirst = true;
    m_Bext = VectorXd::Zero(3 * this->getNumParticle());
    m_guess = VectorXd::Zero(3 * this->getNumParticle());
//    Binit.addConstField(Vector3d(0.0, 2e-4, 0.0));
    m_Binit.addPointSource(Vector3d(0.5, -1, 0.5), Vector3d(0, 2, 0));

    for (int particle = 0; particle < getNumParticle(); ++particle){
        Vector3d bExt = m_Binit.getMagneticField(getPos(particle));
        for (int i : {0, 1, 2}){
            m_Bext(3 * particle + i) = bExt(i);
        }
    }


    m_magneticForce = VectorXd::Zero(3 * this->getNumParticle());
    calculateMagneticForce(m_magneticForce);


}

void MagneticWCSPH::update_velocity_position()
{
    if (m_t % m_subupdate == 0){
        //cout << "update Magnet" << endl;
        calculateMagneticForce(m_magneticForce);
    }
    //#pragma omp parallel for
    for(int i = 0; i < _fluid_ptcl_list.size(); i++)
    {
        auto cur_fluid = _fluid_ptcl_list[i];
        //if(!cur_fluid->shouldSim)
        //    continue;
        Vector3d f = cur_fluid->netForce;
        //cout << "f: " << f(0) <<", " << f(1) <<", " << f(2) <<endl;
        Vector3d totalForce = cur_fluid->netForce+ 1000 * Vector3d(m_magneticForce(3*i), m_magneticForce(3*i+1), m_magneticForce(3*i+2));
        Vector3d dvdt = totalForce / fluid_ptcl_mass;
        cur_fluid->velocity += dvdt * dt;
        updateParticlePos(i, cur_fluid->position + cur_fluid->velocity * dt);
    }
    m_t++;
}

 VectorXd MagneticWCSPH::calculateMagneticField(const MatrixXd& A){
     VectorXd h_ext = getExternalB();
     m_cg.compute(A);
     if (m_isFirst){
         m_isFirst = false;
         m_guess = m_cg.solve(h_ext);
     } else {
         m_guess = m_cg.solveWithGuess(h_ext, m_guess);
     }
     return m_guess;
 }


void MagneticWCSPH::buildProblem(MatrixXd& mat){
    double Gamma = getGamma();
    //cout << "getGamma" << endl;
    // room for paralellization
    //cout << "num:" << getNumParticle() << endl;
    #pragma omp parallel for
    for (int particle = 0; particle <  getNumParticle(); ++particle) {
        #pragma omp parallel for
        for (int neighbor = 0; neighbor <  getNumParticle(); ++neighbor) {

            Vector3d r_ik = getPos(particle) -  getPos(neighbor);
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

//                        cout << "index" <<idx1 << " " << idx2 << endl;
//                        cout << "delta" << delta(idx1, idx2) << endl;
//                        cout << "G" << Gamma * G << endl;
                        mat(idx1, idx2) = delta(idx1, idx2) - Gamma * G;
//                        cout << "assigned" << endl;
                    }
                }
            }
        }
    }
}

VectorXd MagneticWCSPH::calculateMagneticForce(VectorXd &F) {

//    cout << "calculateMagneticForce" << endl;
    MatrixXd mat(3 * getNumParticle(), 3 * getNumParticle());
//    cout << "insytantiate matrix" << endl;
    buildProblem(mat);
//     cout << "build Problem" << endl;
    VectorXd b = calculateMagneticField(mat);
    VectorXd m = getGamma() * b;
    for (int i = 0; i < getNumParticle(); ++i){
//        cout << "b: (" << b(3*i) << ", " << b(3*i+1) << ", " << b(3*i +2) << ")" << endl;
    }
//     cout << "build m" << endl;
    double mu_0 = getPermeability();
    // room for paralellization
    for (int target = 0; target < getNumParticle(); ++target) {
        Matrix3d U = Matrix3d::Zero();
        Vector3d m_target = Vector3d(m(3*target), m(3*target + 1), m(3*target + 2));
        #pragma omp parallel for
        for (int source = 0; source < getNumParticle(); ++source) {
            Vector3d r = getPos(target) - getPos(source);
            if(r.norm() < 0.000001)
                continue;
            double l = r.norm();
            Vector3d m_source = Vector3d(m(3*source), m(3*source + 1), m(3*source + 2));
            if (true) {
//            if (l >= 4 * m_h) {
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

    //External Force
    for (int target = 0; target < getNumParticle(); ++target) {
        Vector3d m_target(m(3 * target), m(3 * target + 1), m(3 * target+ 2));
        Vector3d externalForce = mu_0 * m_Binit.getMagneticFieldGrad(getPos(target)) * m_target;

        //cout << "f:" << externalForce(0) << "," << externalForce(1) << "," << externalForce(2) << endl;
        F(3*target) += externalForce(0);
        F(3*target+1) += externalForce(1);
        F(3*target+2) += externalForce(2);
        //cout << "f: " << F(3*target) <<", " << F(3*target+1) <<", " << F(3*target + 2) << endl;
    }



    return F;
}

Matrix3d MagneticWCSPH::objectToWorld (const Vector3d r) const {
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

Matrix3d MagneticWCSPH::delH (const Vector3d r, const Vector3d m) const {
    double l = r.norm();
    return (r.transpose() * m * Matrix3d::Identity() + r * m.transpose() + m * r.transpose()) * A(l) + (r * (r.transpose() * m * r.transpose())/r.norm()) * Aprime(l);
}

double MagneticWCSPH::w (const double q) const {
    if (0 <= q && q < 1){
        return 0.25 * pow(2.0 - q, 3.0) - pow(1.0 - q, 3.0) / M_PI;
    } else if (1 <= q && q < 2){
        return 0.25 * pow(2.0 - q, 3.0) / M_PI;
    } else {
        return 0.0;
    }
}

double MagneticWCSPH::w_avr (const double q) const {
    if (0 <= q && q < 1){
        return (15.0 * pow(q, 3) - 36.0 * pow(q, 1.0) + 40.0)/ (40.0 * M_PI) ;
    } else if (1 <= q && q < 2){
        return -3.0/(4.0 * M_PI * pow(q, 3.0)) * (pow(q, 6.0)/6.0 - 6.0 * pow(q, 5.0)/5.0
        + 3.0 * pow(q, 4.0) - 8.0 * pow(q, 3.0) / 3.0+ 1.0/15.0);
    } else {
        return 3.0/(4.0 * pow(q, 3.0)* M_PI);
    }
}

double MagneticWCSPH::W (const double q) const {
    return w(q/ m_h) / pow(m_h, 3.0);
}

double MagneticWCSPH::delta(const int i,const int j) const {
    return i == j ? 1.0 : 0.0;
}

double MagneticWCSPH::W_avr (const double q) const{
    return w_avr(q/ m_h) / pow(m_h, 3.0);
}

double MagneticWCSPH::wprime (const double q) const {
    if (0 <= q && q < 1) {
        return (-0.75 * pow(2.0-q, 2.0) + 3.0 * pow(1.0-q, 2.0))/M_PI;
    } else if (1 <= q && q < 2) {
        return (-0.75 * pow(2.0-q, 2.0))/M_PI;
    } else {
        return 0.0;
    }
}

double MagneticWCSPH::Wprime (const double q) const  {
    return wprime(q/m_h)/pow(m_h, 4.0);
}

double MagneticWCSPH::A (const double q) const  {
    return (W(q) - W_avr(q))/pow(q, 2.0);
}

double MagneticWCSPH::Aprime (const double q) const  {
    return 5.0 * (W(q) - W_avr(q))/pow(q, 3.0) - Wprime(q)/pow(q, 2.0);
}

double MagneticWCSPH::C1 (const double q) const  {
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
    } else {
        assert(false);
    }
    return qexp.dot(coeff);
}

double MagneticWCSPH::C2 (const double q) const  {
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
    } else {
        assert(false);
    }
    return qexp.dot(coeff);
}

double MagneticWCSPH::getMagneticSusceptibility() const{
    return (double) _chi;
}

VectorXd MagneticWCSPH::getExternalB() const{
    return m_Bext;
}

double MagneticWCSPH::getPermeability() const{
    return (double) _mu_0;

}
double MagneticWCSPH::getGamma() const{
    double V = getVolume();
    double chi = getMagneticSusceptibility();
    return V * chi / (1 + chi);
}

void MagneticWCSPH::draw(Shader *shader) {
    WCSPH::draw(shader);
    static Shape shape = getLineShape();
    if(m_guess.size() == 0) {
        return;
    } else {
        assert(m_guess.size() == _fluid_ptcl_list.size() * 3);
    }
    for(int i = 0; i < _fluid_ptcl_list.size(); i++) {
        Eigen::Vector3f guess(m_magneticForce(3*i), m_magneticForce(3*i+1), m_magneticForce(3*i+2));//Vector3f::Zero();//(m_guess[3 * i], m_guess[3 * i + 1], m_guess[3 * i + 2]);
//        if(guess.norm() < 0.000001)
//            continue;
        auto& ptcl = _fluid_ptcl_list[i];
        Eigen::Affine3f mat = Eigen::Affine3f::Identity();
        Eigen::AngleAxis<float> aa;
        Quaternionf q;
        q = Quaternionf().setFromTwoVectors(Eigen::Vector3f(1, 0, 0), guess);
        aa = q;

        mat.translate(Eigen::Vector3f(ptcl->position[0], ptcl->position[1], ptcl->position[2]));
        mat.rotate(aa);
        mat.scale(Eigen::Vector3f(0.1, 0.1, 0.1));

        shape.setModelMatrix(mat);
        shape.draw(shader);
    }
}
