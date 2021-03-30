

#include "Eigen/Dense"
#include "Eigen/Sparse"
#include "Eigen/SparseQR"


class _SPHApi
{
private:
     Eigen::Vector3d m_pos = Eigen::Vector3d(0.0, 0.0, 0.0);
     int m_num_particles = 0;
     Eigen::VectorXd m_B = Eigen::VectorXd::Zero(10);
public:
    _SPHApi(){
     }

    Eigen::Vector3d getPos(int particle) const {
        return m_pos;
    };

    int getParticles() const{
        return m_num_particles;
    };

    double getVolume() const {
        return 0.0;
    };
    double getMagneticSusceptibility() const {
        return 0.0;
    };
    Eigen::VectorXd getExternalB() const{
        return m_B;
    };
};

class MagneticField
{

public:
    
    MagneticField(double h);

private:

    double m_h;

    Eigen::VectorXd m_guess;
    bool m_isFirst = true;

    Eigen::ConjugateGradient<Eigen::Matrix4d, Eigen::Lower| Eigen::Upper> m_cg;

    double w (const double q) const;
    double w_avr (const double q) const;
    
    double W (const double q) const;
    double W_avr (const double q) const;

    double delta(const int i,const int j) const;

    Eigen::VectorXd calculateMagneticField(const _SPHApi &particles);

    void buildProblem(Eigen::MatrixXd &mat,const _SPHApi &particles);
    
    
    float getG (const _SPHApi &particles, int i, int j, int k, int l);

};

    
