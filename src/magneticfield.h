

#include "Eigen/Dense"
#include "Eigen/Sparse"
#include "Eigen/SparseQR"

typedef Eigen::SparseMatrix<double> SpMat; // declares a column-major sparse matrix type of double
typedef Eigen::Triplet<double> T;

class _SPHApi
{
public:
    Eigen::Vector3f getPos(int particle) const ;

    std::vector<int>& getParticles() const;
    std::vector<int>& getNeighbors(int particle, float dist) const;

    double getVolume() const ;
    double getMagneticSusceptibility() const;
    Eigen::VectorXd& getExternalB() const;
};

class MagneticField
{

public:
    
    MagneticField(double h);

private:

    double m_h;

    const double w (const double q) const;
    const double w_avr (const double q) const;
    
    const double W (const double q) const;
    const double W_avr (const double q) const;

    const double delta(const int i,const int j) const;

    Eigen::VectorXd calculateMagneticField(const _SPHApi &particles);

    void buildProblem(Eigen::MatrixXd &mat,const _SPHApi &particles);
    
    
    float getG (const _SPHApi &particles, int i, int j, int k, int l);

};

    
