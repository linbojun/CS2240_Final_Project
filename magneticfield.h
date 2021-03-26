#include "Eigen/Dense"
#include "Eigen/Sparse"
#include "Eigen/SparseQR"

typedef Eigen::SparseMatrix<double> SpMat; // declares a column-major sparse matrix type of double
typedef Eigen::Triplet<double> T;

class SPH 
{
public:
    Eigen::Vector3f getPos(int particle) const ;

    std::vector<int>& getParticles() const;
    std::vector<int>& getNeighbors(int particle, float dist) const;
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

    void buildProblem(Eigen::MatrixXd &mat,const SPH &particles);
    float getG (const SPH &particles, int i, int j, int k, int l);

};

    