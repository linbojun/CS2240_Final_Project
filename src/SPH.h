#ifndef SPH_H
#define SPH_H

#include <string>
#include <Eigen/Dense>
#include <Eigen/StdVector>
#include <vector>
#include <graphics/shape.h>
#include "positioninit.h"
#define _USE_MATH_DEFINES
#include <math.h>

const float _rho0 = 1000;
const float _dt = 0.0001;
const float _mu = 5;
const float _gamma = 7.0f;
const float _c = 9.0f;

// the best number for grid_segs is probably that which makes max_grid_search equal to 1


//#define dot(a,b) a.dot(b)
#define _W(rab, h)   315.0f/(64.0f * M_PI * pow(h, 9)) * pow((h*h -  pow(rab.norm(), 2)),3)
#define _delta_W(rab, h)  -45.0f/(M_PI * pow(h, 6)) * (h - rab.norm()) * (h - rab.norm()) * rab / rab.norm()
#define _delta_square_W(rab, h) 45.0f/(M_PI * pow(h, 6)) * (h-rab.norm())




//using namespace Eigen;
//using namespace std;

typedef struct particle{
    Eigen::Vector3d position;
//    Vector3d mat_position;
    Eigen::Vector3d velocity;
    double pressure;
    double density;
    float mass;
    std::vector<std::shared_ptr<particle>> neighs;
    Eigen::Vector3d dvdt;
    double drhodt;
    Eigen::Vector3d normal;
} particle;

class SPH{
public:
    SPH(int n, float radius);
    virtual void update(float time_step);
    virtual void draw(Shader *shader);
    const std::vector<std::shared_ptr<particle>> &getParticles() { return m_particle_list; }
    int getNumParticle() const {return m_particle_list.size(); }

    Eigen::Vector3d getPos(int particle) const;
    double getVolume() const;

protected:
    std::vector<std::shared_ptr<particle>> m_particle_list;
    std::vector<std::shared_ptr<particle>> find_neighs(int pi);
    Eigen::Vector3d total_dvdt(std::shared_ptr<particle> cur);
    double single_drhodt(std::shared_ptr<particle> cur);
    void updateParticlePos(int i, Eigen::Vector3d newPos, bool initializing=false);
    double single_pressure(std::shared_ptr<particle> cur);
    void boundry_collision();

private:
	int size;
    int m_numParticles;
    float m_radius;
    float _neighbor_radius;
    int _grid_segs;
    float _voxel_len;
    int _max_grid_search;
    float _dh;

    PositionInit m_posInit;

    double m_time;
    //std::vector<Shape> m_p_shapes;
    Shape m_shape;
    std::vector<std::vector<int>> m_grid;

//    Vector3d single_Viscousity(shared_ptr<particle> cur);
//    Vector3d single_Momentum(shared_ptr<particle> cur);
//    Vector3d single_dvdt(shared_ptr<particle> cur);
//    double single_pressure(shared_ptr<particle> cur);
//    double single_density(shared_ptr<particle> cur);


    void euler_step();
    Eigen::Vector3d momentum_dvdt(std::shared_ptr<particle> cur);
    Eigen::Vector3d viscosity_dvdt(std::shared_ptr<particle> cur);
    Eigen::Vector3d tension_dvdt(std::shared_ptr<particle> cur);
    Eigen::Vector3d adhesion_dvdt(std::shared_ptr<particle> cur);

    Eigen::Vector3i gridPlace(const Eigen::Vector3d& pos);

    Eigen::Vector3d single_normal(std::shared_ptr<particle> cur);



    int gridPlaceIndex(const Eigen::Vector3i& place);

    int gridIndex(Eigen::Vector3d& pos);

};

static inline double W(double r, double h);
static inline double dW(double r, double h);
static inline Eigen::Vector3d gradW(Eigen::Vector3d r, double h);
static inline double W_poly6(Eigen::Vector3d r, double h)
{

  double coefficient = 315.0/(64.0*M_PI*pow(h,9));
  double r_squared = r.dot(r);

  return coefficient * pow(h*h-r_squared, 3);
}

static inline Eigen::Vector3d gradW_poly6(Eigen::Vector3d r, double h)
{

  static double coefficient = -945.0/(32.0*M_PI*pow(h,9));
  double r_squared = r.dot(r);
  return coefficient * pow(h*h - r_squared, 2) * r;
}

static inline double LaplacianW_poly6(Eigen::Vector3d r, double h)
{
  static double coefficient = -945.0/(32.0*M_PI*pow(h,9));
  double r_squared = r.dot(r);

  return coefficient * (h*h-r_squared) * (3.0*h*h - 7.0*r_squared);
}

static inline Eigen::Vector3d gradW_spiky(Eigen::Vector3d r, double h)
{

  static double coefficient = -45.0/(M_PI*pow(h,6));
  double radius = r.norm();
  return coefficient * pow(h-radius, 2) * r / radius;
}


static inline double LaplacianW_viscosity(Eigen::Vector3d r, double h)
{
  static double coefficient = 45.0/(M_PI*pow(h,6));
  double radius = r.norm();
  return coefficient * (h - radius);
}

static inline double W_cohesion(Eigen::Vector3d r, double h)
{
    double coefficient = 32.0/(M_PI * pow(h, 9));
    double ans = 0;
    if(2*r.norm() > h && r.norm() <= h)
    {
        ans = coefficient * pow(h - r.norm(), 3) * pow(r.norm(), 3);
    }
    else if(r.norm() >0 && 2*r.norm() <= h)
    {
        ans = coefficient *2.0 * pow(h - r.norm(), 3) * pow(r.norm(), 3); - pow(h, 6)/64.0;
    }
    return ans;
}

static inline double W_adhesion(Eigen::Vector3d r, double h)
{
    double coefficient = 0.007/pow(h, 3.25);
    double ans = 0;
    if(2*r.norm() > h && r.norm() <= h)
    {
        auto token = -4*r.norm() *r.norm()/h +6*r.norm()-2*h;
        ans = coefficient * pow( token, 0.25);
    }
    return ans;

}



#endif
