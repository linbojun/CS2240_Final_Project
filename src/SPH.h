#ifndef SPH_H
#define SPH_H

#include <string>
#include <Eigen/Dense>
#include <Eigen/StdVector>
#include <vector>
#include <graphics/shape.h>

const float _particle_radius = 0.004f;
const float _neighbor_radius = _particle_radius*1.3f;
const float _rho0 = 1000;
const float _dt = 0.0001;
const float _mu = 1;
const float _gamma = 7.0f;
const float _c = 9.0f;
const float _dh =_neighbor_radius;
// the best number for grid_segs is probably that which makes max_grid_search equal to 1
const int _grid_segs = 1 / _neighbor_radius;
const float _voxel_len = 1.f/_grid_segs;
const int _max_grid_search = ceil(_neighbor_radius / _voxel_len);

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


}particle;

class SPH{
public:
	SPH(int n);
    void update(float time_step);
    void draw(Shader *shader);


private:
	int size;
	double m_time;
    std::vector<std::shared_ptr<particle>> m_particle_list;
    //std::vector<Shape> m_p_shapes;
    Shape m_shape;
    std::vector<std::vector<int>> m_grid;

//    Vector3d single_Viscousity(shared_ptr<particle> cur);
//    Vector3d single_Momentum(shared_ptr<particle> cur);
//    Vector3d single_dvdt(shared_ptr<particle> cur);
//    double single_pressure(shared_ptr<particle> cur);
//    double single_density(shared_ptr<particle> cur);


    void euler_step();
    std::vector<std::shared_ptr<particle>> find_neighs(int pi);
    Eigen::Vector3d total_dvdt(std::shared_ptr<particle> cur);
    Eigen::Vector3d momentum_dvdt(std::shared_ptr<particle> cur);
    Eigen::Vector3d viscosity_dvdt(std::shared_ptr<particle> cur);
    double single_drhodt(std::shared_ptr<particle> cur);
    double single_pressure(std::shared_ptr<particle> cur);
    void boundry_collision();
    void updateParticlePos(int i, Eigen::Vector3d newPos, bool initializing=false);









};

static inline double W(double r, double h);
static inline double dW(double r, double h);
static inline Eigen::Vector3d gradW(Eigen::Vector3d r, double h);


#endif
