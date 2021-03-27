#ifndef SPH_H
#define SPH_H

#include <string>
#include <Eigen/Dense>
#include <Eigen/StdVector>
#include <vector>

#define _particle_radius 0.02f
#define _neighbor_radius _particle_radius*1.3f
#define _rho0 1000
#define _dt 0.0001
#define _mu 1
#define _gamma 7.0f
#define _c 9.0f
#define _dh _neighbor_radius

//#define dot(a,b) a.dot(b)
#define _W(rab, h)   315.0f/(64.0f * M_PI * pow(h, 9)) * pow((h*h -  pow(rab.norm(), 2)),3)
#define _delta_W(rab, h)  -45.0f/(M_PI * pow(h, 6)) * (h - rab.norm()) * (h - rab.norm()) * rab / rab.norm()
#define _delta_square_W(rab, h) 45.0f/(M_PI * pow(h, 6)) * (h-rab.norm())




using namespace Eigen;
using namespace std;

typedef struct particle{
    Vector3d position;
//    Vector3d mat_position;
    Vector3d velocity;
    double pressure;
    double density;
    float mass;
    vector<shared_ptr<particle>> neighs;
    Vector3d dvdt;
    double drhodt;


}particle;

class SPH{
public:
	SPH(int n);
	void update(double time_step);


private:
	int size;
	double m_time;
	vector<shared_ptr<particle>> m_particle_list;

//    Vector3d single_Viscousity(shared_ptr<particle> cur);
//    Vector3d single_Momentum(shared_ptr<particle> cur);
//    Vector3d single_dvdt(shared_ptr<particle> cur);
//    double single_pressure(shared_ptr<particle> cur);
//    double single_density(shared_ptr<particle> cur);


    void euler_step();
    vector<shared_ptr<particle>> find_neighs(shared_ptr<particle> cur);
    Vector3d total_dvdt(shared_ptr<particle> cur);
    Vector3d momentum_dvdt(shared_ptr<particle> cur);
    Vector3d viscosity_dvdt(shared_ptr<particle> cur);
    double single_drhodt(shared_ptr<particle> cur);
    double single_pressure(shared_ptr<particle> cur);
    void boundry_collision();








};

static inline double W(double r, double h);
static inline double dW(double r, double h);
static inline Vector3d gradW(Vector3d r, double h);


#endif
