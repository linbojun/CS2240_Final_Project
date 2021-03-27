#ifndef SPH_H
#define SPH_H

#include <string>
#include <Eigen/Dense>
#include <Eigen/StdVector>

#define _neighbor_radius 1
#define _rho0 1000
#define _particle_radius 0.02f
#define _dt 0.0001
#define _mu 1
#define _gamma 7.0f
#define _c 9.0f

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


}particle;

class SPH{
public:
	SPH(int n);
	void update(double time_step);

private:
	int size;
	double m_time;
	vector<shared_ptr<particle>> m_particle_list;

//    double single_density(shared_ptr<particle> cur);

//    Vector3d single_accerlation_momentum(shared_ptr<particle> cur);

//    void computeDensities();
//    void computePressure();
//    void computeViscousity();
//    void compute_total_force();
    Vector3d single_Viscousity(shared_ptr<particle> cur);
    Vector3d single_Momentum(shared_ptr<particle> cur);
    Vector3d single_dvdt(shared_ptr<particle> cur);
    double single_pressure(shared_ptr<particle> cur);
    double single_density(shared_ptr<particle> cur);
    void euler_step();








};

static inline double W(double r, double h);
static inline double dW(double r, double h);
static inline Vector3d gradW(Vector3d r, double h);


#endif
