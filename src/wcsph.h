#ifndef WCSPH_H
#define WCSPH_H


#include <string>
#include <Eigen/Dense>
#include <Eigen/StdVector>
#include <vector>
#include <graphics/shape.h>
#include "Kernel.h"
#include "positioninit.h"
#define _USE_MATH_DEFINES
#include <math.h>

//using namespace Eigen;
//using namespace std;

typedef struct fluid_ptcl fluid_ptcl ;

typedef struct wall_ptcl
{
    Eigen::Vector3d position;
    Eigen::Vector3d normal;
    double density;
    double pressure;
    double mass;
    // std::vector<int> _boundaryActive;
    bool active;
    std::vector<int> fluid_neighs;
    std::vector<std::shared_ptr<wall_ptcl>> wall_neighs;
    
} wall_ptcl;


typedef struct fluid_ptcl{
    Eigen::Vector3d position;
    Eigen::Vector3d velocity;
    // Vector3d positionsNew;
    // Vector3d velocitiesNew;
    // Vector3d positionsPreShock;
    // Vector3d velocitiesPreShock;
    Eigen::Vector3d normal;
    Eigen::Vector3d netForce;
    // Vector3d pressureForces;
    double density;
    double pressure;
//    double drhodt;
    std::vector<int> fluid_neighs;
    std::vector<std::shared_ptr<wall_ptcl>> wall_neighs;
    bool shouldSim;
}fluid_ptcl;


#define _GAMMA 7.0
#define _VISCOSITY 0.01f
#define _C 9.0f


class WCSPH{

public:
    WCSPH(double radius, double spacing);
    virtual void draw(Shader *shader);
    virtual void update(double time_step);

    int getNumParticle() { return _fluid_ptcl_list.size(); }
    const std::vector<std::shared_ptr<fluid_ptcl>> &getParticles() { return _fluid_ptcl_list; }

    Eigen::Vector3d getPos(int particle) { return _fluid_ptcl_list[particle]->position; }
    double getVolume() const { return pow(kernel_radius, 3.0); }
protected:
    double fluid_ptcl_mass;
    double dt = 0.003;
    double ptcl_radius;
    double rho0 = 300;
    double surface_tension = 5;
    int _num_part_sim;
    float simwait_secs = 10;
    double t = 0;
    double alpha = 0.2;
//    double alpha = 1;
    int _grid_segs;
    float _voxel_len;
    int _max_grid_search;
    std::vector<std::vector<int>> m_grid;



    double kernel_factor = 3.0; //kernel_radius = kernel_factor * ptcl_radius
    double kernel_radius; //h, smoothing length
    double compression_threshold;

    PositionInit m_posInit;
    Kernel kernel;
    std::vector<std::shared_ptr<wall_ptcl>> _wall_ptcl_list;
    std::vector<std::shared_ptr<fluid_ptcl>> _fluid_ptcl_list;

    // Box3d bounds;


    void update_all_neighs();
    void find_fluid_neighs(int pi);
    void find_wall_neighs(int pi);

    void update_all_density_and_pressure();
    void update_all_density_and_pressure_old();

    void update_all_normal();
    void update_net_force();
    virtual void update_velocity_position();
    void boundary_collision();
    void single_pressure(std::shared_ptr<fluid_ptcl> cur);
    void single_drhodt(std::shared_ptr<fluid_ptcl> cur);
    Eigen::Vector3i gridPlace(const Eigen::Vector3d& pos);
    void updateParticlePos(int i, Eigen::Vector3d newPos, bool initializing=false);

    int gridPlaceIndex(const Eigen::Vector3i& place);

    int gridIndex(Eigen::Vector3d& pos);



    Shape get_sphere_shape(float r, int res);


};
//Shape get_sphere_shape(double r, int res);





#endif
