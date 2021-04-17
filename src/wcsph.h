#ifndef WCSPH_H
#define WCSPH_H


#include <string>
#include <Eigen/Dense>
#include <Eigen/StdVector>
#include <vector>
#include <graphics/shape.h>
#include "Kernel.h"
#include "positioninit.h"

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
    std::vector<std::shared_ptr<fluid_ptcl>> fluid_neighs;
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
    std::vector<std::shared_ptr<fluid_ptcl>> fluid_neighs;
    std::vector<std::shared_ptr<wall_ptcl>> wall_neighs;
}fluid_ptcl;


#define _GAMMA 7.0
#define _VISCOSITY 0.05f
#define _C 9.0f


class WCSPH{

public:
    WCSPH();
    void draw(Shader *shader);
    void update(double time_step);




private:
    double fluid_ptcl_mass;
    double dt = 0.0625;
    //double dt = 0.00625;
    double ptcl_radius = 0.01;
    double rho0 = 10000;
    double surface_tension = 2.0;
    double alpha = 0.5;
//    double alpha = 1;


    double kernel_factor = 15.0; //kernel_radius = kernel_factor * ptcl_radius
    double kernel_radius; //h, smoothing length
    double compression_threshold;

    PositionInit m_posInit;
    Kernel kernel;

    std::vector<std::shared_ptr<wall_ptcl>> _wall_ptcl_list;
    std::vector<std::shared_ptr<fluid_ptcl>> _fluid_ptcl_list;


    // Box3d bounds;


    void update_all_neighs();
    void find_fluid_neighs(std::shared_ptr<fluid_ptcl> cur);
    void find_wall_neighs(std::shared_ptr<wall_ptcl> cur);

    void update_all_density_and_pressure();
    void update_all_density_and_pressure_old();

    void update_all_normal();
    void update_net_force();
    void update_velocity_position();
    void boundary_collision();
    void single_pressure(std::shared_ptr<fluid_ptcl> cur);
    void single_drhodt(std::shared_ptr<fluid_ptcl> cur);


    Shape get_sphere_shape(float r, int res);


};
//Shape get_sphere_shape(double r, int res);





#endif
