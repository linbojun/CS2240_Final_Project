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
    Eigen::Vector3f position;
    Eigen::Vector3f normal;
    float density;
    float pressure;
    float mass;
    // std::vector<int> _boundaryActive;
    bool active;
    std::vector<std::shared_ptr<fluid_ptcl>> fluid_neighs;
    std::vector<std::shared_ptr<wall_ptcl>> wall_neighs;
    
} wall_ptcl;


typedef struct fluid_ptcl{
    Eigen::Vector3f position;
    Eigen::Vector3f velocity;
    // Vector3f positionsNew;
    // Vector3f velocitiesNew;
    // Vector3f positionsPreShock;
    // Vector3f velocitiesPreShock;
    Eigen::Vector3f normal;
    Eigen::Vector3f netForce;
    // Vector3f pressureForces;
    float density;
    float pressure;
//    float drhodt;
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
    void update(float time_step);




private:
	float fluid_ptcl_mass;
    float dt = 0.005f;
	float ptcl_radius = 0.01;
	float rho0 = 1000;
	float surface_tension = 1;
//    float alpha = 0.1;
    float alpha = 1;


	float kernel_factor = 3.0; //kernel_radius = kernel_factor * ptcl_radius 
	float kernel_radius; //h, smoothing length 
	float compression_threshold;

    PositionInit m_posInit;
    Kernel kernel;

    std::vector<std::shared_ptr<wall_ptcl>> _wall_ptcl_list;
    std::vector<std::shared_ptr<fluid_ptcl>> _fluid_ptcl_list;


	// Box3f bounds;


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
//Shape get_sphere_shape(float r, int res);





#endif
