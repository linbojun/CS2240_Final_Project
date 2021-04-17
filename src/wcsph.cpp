#include "wcsph.h"
#include <math.h>
#include <iostream>


using namespace Eigen;
using namespace std;

#define USE_WALL 0
WCSPH::WCSPH():
    m_posInit(ptcl_radius*2)
{
    kernel_radius = kernel_factor * ptcl_radius;
    fluid_ptcl_mass = 4.0/3.0 * M_PI * ptcl_radius * ptcl_radius * ptcl_radius;
    // bounds.reset();

    int num_fluid_particle = 1000;
    kernel.init(kernel_radius);
    m_posInit.addBox(Vector3f(0.5, 0.4, 0.5), 0.2, 0.2, 0.2);


    //assume bound is x = 1,-1
//    for(int i = 0; i < num_fluid_particle; i++)
//    {
//         shared_ptr<fluid_ptcl> new_particle(new fluid_ptcl);
//         _fluid_ptcl_list.push_back(new_particle);
//         double rand_x = static_cast <double> (rand()) / (static_cast <double> (RAND_MAX/2));
//         double rand_y = static_cast <double> (rand()) / (static_cast <double> (RAND_MAX/2));
//         double rand_z = static_cast <double> (rand()) / (static_cast <double> (RAND_MAX/2));
//         Vector3d pos(rand_x, rand_y, rand_z);
//         Vector3d zeros(0,0,0);
//         new_particle->position = pos;
//         new_particle->velocity = zeros;
//         new_particle->pressure = 0;
//         new_particle->density = rho0;
	
//    }
    double numParticles = m_posInit.getNumParticles();
    cout << "num" << numParticles << endl;


    for(int i = 0; i < m_posInit.getNumParticles(); i++){
         shared_ptr<fluid_ptcl> new_particle(new fluid_ptcl);
        _fluid_ptcl_list.push_back(new_particle);
        Vector3d zeros(0,0,0);
//        updateParticlePos(i, Vector3d::Random() * 0.5 + Vector3d(0.5, 0.5, 0.5));
        //updateParticlePos(i, m_posInit.getPt(i).cast<double> ());
        new_particle->position = m_posInit.getPt(i).cast<double>();
        new_particle->velocity = zeros;
        new_particle->pressure = 0;
        new_particle->density = rho0;
        //m_p_shapes.push_back(get_sphere_shape(adius));
    }

    //add the wall particle


}

void WCSPH::update_all_neighs()
{
#if USE_WALL
    for(auto wall_ptcl : _wall_ptcl_list)
    {
        wall_ptcl->active = false;
    }
#endif
    for(auto cur: _fluid_ptcl_list)
    {
        find_fluid_neighs(cur);
        //cout<<"num_neighs:"<<cur->fluid_neighs.size()<<endl;;
    }
#if USE_WALL
    for(auto cur: _wall_ptcl_list)
    {
        find_wall_neighs(cur);
    }
#endif
}

void WCSPH::find_fluid_neighs(shared_ptr<fluid_ptcl> cur)
{
    cur->fluid_neighs.clear();
    cur->wall_neighs.clear();
    auto ra = cur->position;
    for(auto fluid_ptcl: _fluid_ptcl_list)
    {
        auto rb = fluid_ptcl->position;
        auto rab = ra - rb;
        if(rab.norm() < kernel_radius && rab.norm() > 0)
        {
            cur->fluid_neighs.push_back(fluid_ptcl);
        }
    }

    for(auto wall_ptcl: _wall_ptcl_list)
    {
        auto rb = wall_ptcl->position;
        auto rab = ra - rb;
        if(rab.norm() < kernel_radius)
        {
            cur->wall_neighs.push_back(wall_ptcl);
            wall_ptcl->active = true;
        }
    }
}

void WCSPH::find_wall_neighs(shared_ptr<wall_ptcl> cur)
{
    if(!cur->active){
        return;
    }
    cur->wall_neighs.clear();
    cur->fluid_neighs.clear();
    auto ra = cur->position;
    for(auto fluid_ptcl: _fluid_ptcl_list)
    {
        auto rb = fluid_ptcl->position;
        auto rab = ra - rb;
        if(rab.norm() < kernel_radius)
        {
            cur->fluid_neighs.push_back(fluid_ptcl);
        }
    }

    for(auto wall_ptcl: _wall_ptcl_list)
    {
        auto rb = wall_ptcl->position;
        auto rab = ra - rb;
        if(rab.norm() < kernel_radius && rab.norm() > 0)
        {
            cur->wall_neighs.push_back(wall_ptcl);
        }
    }
}

void WCSPH::update(double time_step)
{
	
    update_all_neighs();
//    update_all_density_and_pressure();
    update_all_density_and_pressure_old();
    update_all_normal();
    update_net_force();
    update_velocity_position();
    boundary_collision();
//    exit(0);

}

void WCSPH::update_all_density_and_pressure_old()
{
    cout<<"update_all_density_and_pressure_old()"<<endl;
    for(auto cur: _fluid_ptcl_list)
    {
        single_drhodt(cur);
        single_pressure(cur);
    }
}

void WCSPH::single_drhodt(shared_ptr<fluid_ptcl> cur)
{
    double drhodt = 0;
    auto ra = cur->position;
    auto va = cur->velocity;
    auto neighs = cur->fluid_neighs;
    for(auto neigh: neighs)
    {
        auto rb = neigh->position;
        auto vb = neigh->velocity;
        auto vab = va - vb;
        auto rab = ra - rb;
        auto rab_sqr = rab.norm() * rab.norm();
//        drhodt += fluid_ptcl_mass * vab.dot(gradW(rab, _dh));
        drhodt += fluid_ptcl_mass * vab.dot(kernel.poly6GradConstant * kernel.poly6Grad(rab, rab_sqr));
    }

//    cout<<"====================================="<<endl;
//    cout<<"density_old:"<<cur->density<<endl;
    cur->density += drhodt * dt;
//    cout<<"drhodt:"<<drhodt<<endl;
//    cout<<"density_new:"<<cur->density<<endl;

}

void WCSPH::single_pressure(shared_ptr<fluid_ptcl> cur)
{
    auto p0 = rho0 * _C * _C / _GAMMA;
    auto press= p0 * (pow(cur->density / rho0, _GAMMA) - 1.0);
//    cout<<endl<<"press:"<<press<<endl;
    cur->pressure = press > 0 ? press : 0;

}

void WCSPH::update_all_density_and_pressure()
{
#if USE_WALL
    for(auto wall_ptcl: _wall_ptcl_list)
    {
        if(wall_ptcl->active == false)
            continue;
        double fluid_density = 0;
        auto ra = wall_ptcl->position;
        for(auto fluid_neigh: wall_ptcl->fluid_neighs)
        {
            auto rb = fluid_neigh->position;
            auto rab = ra - rb;
            double rab_sqr = rab.norm() * rab.norm();
            fluid_density += kernel.poly6(rab_sqr);
        }

        double wall_density = 0;
        for(auto wall_neigh: wall_ptcl->wall_neighs)
        {
            auto rb = wall_neigh->position;
            auto rab = ra - rb;
            double rab_sqr = rab.norm() * rab.norm();
            wall_density += kernel.poly6(rab_sqr) * wall_neigh->mass;
        }

        double density = kernel.poly6Constant * fluid_ptcl_mass * fluid_density;
        density += kernel.poly6Constant * wall_density;

        auto B = rho0 * _C * _C / _GAMMA;
        double token = density / rho0;
        auto pressure = B * (pow(token, _GAMMA) - 1.0);
        wall_ptcl->density = density;
        wall_ptcl->pressure = pressure;

    }
#endif

    for(auto fluid_ptcl: _fluid_ptcl_list)
    {
        double fluid_density = 0.f;
        auto ra = fluid_ptcl->position;
        for(auto fluid_neigh: fluid_ptcl->fluid_neighs)
        {
            auto rb = fluid_neigh->position;
            auto rab = ra - rb;
            double rab_sqr = rab.norm() * rab.norm();
            fluid_density += kernel.poly6(rab_sqr);
        }
        double wall_density = 0;
#if USE_WALL
        for(auto wall_neigh: fluid_ptcl->wall_neighs)
        {
            auto rb = wall_neigh->position;
            auto rab = ra - rb;
            double rab_sqr = rab.norm() * rab.norm();
            wall_density += kernel.poly6(rab_sqr) * wall_neigh->mass;
        }
#endif

        double density = kernel.poly6Constant * fluid_ptcl_mass * fluid_density;
        density += kernel.poly6Constant * wall_density;
//        cout<<"=================================="<<endl;
//        cout<<"kernel.poly6Constant"<<kernel.poly6Constant<<endl;
//        cout<<"fluid_ptcl_mass"<<fluid_ptcl_mass<<endl;
//        cout<<"fluid_density: "<<fluid_density<<endl;
//        cout<<"density: "<<fluid_density<<endl;



        auto B = rho0 * _C * _C / _GAMMA;
        double token = density / rho0;
        auto pressure = B * (pow(token, _GAMMA) - 1.0);
        fluid_ptcl->density = density;
        fluid_ptcl->pressure = pressure;
//        cout<<endl;
//        cout<<"B:"<<B<<endl;
//        cout<<"density / rho0:"<<density / rho0<<endl;
//        cout<<"pressure:"<<pressure<<endl;

    }
}

void WCSPH::update_all_normal()
{
    for(auto fluid_ptcl: _fluid_ptcl_list)
    {
//        cout<<"+++++++++++++++++++++++++++++"<<endl;
        Vector3d normal(0, 0, 0);
        auto ra = fluid_ptcl->position;
        for(auto fluid_neigh: fluid_ptcl->fluid_neighs)
        {
            auto rb = fluid_neigh->position;
            auto rho_b = fluid_neigh->density;
            auto rab = ra - rb;
            double rab_sqr = rab.norm() * rab.norm();
            normal += kernel.poly6GradConstant * kernel.poly6Grad(rab, rab_sqr) / rho_b ;
//            cout<<"kernel.poly6GradConstant:"<<kernel.poly6GradConstant<<endl;
//            cout<<"kernel.poly6Grad(rab, rab_sqr):"<<kernel.poly6Grad(rab, rab_sqr)<<endl;
//            cout<<"rho_b:"<<rho_b<<endl;
//            cout<<"product:"<<kernel.poly6GradConstant * kernel.poly6Grad(rab, rab_sqr) / rho_b<<endl;

        }
        normal *= kernel_radius * fluid_ptcl_mass;
        fluid_ptcl->normal = normal;

//        cout<<"normal:"<<normal<<endl;
    }
}


void WCSPH::update_net_force()
{
    for(auto cur_fluid: _fluid_ptcl_list)
    {
        Vector3d net_force(0,0,0);
        Vector3d forcePressure(0,0,0);
        Vector3d forceViscosity(0,0,0);
        Vector3d forceCohesion(0,0,0);
        Vector3d forceCurvature(0,0,0);

        auto va = cur_fluid->velocity;
        auto ra = cur_fluid->position;
        auto na = cur_fluid->normal;
        auto rho_a = cur_fluid->density;
        auto pa = cur_fluid->pressure;

        //force by the fluid particle
        for(auto fluid_neigh: cur_fluid->fluid_neighs)
        {
            auto vb = fluid_neigh->velocity;
            auto rb = fluid_neigh->position;
            auto nb = fluid_neigh->normal;
            auto rho_b = fluid_neigh->density;
            auto pb = fluid_neigh->pressure;
            if(cur_fluid  == fluid_neigh)
                continue;
            auto rab = ra - rb;
            auto vab = va-vb;
            double rab_sqr = rab.norm() * rab.norm();
            if(rab.norm() < kernel_radius && rab.norm() > 0.001f)
            {
                //Momentum
                //forcePressure -= fluid_ptcl_mass * fluid_ptcl_mass  * (pressure_a / (rho_a * rho_a) + pressure_b / (rho_b * rho_b)) * kernel.spikyGradConstant * kernel.spikyGrad(rab, rab.norm());
                forcePressure -= fluid_ptcl_mass * fluid_ptcl_mass * (pa / (rho_a * rho_a) + pb / (rho_b * rho_b)) * kernel.spikyGradConstant * kernel.spikyGrad(rab, rab.norm());


                //Viscosity
                if(vab.dot(rab) < 0.0)
                {
                    double eps = 0.01f;
                    auto v = 2 * alpha * kernel_radius * _C / (rho_a + rho_b);
                    auto pi_ab = -v * vab.dot(rab) / (rab.dot(rab) + eps * kernel_radius * kernel_radius);
                    forceViscosity -= fluid_ptcl_mass *fluid_ptcl_mass * pi_ab * kernel.spikyGradConstant * kernel.spikyGrad(rab, rab.norm()); // minus?
                    //dvdt += mass * pi_ab * gradW(rab, dh);
                }
//                    forceViscosity -= (va - vb) * (kernel.viscosityLaplace(rab.norm()) / rho_b);

                //surface tension
                double correctionFactor = 2.f * rho0 / (rho_a + rho_b);
                forceCohesion += correctionFactor * (rab / rab.norm()) * kernel.surfaceTension(rab.norm());
                forceCurvature += correctionFactor * (na - nb);
            }
            else if (rab_sqr == 0.f)
            {
                    // Avoid collapsing particles
                    fluid_neigh->position += Vector3d(1e-5f, 1e-5f, 1e-5f);
            }


        }
#if USE_WALL
        //force by wall particle
        for(auto wall_neigh: cur_fluid->wall_neighs)
        {
            auto rb = wall_neigh->position;
            auto rho_b = wall_neigh->density;
            auto pressure_b = wall_neigh->pressure;
            auto rab = ra - rb;
            if (rab.norm() < kernel_radius && rab.norm() > 0.001f)
            {
                forcePressure -= fluid_ptcl_mass * wall_neigh->mass * (pressure_a / (rho_a * rho_a) + pressure_b / (rho_b * rho_b)) * kernel.spikyGradConstant * kernel.spikyGrad(rab, rab.norm());
            }


        }
#endif

//        forceViscosity *= _VISCOSITY * fluid_ptcl_mass * kernel.viscosityLaplaceConstant;

        forceCohesion *= -surface_tension * fluid_ptcl_mass * fluid_ptcl_mass * kernel.surfaceTensionConstant;
        forceCurvature *= -surface_tension * fluid_ptcl_mass;

        net_force += forcePressure;
        net_force += forceViscosity;
        net_force += forceCohesion + forceCurvature;
//        cout<<"-------------------------------------"<<endl;
//        cout<<"forcePressure:"<<forcePressure<<endl;
//        cout<<"forceViscosity:"<<forceViscosity<<endl;
//        cout<<"forceCohesion:"<<forceCohesion<<endl;
//        cout<<"forceCurvature:"<<forceCurvature<<endl;
        //force += fluid_ptcl_mass * _gravity;
        cur_fluid->netForce = net_force;
    }

}

void WCSPH::update_velocity_position()
{
    for(auto cur_fluid: _fluid_ptcl_list)
    {
        Vector3d dvdt = cur_fluid->netForce / fluid_ptcl_mass;
        cur_fluid->velocity += dvdt * dt;
        cur_fluid->position += cur_fluid->velocity * dt;
    }
}

void WCSPH::boundary_collision()
{
    for(auto cur_fluid: _fluid_ptcl_list)
    {
        auto r = cur_fluid->position;
        auto v = cur_fluid->velocity;
        auto k = 0.3;
        for (int i = 0; i < 3; i++)
        {
            if (r(i,0) < 0)
            {
                cur_fluid->position(i,0) = 0;
                if (v(i,0) < 0.0) {

//                    cur_fluid->velocity(i,0) += -(1 + k) * v(i,0);
                     cur_fluid->velocity(i,0)  = 0;
                }
            }
            if (r(i,0) > 1)
            {
                cur_fluid->position(i,0) = 1;
                if (v(i,0) > 0.0) {
//                    cur_fluid->velocity(i,0) += -(1 + k) * v(i,0);
                    cur_fluid->velocity(i,0) = 0;
                }
            }
        }
    }

}

void WCSPH::draw(Shader *shader) {
    static Shape shape = get_sphere_shape(ptcl_radius, 10);
    for(int i = 0; i < _fluid_ptcl_list.size(); i++) {
        auto& ptcl = _fluid_ptcl_list[i];
        shape.setModelMatrix(Eigen::Affine3f(Eigen::Translation3f(ptcl->position[0], ptcl->position[1], ptcl->position[2])));
        shape.draw(shader);
    }
}

Shape WCSPH::get_sphere_shape(float r, int res) {
    std::vector<Vector3f> points;
    std::vector<Vector3i> faces;
    for(int thetas = 0; thetas < res; thetas++) {
        for(int phis = 1; phis < res-1; phis++) {
            float x = r * sin(phis*M_PI/res)*cos(thetas*2*M_PI/res);
            float y = r * sin(phis*M_PI/res)*sin(thetas*2*M_PI/res);
            float z = r * cos(phis*M_PI/res);
            points.push_back(Vector3f(x, y, z));
        }
    }
    // idx = thetas*18 + phis - 1
    for(int thetas = 1; thetas < res+1; thetas++) {
        for(int phis = 2; phis < res; phis++) {
            faces.push_back(Vector3i((thetas%res)*(res-2) + phis - 1, (thetas-1)*(res-2) + phis - 1, (thetas-1)*(res-2) + phis - 2));
            faces.push_back(Vector3i((thetas%res)*(res-2) + phis - 1, (thetas-1)*(res-2) + phis - 2, (thetas%res)*(res-2) + phis - 2));
        }
    }
    int topi = points.size();
    points.push_back(Vector3f(0, 0, r));
    int boti = points.size();
    points.push_back(Vector3f(0, 0, -r));
    for(int thetas = 0; thetas < res; thetas++) {
        faces.push_back(Vector3i(topi, thetas*(res-2), ((thetas+1)%res)*(res-2)));
        faces.push_back(Vector3i(boti, thetas*(res-2) + (res-3), ((thetas+1)%res)*(res-2) + (res-3)));
    }
    Shape shape;
    shape.init(points, faces, true);
    shape.setVertices(points);
    return shape;
}

