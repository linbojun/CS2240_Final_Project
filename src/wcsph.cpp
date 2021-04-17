#include "wcsph.h"
#include <math.h>
#include <iostream>


using namespace Eigen;
using namespace std;

inline bool inBounds(const Vector3d& pos) {
    return 1 >= pos[0] && pos[0] >= 0 && 1 >= pos[1] && pos[1] >= 0 && 1 >= pos[2] && pos[2] >= 0;
}

inline Vector3i WCSPH::gridPlace(const Vector3d& pos) {
    return Vector3i(pos[0] * (_grid_segs - 1), pos[1] * (_grid_segs - 1), pos[2] * (_grid_segs - 1));
}

inline int WCSPH::gridPlaceIndex(const Vector3i& place) {
    return place[0] * _grid_segs * _grid_segs + place[1] * _grid_segs + place[2];
}

inline int WCSPH::gridIndex(Vector3d& pos) {
    return gridPlaceIndex(gridPlace(pos));
}

void WCSPH::updateParticlePos(int i, Vector3d newPos, bool initializing) {
    auto particle = _fluid_ptcl_list[i];
    if(!initializing) {
        // delete from old spot in grid
        if(inBounds(particle->position)) {
            auto& vec = m_grid[gridIndex(particle->position)];
            vec.erase(std::remove(vec.begin(), vec.end(), i), vec.end());
        }
    }
    // add to new spot in grid
    if(inBounds(newPos)) {
        m_grid[gridIndex(newPos)].push_back(i);
    }
    particle->position = newPos;
}


#define USE_WALL 0
WCSPH::WCSPH():
    m_posInit(ptcl_radius*1.8)
{

    kernel_radius = kernel_factor * ptcl_radius;
    fluid_ptcl_mass = 4.0/3.0 * M_PI * ptcl_radius * ptcl_radius * ptcl_radius;

    _grid_segs = 1 / kernel_radius;
    _voxel_len = 1.f/_grid_segs;
    _max_grid_search = ceil(kernel_radius / _voxel_len);
    assert(_max_grid_search == 1);

    m_grid.resize(_grid_segs * _grid_segs * _grid_segs);
    // bounds.reset();

    kernel.init(kernel_radius);
    m_posInit.addBox(Vector3f(0.5, 0.3, 0.5), 1, 0.6, 1);
    int npBox = m_posInit.getNumParticles();
    m_posInit.setRadius(ptcl_radius*1.5);
    m_posInit.addSphere(Vector3f(0.5, 0.8, 0.5), 0.08, Vector3f(0, 0, 0));

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
        updateParticlePos(i, m_posInit.getPt(i).cast<double>(), true);
        new_particle->velocity = m_posInit.getVel(i).cast<double>();
        new_particle->pressure = 0;
        new_particle->density = rho0;
        new_particle->shouldSim = i < npBox;
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
    for(int i = 0; i < _fluid_ptcl_list.size(); i++)
    {
        //if(!_fluid_ptcl_list[i]->shouldSim)
        //    continue;
        find_fluid_neighs(i);
        //cout<<"num_neighs:"<<cur->fluid_neighs.size()<<endl;;
    }
#if USE_WALL
    for(int i = 0; i < _wall_ptcl_list.size(); i++)
    {
        find_wall_neighs(i);
    }
#endif
}

void WCSPH::find_fluid_neighs(int pi)
{
    shared_ptr<fluid_ptcl> &cur = _fluid_ptcl_list.at(pi);
    cur->fluid_neighs.clear();
//    assert(inBounds(cur->position));
    Vector3i place = gridPlace(cur->position);
    for(int x = max(0, place[0] - _max_grid_search); x <= min(_grid_segs - 1, place[0] + _max_grid_search); x++) {
        for(int y = max(0, place[1] - _max_grid_search); y <= min(_grid_segs - 1, place[1] + _max_grid_search); y++) {
            for(int z = max(0, place[2] - _max_grid_search); z <= min(_grid_segs - 1, place[2] + _max_grid_search); z++) {
                int idx = gridPlaceIndex(Vector3i(x, y, z));
                for(int poss_neigh : m_grid[idx]) {
                    if(poss_neigh == pi)
                        continue;
                    shared_ptr<fluid_ptcl>& ptcl_other = _fluid_ptcl_list.at(poss_neigh);
                    auto dist = ptcl_other->position - cur->position;
                    if(dist.norm() < kernel_radius) {
                        cur->fluid_neighs.push_back(poss_neigh);
                    }
                }
            }
        }
    }

    cur->wall_neighs.clear();
    auto ra = cur->position;

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

void WCSPH::find_wall_neighs(int pi)
{
    shared_ptr<wall_ptcl> &cur = _wall_ptcl_list.at(pi);
    if(!cur->active){
        return;
    }
    cur->wall_neighs.clear();
    cur->fluid_neighs.clear();
//    assert(inBounds(cur->position));
    Vector3i place = gridPlace(cur->position);
    for(int x = max(0, place[0] - _max_grid_search); x <= min(_grid_segs - 1, place[0] + _max_grid_search); x++) {
        for(int y = max(0, place[1] - _max_grid_search); y <= min(_grid_segs - 1, place[1] + _max_grid_search); y++) {
            for(int z = max(0, place[2] - _max_grid_search); z <= min(_grid_segs - 1, place[2] + _max_grid_search); z++) {
                int idx = gridPlaceIndex(Vector3i(x, y, z));
                for(int poss_neigh : m_grid[idx]) {
                    if(poss_neigh == pi)
                        continue;
                    shared_ptr<fluid_ptcl>& ptcl_other = _fluid_ptcl_list.at(poss_neigh);
                    auto dist = ptcl_other->position - cur->position;
                    if(dist.norm() < kernel_radius) {
                        cur->fluid_neighs.push_back(poss_neigh);
                    }
                }
            }
        }
    }
    auto ra = cur->position;
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
    t += dt;
    if(t > simwait_secs) {
        for(auto cur : _fluid_ptcl_list) {
            if(!cur->shouldSim) {
                cur->shouldSim = true;
                cur->velocity += Eigen::Vector3d(0, -1.5, 0);
            }
        }
    }
//    exit(0);

}

void WCSPH::update_all_density_and_pressure_old()
{
    cout<<"update_all_density_and_pressure_old()"<<endl;
    for(auto cur: _fluid_ptcl_list)
    {
        //if(!cur->shouldSim)
        //    continue;
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
    for(auto i: neighs)
    {
        auto neigh = _fluid_ptcl_list[i];
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
        for(auto i: wall_ptcl->fluid_neighs)
        {
            auto fluid_neigh = _fluid_ptcl_list[i];
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
        //if(!fluid_ptcl->shouldSim)
        //    continue;
        double fluid_density = 0.f;
        auto ra = fluid_ptcl->position;
        for(auto i: fluid_ptcl->fluid_neighs)
        {
            auto fluid_neigh = _fluid_ptcl_list[i];
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
        //if(!fluid_ptcl->shouldSim)
        //   continue;
//        cout<<"+++++++++++++++++++++++++++++"<<endl;
        Vector3d normal(0, 0, 0);
        auto ra = fluid_ptcl->position;
        for(auto i: fluid_ptcl->fluid_neighs)
        {
            auto fluid_neigh = _fluid_ptcl_list[i];
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
        //if(!cur_fluid->shouldSim)
        //    continue;
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
        for(auto i: cur_fluid->fluid_neighs)
        {
            auto fluid_neigh = _fluid_ptcl_list[i];
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
                    updateParticlePos(i, fluid_neigh->position + Vector3d(1e-5f, 1e-5f, 1e-5f));
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
        if(cur_fluid->shouldSim)
            net_force += fluid_ptcl_mass * Vector3d(0, -0.1, 0); // gravity
        cur_fluid->netForce = net_force;
    }

}

void WCSPH::update_velocity_position()
{
    for(int i = 0; i < _fluid_ptcl_list.size(); i++)
    {
        auto cur_fluid = _fluid_ptcl_list[i];
        //if(!cur_fluid->shouldSim)
        //    continue;
        Vector3d dvdt = cur_fluid->netForce / fluid_ptcl_mass;
        cur_fluid->velocity += dvdt * dt;
        updateParticlePos(i, cur_fluid->position + cur_fluid->velocity * dt);
    }
}

void WCSPH::boundary_collision()
{
    for(int i = 0; i < _fluid_ptcl_list.size(); i++)
    {
        auto cur_fluid = _fluid_ptcl_list[i];
        //if(!cur_fluid->shouldSim)
        //    continue;
        auto r = cur_fluid->position;
        auto v = cur_fluid->velocity;
        auto k = 0.4;
        for (int i = 0; i < 3; i++)
        {
            if (r(i,0) < 0)
            {
                r(i,0) = 0;
                if (v(i,0) < 0.0) {

//                    cur_fluid->velocity(i,0) += -(1 + k) * v(i,0);
                     cur_fluid->velocity(i,0)  = -k*cur_fluid->velocity(i,0);
                }
            }
            if (r(i,0) > 1)
            {
                r(i,0) = 1;
                if (v(i,0) > 0.0) {
//                    cur_fluid->velocity(i,0) += -(1 + k) * v(i,0);
                    cur_fluid->velocity(i,0) = -k*cur_fluid->velocity(i, 0);
                }
            }
        }
        updateParticlePos(i, r);
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

