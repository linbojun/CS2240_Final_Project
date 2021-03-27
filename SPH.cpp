#include <SPH.h>
#include <math.h>

using namespace Eigen;
using namespace std;

static inline double _W(double r, double h) {
    auto k = 10. / (7. * M_PI * h * h);
    auto q = r / h;
    auto res = 0.0;
    if (q <= 1.0)
        res = k * (1 - 1.5 * q * q + 0.75 * q * q * q);
    else if (q < 2.0) {
        auto two_m_q = 2 - q;
        res = k * 0.25 * two_m_q * two_m_q * two_m_q;
    }
    return res;
}
static inline double dW(double r, double h) {
    auto k = 10. / (7. * M_PI * h * h);
    auto q = r / h;
    auto res = 0.0;
    if (q <= 1.0)
        res = (k / h) * (-3 * q + 2.25 * q * q);
    else if (q < 2.0) {
        auto two_m_q = 2 - q;
        res = -0.75 * (k / h) * two_m_q * two_m_q;
    }
    return res;
}

static inline Vector3d gradW(Vector3d r, double h) {
    if (r.dot(r) == 0.0)
        return Vector3d(0.0);
    return dW(r.norm(), h) * r / r.norm();
}


SPH::SPH(int n)
{
	for(int i = 0; i < n; i++){
		shared_ptr<particle> new_particle(new particle);
		m_particle_list.push_back(new_particle);
        Vector3d zeros(0,0,0);
        new_particle->position = zeros;
        new_particle->velocity = zeros;
        new_particle->pressure = 0;
        new_particle->density = rho0;
        new_particle->mass = particle_radius * particle_radius * particle_radius * rho0;

	}
}

// void SPH::compute_all_Densities()
// {
//     for(unsigned int i = 0; i < m_particle_list.size(); i++)
//     {
//         shared_ptr<particle> cur = m_particle_list.at(i);
//         cur->density += single_density(cur->position);
//     }
// }


// void SPH::compute_all_Pressure()
// {
//     for(unsigned int i = 0; i < m_particle_list.size(); i++)
//     {
//          shared_ptr<particle> cur = m_particle_list.at(i);
//          cur->pressure = single_pressure(cur);
//     }
// }

Vector3d SPH::single_Viscousity(shared_ptr<particle> cur)
{
    Vector3d accerlation(0,0,0);
    auto ra = cur->position;
    auto va = cur->velocity;
    auto pa = cur->pressure;
    auto rho_a = cur->density;
    for (unsigned int pid = 0; pid < m_particle_list.size(); pid++)
    {
        shared_ptr<particle> neigh = m_particle_list.at(pid);
        if((neigh->position - ra).norm() < neighbor_radius)
        {
            auto rb = neigh->position;
            auto vb = neigh->velocity;
            auto rho_b = neigh->density;
            auto rab = ra - rb;
            accerlation += _mu / rho_a * neigh->mass * (vb-va) / rho_b * delta_square_W(rab, neighbor_radius);
        }
    }
    return accerlation;
}

Vector3d SPH::single_Momentum(shared_ptr<particle> cur)
{
    Vector3d accerlation(0,0,0);
    auto ra = cur->position;
    auto pa = cur->pressure;
    auto rho_a = cur->density;
    for (unsigned int pid = 0; pid < m_particle_list.size(); pid++)
    {
        shared_ptr<particle> neigh = m_particle_list.at(pid);
        if((neigh->position - ra).norm() < neighbor_radius)
        {
            auto rb = neigh->position;
            auto pb = neigh->pressure;
            auto rho_b = neigh->density;
            auto rab = ra - rb;
            accerlation += neigh->mass * (pa / pow(rho_a, 2) + pb / pow(rho_b, 2)) * delta_W(rab, neighbor_radius);
        }
    }
    return accerlation;

}

Vector3d SPH::single_dvdt(shared_ptr<particle> cur)
{
    Vector3d gravity(0,-0.1,0);
    Vector3d tot_acc(0,0,0);
    tot_acc = gravity - single_Momentum(cur) + single_Viscousity(cur);
    return tot_acc;
}

double SPH::single_pressure(shared_ptr<particle> cur)
{
    auto k = rho0 * 81.0f / 7.0f;
    return k * (std::pow(cur->density / rho0, 7) - 1.0);
}


double SPH::single_density(shared_ptr<particle> cur)
{
    auto ra = cur->position;
    auto va = cur->velocity;
    double density = 0;
    for (unsigned int pid = 0; pid < m_particle_list.size(); pid++)
    {
        shared_ptr<particle> neigh = m_particle_list.at(pid);
        if((neigh->position - ra).norm() < neighbor_radius)
        {

            auto rb = neigh->position;
            auto vb = neigh->velocity;
            auto vab = va-vb;
            auto rab = ra- rb;
            density += W(rab, neighbor_radius) * neigh->mass;
        }
    }
    return density;

}

void SPH::euler_step()
{
    //update 
    for (unsigned int pid = 0; pid < m_particle_list.size(); pid++)
    {
        shared_ptr<particle> cur = m_particle_list.at(pid);
        cur->density = single_density(cur);
    }
    for (unsigned int pid = 0; pid < m_particle_list.size(); pid++)
    {
        shared_ptr<particle> cur = m_particle_list.at(pid);
        cur->pressure = single_pressure(cur);
    }
    for (unsigned int pid = 0; pid < m_particle_list.size(); pid++)
    {
        shared_ptr<particle> cur = m_particle_list.at(pid);
        auto tot_acc = single_dvdt(cur);
        cur->velocity += tot_acc * dt;
    }
    for (unsigned int pid = 0; pid < m_particle_list.size(); pid++)
    {
        shared_ptr<particle> cur = m_particle_list.at(pid);
        cur->position += cur->velocity * dt; 
    }


}


