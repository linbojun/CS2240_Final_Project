#include <SPH.h>
#include <math.h>

using namespace Eigen;
using namespace std;

static inline double W(double r, double h) {
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
        new_particle->density = _rho0;
        new_particle->mass = _particle_radius * _particle_radius * _particle_radius * _rho0;

	}
}

vector<shared_ptr<particle>> SPH::find_neighs(shared_ptr<particle> cur)
{
    vector<shared_ptr<particle>> neighs;
    for(auto ptcl: m_particle_list)
    {
        if(ptcl.get() == cur.get())
            continue;
        auto dist = ptcl->position - cur->position;
        if(dist.norm() < _neighbor_radius)
        {
            neighs.push_back(ptcl);
        }
    }
    return neighs;

}


void SPH::euler_step()
{
    for(unsigned int i = 0; i < m_particle_list.size(); i++)
    {
        shared_ptr<particle> cur = m_particle_list.at(i);
        cur->neighs  = find_neighs(cur);
        auto dvdt = total_dvdt(cur);
        cur->drhodt = single_drhodt(cur);
        cur->dvdt = dvdt;
    }
    for(unsigned int i = 0; i < m_particle_list.size(); i++)
    {
        shared_ptr<particle> cur = m_particle_list.at(i);
        cur->position += 0.5 * _dt * cur->velocity;
        cur->velocity += 0.5 * _dt * cur->dvdt;
        cur->density += _dt * cur->drhodt;
        cur->pressure = single_pressure(cur);
    }
    for(unsigned int i = 0; i < m_particle_list.size(); i++)
    {
        shared_ptr<particle> cur = m_particle_list.at(i);
        cur->position += _dt * cur->velocity;
        cur->velocity += _dt * cur->dvdt;
    }
    boundry_collision();
}

Vector3d SPH::total_dvdt(shared_ptr<particle> cur)
{
    Vector3d dvdt(0,0,0);
    dvdt += momentum_dvdt(cur);
    dvdt += viscosity_dvdt(cur);
    return dvdt;
}

Vector3d SPH::momentum_dvdt(shared_ptr<particle> cur)
{
    Vector3d gravity(0,-0.1,0);
    Vector3d dvdt(0,0,0);
    auto ra = cur->position;
    auto pa = cur->pressure;
    auto rho_a  = cur->density;
    for(auto neigh: cur->neighs)
    {
        auto rb = neigh->position;
        auto pb = neigh->pressure;
        auto rho_b = neigh->density;
        dvdt += -neigh->mass * (pa / (rho_a * rho_a) + pb / (rho_b * rho_b)) * gradW(ra - rb, _neighbor_radius);
    }
    dvdt += gravity;
    return dvdt;
}

Vector3d SPH::viscosity_dvdt(shared_ptr<particle> cur)
{
    constexpr float eps = 0.01f;
    Vector3d dvdt(0.0);
    auto ra = cur->position;
    auto va = cur->velocity;
    auto rho_a = cur->density;
    for(auto neigh: cur->neighs)
    {
        auto rho_b = neigh->density;
        auto rb = neigh->position;
        auto vb = neigh->velocity;
        auto vab = va - vb;
        auto rab = ra - rb;
        if(vab.dot(rab) < 0)
        {
            auto v = -2 * _dh * _c / (rho_a +rho_b);
            auto pi_ab = -v * vab.dot(rab) / (rab.dot(rab) + eps * _dh * _dh);
             dvdt += neigh->mass * pi_ab * gradW(rab, _dh);
        }
    }
    return dvdt;
}

double SPH::single_drhodt(shared_ptr<particle> cur)
{
    double drhodt = 0;
    auto ra = cur->position;
    auto va = cur->velocity;
    auto neighs = cur->neighs;
    for(auto neigh: neighs)
    {
        auto rb = neigh->position;
        auto vb = neigh->velocity;
        auto vab = va - vb;
        auto rab = ra - rb;
        drhodt += neigh->mass * vab.dot(gradW(rab, _dh));
       }
       return drhodt;
}

double SPH::single_pressure(shared_ptr<particle> cur)
{
    auto p0 = _rho0 * _c * _c / _gamma;
    return p0 * (pow(cur->density / _rho0, _gamma) - 1.0);
}

void SPH::boundry_collision()
{
    for(auto ptcl: m_particle_list)
    {
        auto pos = ptcl->position;
        auto velo = ptcl->velocity;
        double k = 0.4;
        for(int i = 0; i < 3; i++)
        {
            if (pos(i,0) < 0)
            {
                pos(i,0) = 0;
                if(velo(i,0) < 0)
                    velo(i,0) += (k+1) * velo(i,0);
            }
            if (pos(i,0) > 1)
            {
               pos(i,0) = 1;
                if (velo(i,0) > 0)
                    velo(i,0) += -(1 + k) *velo(i,0);
            }
        }
    }
}
