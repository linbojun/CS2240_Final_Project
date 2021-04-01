#include <SPH.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <graphics/Shader.h>
#include <iostream>
#include <Eigen/Dense>

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
        return Vector3d(0, 0, 0);
    return dW(r.norm(), h) * r / r.norm();
}

inline bool inBounds(const Vector3d& pos) {
    return 1 >= pos[0] && pos[0] >= 0 && 1 >= pos[1] && pos[1] >= 0 && 1 >= pos[2] && pos[2] >= 0;
}

inline Vector3i SPH::gridPlace(const Vector3d& pos) {
    return Vector3i(pos[0] * (_grid_segs - 1), pos[1] * (_grid_segs - 1), pos[2] * (_grid_segs - 1));
}

inline int SPH::gridPlaceIndex(const Vector3i& place) {
    return place[0] * _grid_segs * _grid_segs + place[1] * _grid_segs + place[2];
}

inline int SPH::gridIndex(Vector3d& pos) {
    return gridPlaceIndex(gridPlace(pos));
}


Shape get_sphere_shape(float r, int res) {
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
        for(int phis = 2; phis < res-1; phis++) {
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

SPH::SPH(int n, float radius) : m_radius(radius)
{
    _neighbor_radius = m_radius * 1.3f;
    _grid_segs = 1 / _neighbor_radius;
    _voxel_len = 1.f/_grid_segs;
    _max_grid_search = ceil(_neighbor_radius / _voxel_len);
    _dh =_neighbor_radius;

    m_grid.resize(_grid_segs * _grid_segs * _grid_segs);
	for(int i = 0; i < n; i++){
		shared_ptr<particle> new_particle(new particle);
		m_particle_list.push_back(new_particle);
        Vector3d zeros(0,0,0);
        updateParticlePos(i, Vector3d::Random() * 0.5 + Vector3d(0.5, 0.5, 0.5));
        new_particle->velocity = zeros;
        new_particle->pressure = 0;
        new_particle->density = _rho0;
        new_particle->mass = m_radius * m_radius * m_radius * _rho0;
        //m_p_shapes.push_back(get_sphere_shape(adius));
	}
}

vector<shared_ptr<particle>> SPH::find_neighs(int pi)
{
    shared_ptr<particle> &cur = m_particle_list.at(pi);
    vector<shared_ptr<particle>> neighs;
    assert(inBounds(cur->position));
    Vector3i place = gridPlace(cur->position);
    for(int x = max(0, place[0] - _max_grid_search); x <= min(_grid_segs - 1, place[0] + _max_grid_search); x++) {
        for(int y = max(0, place[1] - _max_grid_search); y <= min(_grid_segs - 1, place[1] + _max_grid_search); y++) {
            for(int z = max(0, place[2] - _max_grid_search); z <= min(_grid_segs - 1, place[2] + _max_grid_search); z++) {
                int idx = gridPlaceIndex(Vector3i(x, y, z));
                for(int poss_neigh : m_grid[idx]) {
                    if(poss_neigh == pi)
                        continue;
                    shared_ptr<particle>& ptcl_other = m_particle_list.at(poss_neigh);
                    auto dist = ptcl_other->position - cur->position;
                    if(dist.norm() < _neighbor_radius) {
                        neighs.push_back(ptcl_other);
                    }
                }
            }
        }
    }

    return neighs;

}


void SPH::update(float seconds)
{
    for(unsigned int i = 0; i < m_particle_list.size(); i++)
    {
        shared_ptr<particle> cur = m_particle_list.at(i);
        cur->neighs  = find_neighs(i);
        auto dvdt = total_dvdt(cur);
        cur->drhodt = single_drhodt(cur);
        cur->dvdt = dvdt;
    }
    for(unsigned int i = 0; i < m_particle_list.size(); i++)
    {
        shared_ptr<particle> cur = m_particle_list.at(i);
        updateParticlePos(i, cur->position + 0.5 * seconds * cur->velocity);
        cur->velocity += 0.5 * seconds * cur->dvdt;
        cur->density += seconds * cur->drhodt;
        cur->pressure = single_pressure(cur);
    }
    for(unsigned int i = 0; i < m_particle_list.size(); i++)
    {
        shared_ptr<particle> cur = m_particle_list.at(i);
        updateParticlePos(i, cur->position + seconds * cur->velocity);
        cur->velocity += seconds * cur->dvdt;
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
    Vector3d dvdt(0, 0, 0);
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
    for(int i = 0; i < m_particle_list.size(); i++)
    {
        auto &ptcl = m_particle_list[i];
        auto pos = ptcl->position;
        auto velo = ptcl->velocity;
        double k = 0.4;
        for(int i = 0; i < 3; i++)
        {
            if (pos[i] < 0)
            {
                pos(i,0) = 0;
                if(velo(i,0) < 0)
                    velo(i,0) = -k * velo(i,0);
            }
            if (pos(i,0) > 1)
            {
               pos(i,0) = 1;
                if (velo(i,0) > 0)
                    velo(i,0) = -k *velo(i,0);
            }
        }
        updateParticlePos(i, pos);
        ptcl->velocity = velo;
    }
}


void SPH::updateParticlePos(int i, Vector3d newPos, bool initializing) {
    auto particle = m_particle_list[i];
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

void SPH::draw(Shader *shader) {
    static Shape shape = get_sphere_shape(m_radius, 4);
    for(int i = 0; i < m_particle_list.size(); i++) {
        auto& ptcl = m_particle_list[i];
        //cout << "pos=" << ptcl->position;
        shape.setModelMatrix(Eigen::Affine3f(Eigen::Translation3f(ptcl->position[0], ptcl->position[1], ptcl->position[2])));
        shape.draw(shader);
    }
}
