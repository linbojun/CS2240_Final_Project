#include "simulation.h"

#include <iostream>

#include "graphics/MeshLoader.h"

#include "magneticwcsph.h"
#include <QApplication>
#include <memory>

using namespace Eigen;

const double H_VALUE = 1.0;

Simulation::Simulation()
{
}

void Simulation::init(int num_particles, float timestep, float radius)
{

//    m_sph = std::make_shared<MagneticSPH>(num_particles, radius, H_VALUE);
//    m_sph = std::make_shared<SPH>(num_particles, radius);
    m_sph = std::make_shared<MagneticWCSPH>(0, radius, H_VALUE);
    std::cout << "Num particles: " << m_sph->getNumParticle() << "\n";
    m_timestep = timestep;
    std::vector<Vector3f> groundVerts;
    std::vector<Vector3i> groundFaces;
    groundVerts.emplace_back(0, -0.1, 0);
    groundVerts.emplace_back(0, -0.1, 1);
    groundVerts.emplace_back(1, -0.1, 1);
    groundVerts.emplace_back(1, -0.1, 0);
    groundFaces.emplace_back(0, 1, 2);
    groundFaces.emplace_back(0, 2, 3);
    m_ground.init(groundVerts, groundFaces, true);
}

void Simulation::update(float time)
{
    for(int i = 0; i < 1; i++)
        m_sph->update(m_timestep);
}

void Simulation::draw(Shader *shader)
{
    //std::cout << "Rendered frame!\n";
    m_sph->draw(shader);
    m_ground.draw(shader);
}
