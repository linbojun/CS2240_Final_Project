#include "simulation.h"

#include <iostream>

#include "graphics/MeshLoader.h"

#include "SPH.h"
#include <QApplication>
#include <memory>

using namespace Eigen;

const double H_VALUE = 1.0;

Simulation::Simulation()
{
}

void Simulation::init(int num_particles, float timestep, float radius)
{
    m_sph = std::make_shared<SPH>(num_particles, radius);//, H_VALUE);
    m_timestep = timestep;
    std::vector<Vector3f> groundVerts;
    std::vector<Vector3i> groundFaces;
    groundVerts.emplace_back(-5, -1, -5);
    groundVerts.emplace_back(-5, -1, 5);
    groundVerts.emplace_back(5, -1, 5);
    groundVerts.emplace_back(5, -1, -5);
    groundFaces.emplace_back(0, 1, 2);
    groundFaces.emplace_back(0, 2, 3);
    m_ground.init(groundVerts, groundFaces, true);
}

void Simulation::update(float time)
{
    m_sph->update(m_timestep);
}

void Simulation::draw(Shader *shader)
{
    std::cout << "Rendered frame!\n";
    m_sph->draw(shader);
    m_ground.draw(shader);
}
