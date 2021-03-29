#include "simulation.h"

#include <iostream>

#include "graphics/MeshLoader.h"

#include "SPH.h"
#include <QApplication>
#include <memory>

using namespace Eigen;

Simulation::Simulation()
{
}

void Simulation::init()
{
    m_sph = std::make_shared<SPH>(1000);
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

void Simulation::update(float seconds)
{
    m_sph->update(seconds);
}

void Simulation::draw(Shader *shader)
{
    std::cout << "Rendered frame!\n";
    m_sph->draw(shader);
    m_ground.draw(shader);
}
