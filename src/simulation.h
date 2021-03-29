#ifndef SIMULATION_H
#define SIMULATION_H

#include "graphics/shape.h"
#include "SPH.h"
#include <memory>

class Shader;

class Simulation
{
public:
    Simulation();

    void init();

    void update(float seconds);

    void draw(Shader *shader);


private:
    std::shared_ptr<SPH> m_sph;
    Shape m_ground;
};

#endif // SIMULATION_H
