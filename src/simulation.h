#ifndef SIMULATION_H
#define SIMULATION_H

#include "graphics/shape.h"
#include "magneticSPH.h"
#include <memory>
#include "Kernel.h"
#include "wcsph.h"
class Shader;

class Simulation
{
public:
    Simulation();

    void init(int num_particles, float timestep, float radius);

    void update(float seconds);

    void draw(Shader *shader);


private:
    std::shared_ptr<WCSPH> m_sph;
    Shape m_ground;
    float m_timestep;
};

#endif // SIMULATION_H
