#ifndef SIMDUMP_H
#define SIMDUMP_H
#include <string>

namespace SimDump {
    void dump_full_sim(std::string outfile, int nparticles, float radius, int nframes, float timestep, int nevery);
}

#endif // SIMDUMP_H
