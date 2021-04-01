#include "simdump.h"
#include <fstream>
#include <SPH.h>
#include <iostream>
#include <mainwindow.h>
#include "view.h"

#include "viewformat.h"

#include <QApplication>
#include <QKeyEvent>

#include <iostream>
namespace SimDump {
void dump_next_frame(std::ofstream &out, SPH &sim, float timestep, bool onlyUpdate) {
    if(!onlyUpdate) {
        auto &vec = sim.getParticles();
        for (auto& ptcl : vec) {
            out << ptcl->position[0] << "," << ptcl->position[1] << "," << ptcl->position[2] << "\n";
        }
        out << "---DONE---\n";
    }
    sim.update(timestep);
}

void dump_full_sim(std::string outfile, int nparticles, float radius, int nframes, float timestep, int nevery)
{
    std::ofstream out;
    out.open(outfile);
    SPH sim(nparticles, radius);
    out << nparticles << "\n";
    for(int i = 0; i < nframes; i++) {
        dump_next_frame(out, sim, timestep, i % nevery != 0);
        std::cout << "Frame " << i << " dumped.\n";
        fflush(stdout);
    }
    out.close();
}
}
