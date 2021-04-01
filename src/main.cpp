#include <QApplication>
#include "mainwindow.h"
#include <cstdlib>
#include <ctime>
#include <iostream>
#include "simdump.h"

void usage() {
    std::cout << "Usage: <simulation.exe> <draw/output> [num-particles=10000] [particle-rad=0.03] [timestep=0.01] [mode-dependent-params...]";
}

void usage_output() {
    std::cout << "Usage: <simulation.exe> output [num-particles=10000] [particle-rad=0.03] [timestep=0.01] [frameskip=10] [nframes=500] [out-file=\"out\\sim-out.txt\"]";
}

int main(int argc, char *argv[])
{
    if(argc < 2) {
        usage();
        exit(1);
    }
    char *mode = argv[1];
    int num_particles = 10000;
    float timestep = 0.01;
    float radius = 0.03;
    if(argc > 2) {
        num_particles = std::stoi(argv[2]);
        if(num_particles < 1) {
            usage();
            exit(1);
        }
    }
    if(argc > 3) {
        timestep = std::stof(argv[3]);
        if(timestep <= 0) {
            usage();
            exit(1);
        }
    }
    if(argc > 4) {
        radius = std::stof(argv[4]);
        if(radius <= 0) {
            usage();
            exit(1);
        }
    }
    if(!strcmp(mode, "draw")) {
        QApplication a(argc, argv);
        MainWindow w;
        srand (static_cast <unsigned> (time(0)));
        // We cannot use w.showFullscreen() here because on Linux that creates the
        // window behind all other windows, so we have to set it to fullscreen after
        // it has been shown.
        w.show();
        //w.setWindowState(w.windowState() | Qt::WindowFullScreen); // Comment out this line to have a windowed 800x600 game on startup.

        return a.exec();
    } else if(!strcmp(mode, "output")) {
        int frameskip = 10;
        std::string out_file = "out\\sim-out.txt";
        int nframes = 500;
        if(argc > 5) {
            frameskip = std::stoi(argv[5]);
            if(frameskip <= 0) {
                usage_output();
                exit(1);
            }
        }
        if (argc > 6) {
            nframes = std::stoi(argv[6]);
            if(nframes <= 0) {
                usage_output();
                exit(1);
            }
        }
        if (argc > 7) {
            out_file = argv[7];
        }
        SimDump::dump_full_sim(out_file, num_particles, radius, nframes, timestep, frameskip);
    } else {
        usage();
        exit(1);
    }
}

