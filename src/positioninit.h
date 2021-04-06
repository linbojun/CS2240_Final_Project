#ifndef POSITIONINIT_H
#define POSITIONINIT_H

#include "Eigen/Dense"
#include <vector>

class PositionInit
{
public:
    PositionInit(float interParticleDist);

    void addSphere(Eigen::Vector3f center, float radius);

    void addBox(Eigen::Vector3f center, float width, float height, float thick);

    Eigen::Vector3f& getPt(int i);

    int getNumParticles();



private:
    float m_dx;

    std::vector<Eigen::Vector3f> m_pts;

};

#endif // POSITIONINIT_H
