#include "positioninit.h"
#include "math.h"

using namespace std;
using namespace Eigen;

PositionInit::PositionInit(float interParticleDist)
{
    m_dx = interParticleDist;
}

Eigen::Vector3f& PositionInit::getPt(int i ){
    return m_pts[i];
}

int PositionInit::getNumParticles(){
    return m_pts.size();
}

void PositionInit::addSphere(Eigen::Vector3f center, float radius){
    int boxSz = (int) ceil(radius / m_dx);

    for (int i = - boxSz; i < boxSz; ++i){
        for (int j = - boxSz; j < boxSz; ++j){
            for (int k = - boxSz; k < boxSz; ++k){
                if ((i*i + j*j + k*k)*m_dx < radius * radius){
                    m_pts.emplace_back(i * m_dx + center(0), j * m_dx + center(1), k * m_dx + center(2));
                }
            }
        }
    }
}

void PositionInit::addBox(Eigen::Vector3f center, float width, float height, float thick){
    int widthInt = (int) ceil(width / m_dx / 2.f);
    int heightInt = (int) ceil(height / m_dx / 2.f);
    int thickInt = (int) ceil(thick / m_dx / 2.f);
    for (int i = - widthInt; i < widthInt; ++i){
        for (int j = - heightInt; j < heightInt; ++j){
            for (int k = - thickInt; k < thickInt; ++k){
                m_pts.emplace_back(i * m_dx + center(0), j * m_dx + center(1), k * m_dx + center(2));
            }
        }
    }
}

