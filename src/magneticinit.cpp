#include "magneticinit.h"
#include <math.h>

using namespace Eigen;
using namespace std;

MagneticInit::MagneticInit()
{

}

Eigen::Vector3d MagneticInit::getMagneticField(const Eigen::Vector3d &position) const{
    Vector3d B = Vector3d::Zero();
    for (const auto &src: m_src){
        B += getMagneticField(src, position);
    }
    return B;
}

void MagneticInit::addPointSource(Eigen::Vector3d position, Eigen::Vector3d m){
     m_src.emplace_back(POINT_SOURCE, position, m, Vector3d::Zero());
}

void MagneticInit::addConstField(Eigen::Vector3d B){
    m_src.emplace_back(CONST_SOURCE, Vector3d::Zero(), Vector3d::Zero(), B);
}


Eigen::Vector3d MagneticInit::getMagneticField(const MagneticSource_t& src, const Eigen::Vector3d &position) const{
    switch (src.type){
        case POINT_SOURCE:
            return getDipoleField(src.position, position, src.dipole);
        case CONST_SOURCE:
            return src.B;
    }
}

Eigen::Vector3d MagneticInit::getDipoleField(const Eigen::Vector3d &srcPos, const Eigen::Vector3d &targetPos,
                                             const Eigen::Vector3d dipole) const{
    Vector3d r = targetPos - srcPos;
    double R = r.norm();
    return 1 / 4.0 / M_PI * (3 * r * dipole.dot(r) / pow(R, 5.0)  - dipole / pow(R, 3.0));

}

Eigen::Matrix3d MagneticInit::getMagneticFieldGrad (const Eigen::Vector3d &position) const{
    Matrix3d grad;
    for (const auto & src: m_src){
        if (src.type == POINT_SOURCE){
            double rm = position.dot(src.dipole);
            double r = position.norm();
            Matrix3d mrT = src.dipole * position.transpose();
            Matrix3d rrT  = position * position.transpose();
            grad += 3.0/ 4.0 * M_PI * (rm / pow(r, 5.0) * Matrix3d::Identity()
                                       + (mrT + mrT.transpose())/pow(r, 5.0)
                                       - 5 * rrT * rm / pow(r, 7.0));
        }
    }

    return grad;

}


