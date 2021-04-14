#ifndef MAGNETICINIT_H
#define MAGNETICINIT_H

#include "Eigen/Dense"
#include <vector>

typedef enum {POINT_SOURCE, CONST_SOURCE} MagneticType;

typedef struct MagneticSource {
    MagneticType type;
    // for point source
    Eigen::Vector3d position;
    Eigen::Vector3d dipole;
    //for constant field
   Eigen::Vector3d B;

    MagneticSource(MagneticType type, Eigen::Vector3d position, Eigen::Vector3d dipole, Eigen::Vector3d B):
        type(type),
        position(position),
        dipole(dipole),
        B(B)

    {
    //inside parameterized constructor
    }

} MagneticSource_t;


class MagneticInit
{
public:
    MagneticInit();

    Eigen::Vector3d getMagneticField(const Eigen::Vector3d &position) const;

    void addPointSource(Eigen::Vector3d position, Eigen::Vector3d m);

    void addConstField(Eigen::Vector3d B);
private:
    std::vector <MagneticSource_t> m_src;

    Eigen::Vector3d getMagneticField(const MagneticSource_t& src, const Eigen::Vector3d &position) const;

    Eigen::Vector3d getDipoleField(const Eigen::Vector3d &srcPos,
                                   const Eigen::Vector3d &targetPos, const Eigen::Vector3d dipole) const;
};

#endif // MAGNETICINIT_H
