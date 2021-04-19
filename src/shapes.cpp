#include "shapes.h"
#include "graphics/shape.h"
#include "graphics/Shader.h"
#define _USE_MATH_DEFINES
#include <math.h>
using namespace Eigen;

Shape get_sphere_shape(float r, int res) {
    std::vector<Vector3f> points;
    std::vector<Vector3i> faces;
    for(int thetas = 0; thetas < res; thetas++) {
        for(int phis = 1; phis < res-1; phis++) {
            float x = r * sin(phis*M_PI/res)*cos(thetas*2*M_PI/res);
            float y = r * sin(phis*M_PI/res)*sin(thetas*2*M_PI/res);
            float z = r * cos(phis*M_PI/res);
            points.push_back(Vector3f(x, y, z));
        }
    }
    // idx = thetas*18 + phis - 1
    for(int thetas = 1; thetas < res+1; thetas++) {
        for(int phis = 2; phis < res; phis++) {
            faces.push_back(Vector3i((thetas%res)*(res-2) + phis - 1, (thetas-1)*(res-2) + phis - 1, (thetas-1)*(res-2) + phis - 2));
            faces.push_back(Vector3i((thetas%res)*(res-2) + phis - 1, (thetas-1)*(res-2) + phis - 2, (thetas%res)*(res-2) + phis - 2));
        }
    }
    int topi = points.size();
    points.push_back(Vector3f(0, 0, r));
    int boti = points.size();
    points.push_back(Vector3f(0, 0, -r));
    for(int thetas = 0; thetas < res; thetas++) {
        faces.push_back(Vector3i(topi, thetas*(res-2), ((thetas+1)%res)*(res-2)));
        faces.push_back(Vector3i(boti, thetas*(res-2) + (res-3), ((thetas+1)%res)*(res-2) + (res-3)));
    }
    Shape shape;
    shape.init(points, faces, true);
    shape.setVertices(points);
    return shape;
}

Shape getLineShape() {
    // center of coordinate sys = one end of line, line pointing in +x direction
    // triangular prism
    std::vector<Eigen::Vector3f> points = {
        Eigen::Vector3f(0, cos(0)*0.05, sin(0)*0.05),
        Eigen::Vector3f(0, cos(M_PI/3)*0.05, sin(M_PI/3)*0.05),
        Eigen::Vector3f(0, cos(2*M_PI/3)*0.05, sin(2*M_PI/3)*0.05),
        Eigen::Vector3f(1, cos(0)*0.05, sin(0)*0.05),
        Eigen::Vector3f(1, cos(M_PI/3)*0.05, sin(M_PI/3)*0.05),
        Eigen::Vector3f(1, cos(2*M_PI/3)*0.05, sin(2*M_PI/3)*0.05),
    };
    std::vector<Eigen::Vector3i> faces(8);
    faces[0] =
        Eigen::Vector3i(0, 1, 2);
            faces[1] =
        Eigen::Vector3i(3, 4, 5);
            faces[2] =
        Eigen::Vector3i(0, 3, 1);
            faces[3] =
        Eigen::Vector3i(1, 3, 4);
            faces[4] =
        Eigen::Vector3i(1, 4, 2);
            faces[5] =
        Eigen::Vector3i(2, 4, 5);
            faces[6] =
        Eigen::Vector3i(2, 5, 0);
            faces[7] =
        Eigen::Vector3i(0, 5, 3);
    Shape shape;
    shape.init(points, faces, true);
    shape.setVertices(points);
    return shape;
}


