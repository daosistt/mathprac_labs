#ifndef CALC_MESH_H
#define CALC_MESH_H

#include <vector>
#include <string>
#include "CalcNode.h"
#include "Element.h"
#include <functional>

struct CenterOfMass {
    double centerX, centerY, centerZ;
    double vx, vy, vz;
    
    double yaw, pitch, roll;
   
    std::function<double(double)> yaw_rate_func;
    std::function<double(double)> pitch_rate_func;
    std::function<double(double)> roll_rate_func;

    CenterOfMass() 
        : centerX(0.0), centerY(0.0), centerZ(0.0)
        , vx(0.0), vy(0.0), vz(0.0)
        , yaw(0.0), pitch(0.0), roll(0.0) {}
};

class CalcMesh {
  protected:
    std::vector<CalcNode> nodes;
    std::vector<Element> elements;
    double currentTime;

    CenterOfMass com;

    std::function<double(double,double,double,double)> vxFunc;
    std::function<double(double,double,double,double)> vyFunc;
    std::function<double(double,double,double,double)> vzFunc;

  public:
    CalcMesh(const std::vector<double>& nodesCoords, const std::vector<std::size_t>& tetrsPoints);

    void doTimeStep(double tau);
    void snapshot(unsigned int snap_number);
    void setVelocities(std::function<double(double,double,double,double)> fx,
                       std::function<double(double,double,double,double)> fy,
                       std::function<double(double,double,double,double)> fz);

    void updateCenterOfMass();

    void setAngularVelocity(
        std::function<double(double)> yaw_rate,
        std::function<double(double)> pitch_rate,
        std::function<double(double)> roll_rate
    );
};

#endif // CALC_MESH_H
