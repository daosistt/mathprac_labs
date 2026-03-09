#ifndef CALC_MESH_H
#define CALC_MESH_H

#include <vector>
#include <string>
#include "CalcNode.h"
#include "Element.h"
#include <functional>

class CalcMesh {
  protected:
    std::vector<CalcNode> nodes;
    std::vector<Element> elements;
    double currentTime;

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
};

#endif // CALC_MESH_H
