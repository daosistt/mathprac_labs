#ifndef CALC_NODE_H
#define CALC_NODE_H

class CalcNode {
    friend class CalcMesh;

  protected:
    double x;
    double y;
    double z;

    double smth;
    
    double vx;
    double vy;
    double vz;

  public:
    CalcNode();
    CalcNode(double x, double y, double z, double smth, double vx, double vy, double vz);

    void move(double tau);
};

#endif // CALC_NODE_H
