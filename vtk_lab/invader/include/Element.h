#ifndef ELEMENT_H
#define ELEMENT_H

class Element {
    friend class CalcMesh;

  protected:
    unsigned long nodesIds[4];
};

#endif // ELEMENT_H
