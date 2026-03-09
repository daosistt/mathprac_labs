#include <set>
#include <string>
#include <gmsh.h>

int main(int argc, char **argv)
{
    gmsh::initialize();

    gmsh::model::add("circle");

    double lc = 1e-1;
    double R  = 1.0;   // радиус

    // Центр
    gmsh::model::geo::addPoint(0.0,  0.0, 0.0, lc, 1);

    // Точки на окружности (для дуг)
    gmsh::model::geo::addPoint( R,   0.0, 0.0, lc, 2);
    gmsh::model::geo::addPoint( 0.0,  R,  0.0, lc, 3);
    gmsh::model::geo::addPoint(-R,   0.0, 0.0, lc, 4);
    gmsh::model::geo::addPoint( 0.0, -R,  0.0, lc, 5);

    // Четыре дуги, образующие полную окружность
    gmsh::model::geo::addCircleArc(2, 1, 3, 1);
    gmsh::model::geo::addCircleArc(3, 1, 4, 2);
    gmsh::model::geo::addCircleArc(4, 1, 5, 3);
    gmsh::model::geo::addCircleArc(5, 1, 2, 4);

    // Контур и поверхность
    gmsh::model::geo::addCurveLoop({1, 2, 3, 4}, 1);
    gmsh::model::geo::addPlaneSurface({1}, 1);

    gmsh::model::geo::synchronize();

    // Граница (для BC Дирихле)
    gmsh::model::addPhysicalGroup(1, {1, 2, 3, 4}, 1);
    gmsh::model::setPhysicalName(1, 1, "boundary");

    // Объём (поверхность)
    gmsh::model::addPhysicalGroup(2, {1}, 2);
    gmsh::model::setPhysicalName(2, 2, "domain");

    gmsh::model::mesh::generate(2);

    gmsh::write("circle.msh");

    gmsh::finalize();
    return 0;
}
