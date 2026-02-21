#include <set>
#include <gmsh.h>

int main(int argc, char **argv)
{
    gmsh::initialize();

    gmsh::model::add("circle");

    double lc = 1e-2;

    // Центр
    gmsh::model::geo::addPoint(0, 0, 0, lc, 1);

    // Полуокружность
    gmsh::model::geo::addPoint(0, .1, 0, lc, 2);
    gmsh::model::geo::addPoint(0, -.1, 0, lc, 3);

    // Дуги
    gmsh::model::geo::addCircleArc(2, 1, 3, 1);
    gmsh::model::geo::addCircleArc(3, 1, 2, 2);

    // Контур и сетка
    gmsh::model::geo::addCurveLoop({1, 2}, 1);
    gmsh::model::geo::addPlaneSurface({1}, 1);

    // Синхронизация, генерация сетки и отображение
    gmsh::model::geo::synchronize();

    gmsh::model::mesh::generate(2);

    std::set<std::string> args(argv, argv + argc);
    if(!args.count("-nopopup")) gmsh::fltk::run();

    gmsh::finalize();

    return 0;
}

