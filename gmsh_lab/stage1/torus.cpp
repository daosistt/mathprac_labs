#include <set>
#include <gmsh.h>

int main(int argc, char **argv)
{
    gmsh::initialize();

    gmsh::model::add("cylinder");

    double lc = 0.1;
    double R = 0.5;
    double r_1 = 0.2;
    double r_2 = 0.15;

    gmsh::model::geo::addPoint(0, R, 0, lc, 1);

    // Окружность (внешняя) в поперечном сечении тора (YZ)
    gmsh::model::geo::addPoint(0, R + r_1, 0, lc, 2);
    gmsh::model::geo::addPoint(0, R, r_1, lc, 3);
    gmsh::model::geo::addPoint(0, R - r_1, 0, lc, 4);
    gmsh::model::geo::addPoint(0, R, -r_1, lc, 5);

    gmsh::model::geo::addCircleArc(2, 1, 3, 1);
    gmsh::model::geo::addCircleArc(3, 1, 4, 2);
    gmsh::model::geo::addCircleArc(4, 1, 5, 3);
    gmsh::model::geo::addCircleArc(5, 1, 2, 4);

    // Окружность (внутренняя) в поперечном сечении тора (YZ)
    gmsh::model::geo::addPoint(0, R + r_2, 0, lc, 6);
    gmsh::model::geo::addPoint(0, R, r_2, lc, 7);
    gmsh::model::geo::addPoint(0, R - r_2, 0, lc, 8);
    gmsh::model::geo::addPoint(0, R, -r_2, lc, 9);

    gmsh::model::geo::addCircleArc(6, 1, 7, 5);
    gmsh::model::geo::addCircleArc(7, 1, 8, 6);
    gmsh::model::geo::addCircleArc(8, 1, 9, 7);
    gmsh::model::geo::addCircleArc(9, 1, 6, 8);

    // Поверхность (кольцо)
    gmsh::model::geo::addCurveLoop({1, 2, 3, 4, -5, -6, -7, -8}, 1);
    gmsh::model::geo::addPlaneSurface({1}, 1);

    gmsh::model::geo::synchronize();


    double angle = 2.0 * M_PI / 3.0; // 2pi / 3
    int current_surf = 1;

    // Вращение кольца 3 раза на 120 градусов
    for (int i = 0; i < 3; i++) {
        std::vector<std::pair<int,int>> dim_tags;

        gmsh::model::geo::revolve({{2, current_surf}},
            .0, .0, .0,
            .0, .0, 1.0,
            angle, 
            dim_tags
        );

        gmsh::model::geo::synchronize();

        // Ищем торец (dim == 2)
        for (auto& t : dim_tags) {
            if (t.first == 2) {
                current_surf = t.second;
                break;
            }
        }
    }

    // Синхронизация, генерация сетки и отображение
    gmsh::model::geo::synchronize();

    gmsh::model::mesh::generate(3);

    std::set<std::string> args(argv, argv + argc);
    if(!args.count("-nopopup")) gmsh::fltk::run();

    gmsh::finalize();

    return 0;
}
