#include <iostream>
#include <vector>
#include <cassert>
#include <cmath>
#include <gmsh.h>

#include "CalcMesh.h"

int main() {
    const unsigned int GMSH_TETR_CODE = 4;

    gmsh::initialize();
    gmsh::model::add("invader");

    // Считаем STL
    try {
        gmsh::merge("invader.stl"); 
    } catch(...) {
        gmsh::logger::write("Could not load STL mesh: bye!");
        gmsh::finalize();
        return -1;
    }

    // Восстановим геометрию
    double angle = 40;
    bool forceParametrizablePatches = false;
    bool includeBoundary = true;
    double curveAngle = 180;
    gmsh::model::mesh::classifySurfaces(angle * M_PI / 180., includeBoundary, forceParametrizablePatches, curveAngle * M_PI / 180.);
    gmsh::model::mesh::createGeometry();

    // Зададим объём по считанной поверхности
    std::vector<std::pair<int, int> > s;
    gmsh::model::getEntities(s, 2);
    std::vector<int> sl;
    for(auto surf : s) sl.push_back(surf.second);
    int l = gmsh::model::geo::addSurfaceLoop(sl);
    gmsh::model::geo::addVolume({l});

    gmsh::model::geo::synchronize();

    // Зададим мелкость желаемой сетки
    int f = gmsh::model::mesh::field::add("MathEval");
    gmsh::model::mesh::field::setString(f, "F", "3");
    gmsh::model::mesh::field::setAsBackgroundMesh(f);

    // Построим сетку
    gmsh::model::mesh::generate(3);

    // Теперь извлечём из gmsh данные об узлах сетки
    std::vector<double> nodesCoord;
    std::vector<std::size_t> nodeTags;
    std::vector<double> parametricCoord;
    gmsh::model::mesh::getNodes(nodeTags, nodesCoord, parametricCoord);

    // И данные об элементах сетки тоже извлечём, нам среди них нужны только тетраэдры, которыми залит объём
    std::vector<std::size_t>* tetrsNodesTags = nullptr;
    std::vector<int> elementTypes;
    std::vector<std::vector<std::size_t>> elementTags;
    std::vector<std::vector<std::size_t>> elementNodeTags;
    gmsh::model::mesh::getElements(elementTypes, elementTags, elementNodeTags);
    for(unsigned int i = 0; i < elementTypes.size(); i++) {
        if(elementTypes[i] != GMSH_TETR_CODE)
            continue;
        tetrsNodesTags = &elementNodeTags[i];
    }

    if(tetrsNodesTags == nullptr) {
        std::cout << "Can not find tetra data. Exiting." << std::endl;
        gmsh::finalize();
        return -2;
    }

    std::cout << "The model has " <<  nodeTags.size() << " nodes and " << tetrsNodesTags->size() / 4 << " tetrs." << std::endl;

    // На всякий случай проверим, что номера узлов идут подряд и без пробелов
    for(unsigned int i = 0; i < nodeTags.size(); ++i) {
        // Индексация в gmsh начинается с 1, а не с нуля. Ну штош, значит так.
        assert(i == nodeTags[i] - 1);
    }
    // И ещё проверим, что в тетраэдрах что-то похожее на правду лежит.
    assert(tetrsNodesTags->size() % 4 == 0);

    // TODO: неплохо бы полноценно данные сетки проверять, да

    CalcMesh mesh(nodesCoord, *tetrsNodesTags);

    std::string moving = "combined"; // тут можно ихменить тип движения space invader'а
    
    mesh.setVelocities(

        (moving == "ears" || moving == "combined") ? [](double x, double y, double z, double t) { 
            if (y >= ((t < 0.24) ? 7.1 + 5 * t : 8.3 - 5 * (t - 0.24))) {
                if (t < 0.12) {
                    return 20.0;
                } else if (t >= 0.12 && t < 0.36) {
                    return -20.0;
                } else if (t >= 0.36 && t < 0.48) {
                    return 20.0;
                } else {
                    return 0.0;
                }
            } else {
                return 0.0;
            }
        } : [](double x, double y, double z, double t) { return 0.0; },

        (moving == "hands" || moving == "combined") ? [](double x, double y, double z, double t) {
            if ((x < -25.96 || x > 25.96) && (y < ((t < 0.24) ? -7.1 + 5 * t : -5.9 - 5 * (t - 0.24)))) {
                if (t < 0.12) {
                    return 20.0 + 5;
                } else if (t >= 0.12 && t < 0.24) {
                    return -20.0 + 5;
                } else if (t >= 0.24 && t < 0.36) {
                    return -20.0 - 5;
                } else if (t >= 0.36 && t < 0.48) {
                    return 20.0 - 5;
                } else {
                    return 5.0;
                }
            } else {
                if (t < 0.24) {
                    return 5.0;
                } else {
                    return -5.0;
                }
            } 
        } : [](double x, double y, double z, double t) { return 0.0; },

        [](double x, double y, double z, double t) { return 0.0; }
    );

    gmsh::finalize();

    mesh.snapshot(0);
    double tau = 0.01;
    for(unsigned int step = 1; step < 48; step++) {
        mesh.doTimeStep(tau);
        mesh.snapshot(step);
    }

    return 0;
}
