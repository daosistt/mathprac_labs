#include <iostream>
#include <vector>
#include <cassert>
#include <cmath>
#include <gmsh.h>

#include "CalcMesh.h"

int main() {
    const unsigned int GMSH_TETR_CODE = 4;

    gmsh::initialize();
    gmsh::model::add("jet_fighter");

    // Считаем STL
    try {
        gmsh::merge("jet_fighter.stl"); 
    } catch(...) {
        gmsh::logger::write("Could not load STL mesh: bye!");
        gmsh::finalize();
        return -1;
    }

    // Восстановим геометрию
    gmsh::option::setNumber("Mesh.StlRemoveDuplicateTriangles", 1);
    double angle = 40;
    bool forceParametrizablePatches = true;
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

  // some Claude advices ;)
    gmsh::option::setNumber("Mesh.Algorithm3D", 4);
    gmsh::option::setNumber("Mesh.Algorithm", 6);
    gmsh::option::setNumber("Mesh.OptimizeNetgen", 1);
    gmsh::option::setNumber("Mesh.Optimize", 1);
    gmsh::option::setNumber("Mesh.AnisoMax", 1e33);
    gmsh::option::setNumber("Mesh.CharacteristicLengthFromCurvature", 1);
    gmsh::option::setNumber("Mesh.MinimumCirclePoints", 20);

    // Зададим мелкость желаемой сетки
    int f = gmsh::model::mesh::field::add("MathEval");
    gmsh::model::mesh::field::setString(f, "F", "16");
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

    // TODO чуть поменять CalcMesh и doTimeStep учётом центра масс
    // чтобы добавить pitch, roll, yaw  или петли

    mesh.setVelocities(
        [](double x, double y, double z, double t) { return 0.0; },

        [](double x, double y, double z, double t) { return -2000.0; },

        [](double x, double y, double z, double t) { return 0.0; }
    );

    gmsh::finalize();

    mesh.snapshot(0);
    double tau = 0.01;
    for(unsigned int step = 1; step < 150; step++) {
        mesh.doTimeStep(tau);
        mesh.snapshot(step);
    }

    return 0;
}
