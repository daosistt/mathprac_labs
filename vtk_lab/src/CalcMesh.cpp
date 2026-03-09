#include "CalcMesh.h"
#include <vtkDoubleArray.h>
#include <vtkPoints.h>
#include <vtkPointData.h>
#include <vtkTetra.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkUnstructuredGrid.h>
#include <vtkSmartPointer.h>
#include <cmath>
#include <string>

CalcMesh::CalcMesh(const std::vector<double>& nodesCoords, const std::vector<std::size_t>& tetrsPoints) : currentTime(0.0)  {
    // Пройдём по узлам в модели gmsh
    nodes.resize(nodesCoords.size() / 3);
    for(unsigned int i = 0; i < nodesCoords.size() / 3; i++) {
        // Координаты заберём из gmsh
        double pointX = nodesCoords[i*3];
        double pointY = nodesCoords[i*3 + 1];
        double pointZ = nodesCoords[i*3 + 2];
        // Модельная скалярная величина распределена как-то вот так
        double smth = pow(pointX, 2) + pow(pointY, 2) + pow(pointZ, 2);
        nodes[i] = CalcNode(pointX, pointY, pointZ, smth, 0.0, 0.0, 0.0);
    }

    // Пройдём по элементам в модели gmsh
    elements.resize(tetrsPoints.size() / 4);
    for(unsigned int i = 0; i < tetrsPoints.size() / 4; i++) {
        elements[i].nodesIds[0] = tetrsPoints[i*4] - 1;
        elements[i].nodesIds[1] = tetrsPoints[i*4 + 1] - 1;
        elements[i].nodesIds[2] = tetrsPoints[i*4 + 2] - 1;
        elements[i].nodesIds[3] = tetrsPoints[i*4 + 3] - 1;
    }
}

void CalcMesh::doTimeStep(double tau) {
    // Обновляем скорости по текущему времени
    for(auto& node : nodes) {
        node.vx = vxFunc(node.x, node.y, node.z, currentTime);
        node.vy = vyFunc(node.x, node.y, node.z, currentTime);
        node.vz = vzFunc(node.x, node.y, node.z, currentTime);
    }

    // Двигаем точки
    for(unsigned int i = 0; i < nodes.size(); i++) {
        nodes[i].move(tau);
    }

    currentTime += tau;
}

void CalcMesh::snapshot(unsigned int snap_number) {
    // Сетка в терминах VTK
    vtkSmartPointer<vtkUnstructuredGrid> unstructuredGrid = vtkSmartPointer<vtkUnstructuredGrid>::New();
    // Точки сетки в терминах VTK
    vtkSmartPointer<vtkPoints> dumpPoints = vtkSmartPointer<vtkPoints>::New();

    // Скалярное поле на точках сетки
    auto smth = vtkSmartPointer<vtkDoubleArray>::New();
    smth->SetName("smth");

    // Векторное поле на точках сетки
    auto vel = vtkSmartPointer<vtkDoubleArray>::New();
    vel->SetName("velocity");
    vel->SetNumberOfComponents(3);

    // Обходим все точки нашей расчётной сетки
    for(unsigned int i = 0; i < nodes.size(); i++) {
        // Вставляем новую точку в сетку VTK-снапшота
        dumpPoints->InsertNextPoint(nodes[i].x, nodes[i].y, nodes[i].z);

        // Добавляем значение векторного поля в этой точке
        double _vel[3] = {nodes[i].vx, nodes[i].vy, nodes[i].vz};
        vel->InsertNextTuple(_vel);

        // И значение скалярного поля тоже
        smth->InsertNextValue(nodes[i].smth);
    }

    // Грузим точки в сетку
    unstructuredGrid->SetPoints(dumpPoints);

    // Присоединяем векторное и скалярное поля к точкам
    unstructuredGrid->GetPointData()->AddArray(vel);
    unstructuredGrid->GetPointData()->AddArray(smth);

    // А теперь пишем, как наши точки объединены в тетраэдры
    for(unsigned int i = 0; i < elements.size(); i++) {
        auto tetra = vtkSmartPointer<vtkTetra>::New();
        tetra->GetPointIds()->SetId(0, elements[i].nodesIds[0]);
        tetra->GetPointIds()->SetId(1, elements[i].nodesIds[1]);
        tetra->GetPointIds()->SetId(2, elements[i].nodesIds[2]);
        tetra->GetPointIds()->SetId(3, elements[i].nodesIds[3]);
        unstructuredGrid->InsertNextCell(tetra->GetCellType(), tetra->GetPointIds());
    }

    // Создаём снапшот в файле с заданным именем
    std::string fileName = "jet-step-" + std::to_string(snap_number) + ".vtu";
    vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
    writer->SetFileName(fileName.c_str());
    writer->SetInputData(unstructuredGrid);
    writer->Write();
}

void CalcMesh::setVelocities(std::function<double(double,double,double,double)> fx,
                             std::function<double(double,double,double,double)> fy,
                             std::function<double(double,double,double,double)> fz) {
    // Сохраняем вызываемые объекты определяющие векторное поле скоростей
    for(auto& node : nodes) {
        vxFunc = fx;
        vyFunc = fy;
        vzFunc = fz;
    }
}
