"""Convert circle.msh (Gmsh MSH4) → circle.xdmf / circle_boundary.xdmf.

DOLFINx требует:
  - отдельные файлы для 2D-ячеек (треугольники) и 1D-границы (линии)
  - только один тип ячеек в файле (не Mixed)
  - точки только XY (2D)

Usage:
    python convert_mesh.py

Requires:
    pip install meshio h5py numpy
"""

import numpy as np
import meshio

# --- читаем MSH -------------------------------------------------------
msh = meshio.read("circle.msh")

points_2d = msh.points[:, :2]   # убираем Z=0

# --- треугольники (домен) ---------------------------------------------
cells_tri = msh.cells_dict.get("triangle")
if cells_tri is None:
    raise RuntimeError("No triangle cells found in circle.msh")

try:
    tri_tags = msh.cell_data_dict["gmsh:physical"]["triangle"]
except KeyError:
    tri_tags = np.ones(len(cells_tri), dtype=np.int32)

mesh_tri = meshio.Mesh(
    points=points_2d,
    cells=[("triangle", cells_tri)],
    cell_data={"markers": [tri_tags]},
)
meshio.write("circle.xdmf", mesh_tri)
print("Written: circle.xdmf + circle.h5")
