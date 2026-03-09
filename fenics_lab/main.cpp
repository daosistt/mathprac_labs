#include "poisson.h"
#include <basix/finite-element.h>
#include <cmath>
#include <dolfinx.h>
#include <dolfinx/fem/Constant.h>
#include <dolfinx/fem/petsc.h>
#include <dolfinx/la/petsc.h>
#include <petscmat.h>
#include <petscsys.h>
#include <petscsystypes.h>
#include <utility>
#include <vector>

// Загрузка сетки из XDMF
#include <dolfinx/io/XDMFFile.h>

using namespace dolfinx;
using T = PetscScalar;
using U = typename dolfinx::scalar_value_t<T>;

int main(int argc, char* argv[])
{
  dolfinx::init_logging(argc, argv);
  PetscInitialize(&argc, &argv, nullptr, nullptr);
  {
  // Тут xdmf, так как собирал dolfinx без gmshio, а пересобирать пока не оч хочу  
  auto mesh = std::make_shared<mesh::Mesh<U>>(
        [&]()
        {
          auto coord_elem = fem::CoordinateElement<U>(
              mesh::CellType::triangle, 1);
          io::XDMFFile xfile(MPI_COMM_WORLD, "../mesh/circle.xdmf", "r");
          return xfile.read_mesh(coord_elem, mesh::GhostMode::shared_facet,
                                 "Grid", "/Xdmf/Domain");
        }());

    mesh->topology_mutable()->create_connectivity(1, 2); // грани — нужны для BC

    auto basix_elem = basix::create_element<U>(
        basix::element::family::P, basix::cell::type::triangle, 1,
        basix::element::lagrange_variant::unset,
        basix::element::dpc_variant::unset, false);

    auto V = std::make_shared<fem::FunctionSpace<U>>(
        fem::create_functionspace<U>(
            mesh, std::make_shared<fem::FiniteElement<U>>(basix_elem)));

    auto kappa = std::make_shared<fem::Constant<T>>(1.0);

    auto f = std::make_shared<fem::Function<T, U>>(V);
    f->interpolate(
        [](auto x) -> std::pair<std::vector<T>, std::vector<std::size_t>>
        {
          constexpr double A     = 10.0;
          constexpr double sigma = 0.4;   // ширина гауссиана
          std::vector<T> vals;
          vals.reserve(x.extent(1));
          for (std::size_t p = 0; p < x.extent(1); ++p)
          {
            double r2 = x(0, p) * x(0, p) + x(1, p) * x(1, p);
            vals.push_back(A * std::exp(-r2 / (sigma * sigma)));
          }
          return {vals, {vals.size()}};
        });

    auto g = std::make_shared<fem::Function<T, U>>(V);

    fem::Form<T, U> a = fem::create_form<T, U>(
        *form_poisson_a,
        std::vector<std::shared_ptr<const fem::FunctionSpace<U>>>{V, V},
        std::map<std::string, std::shared_ptr<const fem::Function<T, U>>>{},
        std::map<std::string, std::shared_ptr<const fem::Constant<T>>>{{"kappa", kappa}},
        std::map<fem::IntegralType,
                 std::vector<std::pair<int, std::span<const int>>>>{},
        std::vector<std::reference_wrapper<const mesh::EntityMap>>{},
        mesh);

    fem::Form<T, U> L = fem::create_form<T, U>(
        *form_poisson_L,
        std::vector<std::shared_ptr<const fem::FunctionSpace<U>>>{V},
        std::map<std::string, std::shared_ptr<const fem::Function<T, U>>>{
            {"f", f}, {"g", g}},
        std::map<std::string, std::shared_ptr<const fem::Constant<T>>>{},
        std::map<fem::IntegralType,
                 std::vector<std::pair<int, std::span<const int>>>>{},
        std::vector<std::reference_wrapper<const mesh::EntityMap>>{},
        mesh);

    constexpr double R   = 1.0;
    constexpr double eps = 1.0e-6;

    std::vector facets = mesh::locate_entities_boundary(
        *mesh, 1,
        [&](auto x)
        {
          std::vector<std::int8_t> marker(x.extent(1), false);
          for (std::size_t p = 0; p < x.extent(1); ++p)
          {
            double r2 = x(0, p) * x(0, p) + x(1, p) * x(1, p);
            if (std::abs(std::sqrt(r2) - R) < eps)
              marker[p] = true;
          }
          return marker;
        });

    std::vector bdofs = fem::locate_dofs_topological(
        *V->mesh()->topology_mutable(), *V->dofmap(), 1, facets);

    fem::DirichletBC<T, U> bc(T(0), bdofs, V);

    auto u = std::make_shared<fem::Function<T, U>>(V);

    la::petsc::Matrix A(fem::petsc::create_matrix(a), false);
    la::Vector<T> b(L.function_spaces()[0]->dofmap()->index_map,
                    L.function_spaces()[0]->dofmap()->index_map_bs());

    MatZeroEntries(A.mat());
    fem::assemble_matrix(la::petsc::Matrix::set_block_fn(A.mat(), ADD_VALUES),
                         a, {bc});
    MatAssemblyBegin(A.mat(), MAT_FLUSH_ASSEMBLY);
    MatAssemblyEnd(A.mat(), MAT_FLUSH_ASSEMBLY);
    fem::set_diagonal<T>(la::petsc::Matrix::set_fn(A.mat(), INSERT_VALUES), *V,
                         {bc});
    MatAssemblyBegin(A.mat(), MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(A.mat(), MAT_FINAL_ASSEMBLY);

    std::ranges::fill(b.array(), 0);
    fem::assemble_vector(b.array(), L);
    fem::apply_lifting(b.array(), {a}, {{bc}}, {}, T(1));
    b.scatter_rev(std::plus<T>());
    bc.set(b.array(), std::nullopt);

    // ------------------------------------------------------------------
    // 7. Решение: прямой решатель LU
    // ------------------------------------------------------------------
    la::petsc::KrylovSolver lu(MPI_COMM_WORLD);
    la::petsc::options::set("ksp_type", "preonly");
    la::petsc::options::set("pc_type", "lu");
    lu.set_from_options();

    lu.set_operator(A.mat());
    la::petsc::Vector _u(la::petsc::create_vector_wrap(*u->x()), false);
    la::petsc::Vector _b(la::petsc::create_vector_wrap(b), false);
    lu.solve(_u.vec(), _b.vec());

    u->x()->scatter_fwd();

    io::VTKFile file(MPI_COMM_WORLD, "u_circle.pvd", "w");
    file.write<T, U>(
        std::vector<std::reference_wrapper<const fem::Function<T, U>>>{*u},
        0.0);
  }

  PetscFinalize();
  return 0;
}
