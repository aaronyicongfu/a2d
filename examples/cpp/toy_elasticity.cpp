#include <complex>
#include <cstdint>
#include <iomanip>
#include <iostream>

#include "a2dtmp3d.h"
#include "fem/elasticity.h"
#include "fem/helmholtz.h"
#include "fem/model.h"
#include "utils/a2dprofiler.h"
#include "utils/a2dvtk.h"

using namespace A2D;

void main_body(int argc, char* argv[]) {
  static const index_t SPATIAL_DIM = 3;
  typedef index_t I;
  // typedef std::complex<double> T;
  typedef double T;
  typedef BasisOps<SPATIAL_DIM, HexTriLinearBasisFunc, Hex8ptQuadrature> Basis;
  typedef ElasticityPDEInfo<SPATIAL_DIM, I, T> PDE;

  const index_t nx = 64;
  const index_t ny = 64;
  const index_t nz = 64;
  const index_t nnodes = (nx + 1) * (ny + 1) * (nz + 1);
  const index_t nelems = nx * ny * nz;
  const index_t nbcs = (ny + 1) * (nz + 1);

  auto model = std::make_shared<FEModel<I, T, PDE>>(nnodes, nbcs);
  auto element = std::make_shared<NonlinElasticityElement<I, T, Basis>>(nelems);
  model->add_element(element);

  // Set the boundary conditions
  auto bcs = model->get_bcs();
  index_t index = 0;
  for (int k = 0; k < nz + 1; k++) {
    for (int j = 0; j < ny + 1; j++) {
      int i = 0;
      int node = i + (nx + 1) * (j + (ny + 1) * k);

      // Set the boundary conditions
      bcs(index, 0) = node;
      for (int ii = 0; ii < 3; ii++) {
        bcs(index, 1) |= 1U << ii;
      }
      index++;
    }
  }

  // Set the connectivity
  auto conn = element->get_conn();
  for (int k = 0; k < nz; k++) {
    for (int j = 0; j < ny; j++) {
      for (int i = 0; i < nx; i++) {
        int elem = i + nx * (j + ny * k);

        int conn_coord[8];
        for (int kk = 0, index = 0; kk < 2; kk++) {
          for (int jj = 0; jj < 2; jj++) {
            for (int ii = 0; ii < 2; ii++, index++) {
              conn_coord[index] =
                  (i + ii) + (nx + 1) * ((j + jj) + (ny + 1) * (k + kk));
            }
          }
        }

        // Convert to the correct connectivity
        conn(elem, 0) = conn_coord[0];
        conn(elem, 1) = conn_coord[1];
        conn(elem, 2) = conn_coord[3];
        conn(elem, 3) = conn_coord[2];

        conn(elem, 4) = conn_coord[4];
        conn(elem, 5) = conn_coord[5];
        conn(elem, 6) = conn_coord[7];
        conn(elem, 7) = conn_coord[6];
      }
    }
  }

  // Set the node locations
  auto X = model->get_nodes();
  for (int k = 0; k < nz + 1; k++) {
    for (int j = 0; j < ny + 1; j++) {
      for (int i = 0; i < nx + 1; i++) {
        int node = i + (nx + 1) * (j + (ny + 1) * k);

        X(node, 0) = 1.0 * i / nx;
        X(node, 1) = 1.0 * j / ny;
        X(node, 2) = 1.0 * k / nz;
      }
    }
  }

  // Set the node locations - Note: This must be done after setting the
  // connectivity!
  model->init();

  // Set the element
  T q = 5.0, E = 70e3, nu = 0.3;
  T density = 1.0, design_stress = 1e3;
  auto constitutive = std::make_shared<TopoIsoConstitutive<I, T, Basis>>(
      element, q, E, nu, density, design_stress);
  model->add_constitutive(constitutive);

  // Create the design vector
  A2D::CLayout<1> design_layout(model->nnodes);
  auto x = std::make_shared<A2D::MultiArray<T, A2D::CLayout<1>>>(design_layout);

  // Set the design variable values
  x->fill(1.0);
  model->set_design_vars(x);

  // Set up the stress functional
  auto functional = std::make_shared<Functional<I, T, PDE>>();
  auto agg_functional =
      std::make_shared<TopoVonMisesAggregation<I, T, Basis>>(constitutive);
  functional->add_functional(agg_functional);

  // Compute the Jacobian matrix
  Timer* t;
  auto J = model->new_matrix();
  model->jacobian(J);

  int num_levels = 3;
  double omega = 0.6667;
  double epsilon = 0.01;
  bool print_info = true;
  auto amg = model->new_amg(num_levels, omega, epsilon, J, print_info);

  // Set the residuals and apply the boundary conditions
  auto solution = model->new_solution();
  auto residual = model->new_solution();
  solution->fill(1.0);
  model->zero_bcs(solution);

  residual->zero();
  // BSRMatVecMult(*J, *solution, *residual);
  for (int k = nz / 4; k < 3 * nz / 4; k++) {
    int node = nx + (nx + 1) * (0 + (ny + 1) * k);
    (*residual)(node, 1) = -1e2;
  }
  model->zero_bcs(residual);

  // Compute the solution
  index_t monitor = 10;
  index_t max_iters = 80;
  solution->zero();
  amg->cg(*residual, *solution, monitor, max_iters);

  ToVTK<decltype(element->get_conn()), decltype(model->get_nodes())> vtk(conn,
                                                                         X);
  vtk.write_mesh();
  vtk.write_sol("ux", *solution, 0);
  vtk.write_sol("uy", *solution, 1);
  vtk.write_sol("uz", *solution, 2);
}

int main(int argc, char* argv[]) {
  Timer t("main");
  main_body(argc, argv);
}
