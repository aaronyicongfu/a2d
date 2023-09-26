#ifndef A2D_FE_BASE_H
#define A2D_FE_BASE_H

#include <algorithm>
#include <list>
#include <set>

#include "a2ddefs.h"
#include "multiphysics/feelement.h"
#include "multiphysics/femesh.h"

namespace A2D {

/**
 * @brief Element base class.
 *
 * This defines an element that is compatible with the given vector and matrix
 * types.
 */
template <typename T, class VecType, class MatType>
class ElementBase {
 public:
  virtual ~ElementBase() {}

  virtual const ElementMeshBase& get_mesh(FEVarType wrt) const = 0;

  // Add the residual
  virtual void add_residual(T alpha, VecType& data, VecType& geo, VecType& sol,
                            VecType& res) = 0;

  // Add the Jacobian
  virtual void add_jacobian(T alpha, VecType& data, VecType& geo, VecType& sol,
                            MatType& mat) = 0;

  // Add the adjoint-residual product
  virtual void add_adjoint_res_product(FEVarType wrt, T alpha, VecType& data,
                                       VecType& geo, VecType& sol,
                                       VecType& adjoint, VecType& dfdx) = 0;
};

/**
 * @brief Assemble the contributions from a number of elements
 *
 * @tparam VecType The type of vector
 * @tparam MatType The type of matrix to use
 */
template <typename T, class VecType, class MatType>
class ElementAssembler {
 public:
  typedef std::shared_ptr<ElementBase<T, VecType, MatType>> ElemPtr;

  ElementAssembler() {}
  virtual ~ElementAssembler() {}

  // Add an element at the specified index
  void add_element(ElemPtr element) { elements.push_back(element); }

  // Get the CSR structure for all the elements in the mesh
  void get_bsr_data(const index_t block_size, index_t& nrows,
                    std::vector<index_t>& rowp, std::vector<index_t>& cols) {
    std::set<std::pair<index_t, index_t>> pairs;
    typename std::list<ElemPtr>::iterator it;
    for (it = elements.begin(); it != elements.end(); ++it) {
      const ElementMeshBase& mesh = (*it)->get_mesh(FEVarType::STATE);
      mesh.add_matrix_pairs(block_size, pairs);
    }
    ElementMeshBase::create_block_csr(pairs, nrows, rowp, cols);
  }

  // Add the residual from all the elements
  void add_residual(T alpha, VecType& data, VecType& geo, VecType& sol,
                    VecType& res) {
    typename std::list<ElemPtr>::iterator it;
    for (it = elements.begin(); it != elements.end(); ++it) {
      (*it)->add_residual(alpha, data, geo, sol, res);
    }
  }

  // Add the Jacobian contributions from all the elements
  void add_jacobian(T alpha, VecType& data, VecType& geo, VecType& sol,
                    MatType& mat) {
    typename std::list<ElemPtr>::iterator it;
    for (it = elements.begin(); it != elements.end(); ++it) {
      (*it)->add_jacobian(alpha, data, geo, sol, mat);
    }
  }

  // Add the adjoint-residual product
  virtual void add_adjoint_res_product(FEVarType wrt, T alpha, VecType& data,
                                       VecType& geo, VecType& sol,
                                       VecType& adjoint, VecType& dfdx) {
    typename std::list<ElemPtr>::iterator it;
    for (it = elements.begin(); it != elements.end(); ++it) {
      (*it)->add_adjoint_res_product(wrt, alpha, data, geo, sol, adjoint, dfdx);
    }
  }

 private:
  std::list<ElemPtr> elements;
};

/*
  Integrand-type finite-element
*/
template <typename T, class Integrand, class Quadrature, class DataBasis,
          class GeoBasis, class Basis, class VecType, class MatType>
class ElementIntegrand : public ElementBase<T, VecType, MatType> {
 public:
  // Set the element vector to use
  template <typename T0, class BasisType, class VectorType>
  using ElementVector = ElementVector_Serial<T0, BasisType, VectorType>;

  // Set the element matrix to use
  template <typename T0, class BasisType, class VectorType>
  using ElementMatrix = ElementMat_Serial<T0, BasisType, VectorType>;

  virtual const Integrand& get_integrand() = 0;

  // Set the meshes
  void set_meshes(std::shared_ptr<ElementMesh<DataBasis>> data,
                  std::shared_ptr<ElementMesh<GeoBasis>> geo,
                  std::shared_ptr<ElementMesh<Basis>> sol) {
    data_mesh = data;
    geo_mesh = geo;
    sol_mesh = sol;
  }

  const ElementMeshBase& get_mesh(FEVarType wrt) const {
    if (wrt == FEVarType::DATA) {
      return *data_mesh;
    } else if (wrt == FEVarType::GEOMETRY) {
      return *geo_mesh;
    } else {
      return *sol_mesh;
    }
  }

  // Add the residual to the residual vector
  void add_residual(T alpha, VecType& data, VecType& geo, VecType& sol,
                    VecType& res) {
    const Integrand& integrand = this->get_integrand();
    ElementVector<T, DataBasis, VecType> elem_data(*data_mesh, data);
    ElementVector<T, GeoBasis, VecType> elem_geo(*geo_mesh, geo);
    ElementVector<T, Basis, VecType> elem_sol(*sol_mesh, sol);
    ElementVector<T, Basis, VecType> elem_res(*sol_mesh, res);
    element.template add_residual<FEVarType::STATE>(
        integrand, alpha, elem_data, elem_geo, elem_sol, elem_res);
  }

  // Add the Jacobian
  void add_jacobian(T alpha, VecType& data, VecType& geo, VecType& sol,
                    MatType& mat) {
    const Integrand& integrand = this->get_integrand();
    ElementVector<T, DataBasis, VecType> elem_data(*data_mesh, data);
    ElementVector<T, GeoBasis, VecType> elem_geo(*geo_mesh, geo);
    ElementVector<T, Basis, VecType> elem_sol(*sol_mesh, sol);
    ElementMatrix<T, Basis, MatType> elem_jac(*sol_mesh, mat);
    element.template add_jacobian<FEVarType::STATE, FEVarType::STATE>(
        integrand, alpha, elem_data, elem_geo, elem_sol, elem_jac);
  }

  // Add the adjoint-residual product
  void add_adjoint_res_product(FEVarType wrt, T alpha, VecType& data,
                               VecType& geo, VecType& sol, VecType& adjoint,
                               VecType& dfdx) {
    const Integrand& integrand = this->get_integrand();
    ElementVector<T, DataBasis, VecType> elem_data(*data_mesh, data);
    ElementVector<T, GeoBasis, VecType> elem_geo(*geo_mesh, geo);
    ElementVector<T, Basis, VecType> elem_sol(*sol_mesh, sol);
    ElementVector<T, Basis, VecType> elem_adjoint(*sol_mesh, adjoint);

    if (wrt == FEVarType::DATA) {
      ElementVector<T, DataBasis, VecType> elem_dfdx(*data_mesh, dfdx);
      element.template add_jacobian_product<FEVarType::DATA, FEVarType::STATE>(
          integrand, alpha, elem_data, elem_geo, elem_sol, elem_adjoint,
          elem_dfdx);

    } else if (wrt == FEVarType::GEOMETRY) {
      ElementVector<T, GeoBasis, VecType> elem_dfdx(*geo_mesh, dfdx);
      element
          .template add_jacobian_product<FEVarType::GEOMETRY, FEVarType::STATE>(
              integrand, alpha, elem_data, elem_geo, elem_sol, elem_adjoint,
              elem_dfdx);

    } else if (wrt == FEVarType::STATE) {
      ElementVector<T, Basis, VecType> elem_dfdx(*sol_mesh, dfdx);
      element.template add_jacobian_product<FEVarType::STATE, FEVarType::STATE>(
          integrand, alpha, elem_data, elem_geo, elem_sol, elem_adjoint,
          elem_dfdx);
    }
  }

 protected:
  FiniteElement<T, Integrand, Quadrature, DataBasis, GeoBasis, Basis> element;
  std::shared_ptr<ElementMesh<DataBasis>> data_mesh;
  std::shared_ptr<ElementMesh<GeoBasis>> geo_mesh;
  std::shared_ptr<ElementMesh<Basis>> sol_mesh;
};

/**
 * @brief Functional base class
 *
 * This class defines a functional evaluation that is compatible with the
 * given vector and matrix types
 */
template <typename T, class VecType>
class FunctionalBase {
 public:
  virtual ~FunctionalBase() {}

  // Add the residual
  virtual T evaluate(VecType& data, VecType& geo, VecType& sol) = 0;

  // Add to the derivative type
  virtual void add_derivative(FEVarType wrt, T alpha, VecType& data,
                              VecType& geo, VecType& sol, VecType& dfdx) = 0;
};

/*
  Integral type functional
*/
template <typename T, class Integrand, class Quadrature, class DataBasis,
          class GeoBasis, class Basis, class VecType>
class IntegralFunctional : public FunctionalBase<T, VecType> {
 public:
  // Set the element vector to use
  template <typename T0, class BasisType, class VectorType>
  using ElementVector = ElementVector_Serial<T0, BasisType, VectorType>;

  // Set the element matrix to use
  template <typename T0, class BasisType, class VectorType>
  using ElementMatrix = ElementMat_Serial<T0, BasisType, VectorType>;

  virtual const Integrand& get_integrand() = 0;

  // Set the meshes
  void set_meshes(std::shared_ptr<ElementMesh<DataBasis>> data,
                  std::shared_ptr<ElementMesh<GeoBasis>> geo,
                  std::shared_ptr<ElementMesh<Basis>> sol) {
    data_mesh = data;
    geo_mesh = geo;
    sol_mesh = sol;
  }

  // Add the residual to the residual vector
  T evaluate(VecType& data, VecType& geo, VecType& sol) {
    const Integrand& integrand = this->get_integrand();
    ElementVector<T, DataBasis, VecType> elem_data(*data_mesh, data);
    ElementVector<T, GeoBasis, VecType> elem_geo(*geo_mesh, geo);
    ElementVector<T, Basis, VecType> elem_sol(*sol_mesh, sol);
    T value = element.integrate(integrand, elem_data, elem_geo, elem_sol);
    return value;
  }

  // Add the adjoint-residual product
  void add_derivative(FEVarType wrt, T alpha, VecType& data, VecType& geo,
                      VecType& sol, VecType& dfdx) {
    const Integrand& integrand = this->get_integrand();
    ElementVector<T, DataBasis, VecType> elem_data(*data_mesh, data);
    ElementVector<T, GeoBasis, VecType> elem_geo(*geo_mesh, geo);
    ElementVector<T, Basis, VecType> elem_sol(*sol_mesh, sol);

    if (wrt == FEVarType::DATA) {
      ElementVector<T, DataBasis, VecType> elem_dfdx(*data_mesh, dfdx);
      element.template add_residual<FEVarType::DATA>(
          integrand, alpha, elem_data, elem_geo, elem_sol, elem_dfdx);

    } else if (wrt == FEVarType::GEOMETRY) {
      ElementVector<T, GeoBasis, VecType> elem_dfdx(*geo_mesh, dfdx);
      element.template add_residual<FEVarType::GEOMETRY>(
          integrand, alpha, elem_data, elem_geo, elem_sol, elem_dfdx);

    } else if (wrt == FEVarType::STATE) {
      ElementVector<T, Basis, VecType> elem_dfdx(*sol_mesh, dfdx);
      element.template add_residual<FEVarType::STATE>(
          integrand, alpha, elem_data, elem_geo, elem_sol, elem_dfdx);
    }
  }

 protected:
  FiniteElement<T, Integrand, Quadrature, DataBasis, GeoBasis, Basis> element;
  std::shared_ptr<ElementMesh<DataBasis>> data_mesh;
  std::shared_ptr<ElementMesh<GeoBasis>> geo_mesh;
  std::shared_ptr<ElementMesh<Basis>> sol_mesh;
};

}  // namespace A2D

#endif  //  A2D_FE_BASE_H