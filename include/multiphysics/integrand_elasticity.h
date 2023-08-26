#ifndef A2D_ELASTICITY_H
#define A2D_ELASTICITY_H

#include "a2denum.h"
#include "a2dmatops2d.h"
#include "a2dmatops3d.h"
#include "multiphysics/femapping.h"
#include "multiphysics/fespace.h"

namespace A2D {

template <typename T, index_t D>
class IntegrandTopoLinearElasticity {
 public:
  IntegrandTopoLinearElasticity(T E, T nu, T q) : q(q) {
    mu0 = 0.5 * E / (1.0 + nu);
    lambda0 = E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu));
  }

  // Number of dimensions
  static const index_t dim = D;

  // Number of data dimensions
  static const index_t data_dim = 1;

  // Space for the finite-element data
  using DataSpace = FESpace<T, data_dim, L2Space<T, data_dim, dim>>;

  // Space for the element geometry
  using FiniteElementGeometry = FESpace<T, dim, H1Space<T, dim, dim>>;

  // Finite element space
  using FiniteElementSpace = FESpace<T, dim, H1Space<T, dim, dim>>;

  // Mapping of the solution from the reference element to the physical element
  using SolutionMapping = InteriorMapping<T, dim>;

  // The type of matrix used to store data at each quadrature point
  static const index_t ncomp = FiniteElementSpace::ncomp;
  using QMatType = SymMat<T, ncomp>;

  // Data for the element
  T mu0;      // Second Lame parameter
  T lambda0;  // First Lame parameter
  T q;        // The RAMP penalty parameter

  /**
   * @brief Find the integral of the compliance over the entire domain
   *
   * @param wdetJ The determinant of the Jacobian times the quadrature weight
   * @param data The data at the quadrature point
   * @param geo The geometry at the quadrature point
   * @param s The solution at the quadurature point
   * @return T The integrand contribution
   */
  T integrand(T wdetJ, const DataSpace& data, const FiniteElementGeometry& geo,
              const FiniteElementSpace& s) const {
    T rho = data[0];
    T penalty = 1.0 / (1.0 + q * (1.0 - rho));

    // Get the constitutive data at the points
    T mu = penalty * mu0;
    T lambda = penalty * lambda0;

    // Extract the solution
    Mat<T, dim, dim> Ux = (s.template get<0>()).get_grad();

    // The Green-Langrange strain terms
    SymMat<T, dim> E;

    T output;
    MatLinearGreenStrain(Ux, E);
    SymmIsotropicEnergy(mu, lambda, E, output);

    return wdetJ * output;
  }

  /**
   * @brief Evaluate the weak form coefficients for linear elasticity
   *
   * @param wdetJ The quadrature weight times determinant of the Jacobian
   * @param data The data at the quadrature point
   * @param geo The geometry at the quadrature point
   * @param s The trial solution
   * @param coef Output weak form coefficients of the test space
   */
  KOKKOS_FUNCTION void weak(T wdetJ, const DataSpace& data,
                            const FiniteElementGeometry& geo,
                            const FiniteElementSpace& s,
                            FiniteElementSpace& coef) const {
    T rho = data[0];
    T penalty = 1.0 / (1.0 + q * (1.0 - rho));

    // Get the constitutive data at the points
    T mu = penalty * mu0;
    T lambda = penalty * lambda0;

    // Extract the trial solution gradient and the coefficient terms. Here
    // Uxb is the output computed as the derivative of the strain energy
    // w.r.t. Ux
    Mat<T, dim, dim> Ux0 = (s.template get<0>()).get_grad();
    Mat<T, dim, dim>& Uxb = (coef.template get<0>()).get_grad();
    ADMat<Mat<T, dim, dim>> Ux(Ux0, Uxb);

    // The Green-Lagrange strain terms
    SymMat<T, dim> E0, Eb;
    ADMat<SymMat<T, dim>> E(E0, Eb);

    // The strain energy output
    ADScalar<T> output;

    auto strain = MatLinearGreenStrain(Ux, E);
    auto energy = SymmIsotropicEnergy(mu, lambda, E, output);

    // Seed the output value with the wdetJ
    output.bvalue = wdetJ;

    // Reverse the derivatives through the code
    energy.reverse();
    strain.reverse();
  }

  // Evaluate the second order derivatives of the integral
  KOKKOS_FUNCTION void jacobian(T wdetJ, const DataSpace& data,
                                const FiniteElementGeometry& geo,
                                const FiniteElementSpace& s,
                                QMatType& jac) const {
    T rho = data[0];
    T penalty = 1.0 / (1.0 + q * (1.0 - rho));

    // Get the constitutive data at the points
    T mu = penalty * mu0;
    T lambda = penalty * lambda0;

    // Extract displacement gradient
    A2DMat<Mat<T, dim, dim>> Ux(s.template get<0>().get_grad());

    // The Green-Lagrange strain terms
    A2DMat<SymMat<T, dim>> E;

    // The strain energy output
    A2DScalar<T> output;

    // Evaluate the energy
    auto strain = MatLinearGreenStrain(Ux, E);
    auto energy = SymmIsotropicEnergy(mu, lambda, E, output);

    // Seed the output value with the wdetJ
    output.bvalue = wdetJ;

    // Reverse the derivatives through the code
    energy.reverse();
    strain.reverse();

    // Temporary vectors
    FiniteElementSpace p, Jp;

    for (index_t k = 0; k < FiniteElementSpace::ncomp; k++) {
      // Zero Hessian seeds because reverse operations are incremental
      Ux.Ah.zero();
      E.Ah.zero();
      output.hvalue = 0.0;

      // Select k-th column of the Hessian
      p.zero();
      p[k] = T(1.0);

      // Set projection direction
      Ux.set_pvalue((p.template get<0>()).get_grad());

      // Forward sweep to compute dependent directions
      strain.hforward();

      // Reverse sweep to compute Hessian-vector products
      energy.hreverse();
      strain.hreverse();

      // Set k-th column of the Hessian
      Ux.get_hvalue((Jp.template get<0>()).get_grad());
      for (index_t m = 0; m < ncomp; m++) {
        jac(m, k) = Jp[m];
      }
    }
  }

  /**
   * @brief Construct a JacVecProduct functor
   *
   * This functor computes a Jacobian-vector product of the
   *
   * @param wdetJ The quadrature weight times determinant of the Jacobian
   * @param data The data at the quadrature point
   * @param s The solution at the quadrature point
   */
  class JacVecProduct {
   public:
    KOKKOS_FUNCTION JacVecProduct(
        const IntegrandTopoLinearElasticity<T, D>& integrand, T wdetJ,
        const DataSpace& data, const FiniteElementGeometry& geo,
        const FiniteElementSpace& s)
        :  // Initialize constitutive data
          rho(data[0]),
          penalty(1.0 / (1.0 + integrand.q * (1.0 - rho))),
          mu(penalty * integrand.mu0),
          lambda(penalty * integrand.lambda0),

          // Initialize the displacement gradient
          Ux(s.template get<0>().get_grad()),

          // Compute the strain from the displacement gradient
          strain(Ux, E),

          // Compute the strain energy from the strain
          energy(mu, lambda, E, output) {
      // Set the seed on the derivative
      output.bvalue = wdetJ;

      // Reverse mode for the first derivative
      energy.reverse();
      strain.reverse();
    }

    KOKKOS_FUNCTION void operator()(const FiniteElementSpace& p,
                                    FiniteElementSpace& Jp) {
      Ux.set_pvalue((p.template get<0>()).get_grad());

      strain.hforward();
      energy.hreverse();
      strain.hreverse();

      Ux.get_hvalue((Jp.template get<0>()).get_grad());
    }

   private:
    T rho, penalty;
    T mu, lambda;
    A2DMat<Mat<T, dim, dim>> Ux;
    A2DMat<SymMat<T, dim>> E;
    A2DScalar<T> output;

    // Declare types of the operators
    decltype(MatLinearGreenStrain(Ux, E)) strain;
    decltype(SymmIsotropicEnergy(mu, lambda, E, output)) energy;
  };

  /**
   * @brief Construct a JacVecProduct functor
   *
   * @param wdetJ The quadrature weight times determinant of the Jacobian
   * @param data The data at the quadrature point
   * @param s The solution at the quadrature point
   */
  class AdjVecProduct {
   public:
    KOKKOS_FUNCTION AdjVecProduct(
        const IntegrandTopoLinearElasticity<T, D>& integrand, T wdetJ,
        const DataSpace& data, const FiniteElementGeometry& geo,
        const FiniteElementSpace& s)
        :  // Initialize constitutive data
          rho(data[0]),
          q(integrand.q),
          penalty(1.0 / (1.0 + q * (1.0 - rho))),
          mu0(integrand.mu0),
          lambda0(integrand.lambda0),
          mu(penalty * mu0),
          lambda(penalty * lambda0),

          // Initialize the displacement gradient
          Ux(s.template get<0>().get_grad()),

          // Compute the strain from the displacement gradient
          strain(Ux, E),

          // Compute the strain energy from the strain
          energy(mu, lambda, E, output) {
      // Set the seed on the derivative
      output.bvalue = wdetJ;

      // Reverse mode for the first derivative
      energy.reverse();
      strain.reverse();
    }

    KOKKOS_FUNCTION void operator()(const FiniteElementSpace& p,
                                    DataSpace& dfdx) {
      Ux.set_pvalue((p.template get<0>()).get_grad());

      strain.hforward();
      energy.hreverse();
      strain.hreverse();

      T denom = (1.0 + q * (1.0 - rho));
      dfdx[0] +=
          (mu0 * mu.hvalue + lambda0 * lambda.hvalue) * q / (denom * denom);
    }

   private:
    T rho;
    T q;
    T penalty;
    T mu0, lambda0;
    A2DScalar<T> mu, lambda;
    A2DMat<Mat<T, dim, dim>> Ux;
    A2DMat<SymMat<T, dim>> E;
    A2DScalar<T> output;

    // Declare types of the operators
    decltype(MatLinearGreenStrain(Ux, E)) strain;
    decltype(SymmIsotropicEnergy(mu, lambda, E, output)) energy;
  };
};

/*
  Evaluate the volume of the structure, given the constitutive class
*/
template <typename T, index_t C, index_t D, class Integrand>
class IntegrandTopoVolume {
 public:
  // Number of dimensions
  static const index_t dim = D;

  // Number of data dimensions
  static const index_t data_dim = 1;

  // Space for the finite-element data
  using DataSpace = typename Integrand::DataSpace;

  // Space for the element geometry
  using FiniteElementGeometry = typename Integrand::FiniteElementGeometry;

  // Finite element space
  using FiniteElementSpace = typename Integrand::FiniteElementSpace;

  // Mapping of the solution from the reference element to the physical element
  using SolutionMapping = typename Integrand::SolutionMapping;

  IntegrandTopoVolume() = default;

  /**
   * @brief Compute the integrand for this functional
   *
   * @param wdetJ The determinant of the Jacobian times the quadrature weight
   * @param data The data at the quadrature point
   * @param geo The geometry at the quadrature point
   * @param s The solution at the quadurature point
   * @return T The integrand contribution
   */
  T integrand(T wdetJ, const DataSpace& data, const FiniteElementGeometry& geo,
              const FiniteElementSpace& s) const {
    return wdetJ * data[0];
  }

  /**
   * @brief Derivative of the integrand with respect to the data
   *
   * @param wdetJ The determinant of the Jacobian times the quadrature weight
   * @param data The data at the quadrature point
   * @param geo The geometry at the quadrature point
   * @param s The solution at the quadurature point
   * @param dfdx The output derivative value
   */
  void data_derivative(T wdetJ, const DataSpace& data,
                       const FiniteElementGeometry& geo,
                       const FiniteElementSpace& s, DataSpace& dfdx) const {
    dfdx.zero();
    dfdx[0] = wdetJ;
  }
};

template <typename T, index_t D>
class IntegrandTopoBodyForce {
 public:
  // Number of dimensions
  static const index_t dim = D;

  // Number of data dimensions
  static const index_t data_dim = 1;

  // Space for the finite-element data
  using DataSpace = typename IntegrandTopoLinearElasticity<T, D>::DataSpace;

  // Space for the element geometry
  using FiniteElementGeometry =
      typename IntegrandTopoLinearElasticity<T, D>::FiniteElementGeometry;

  // Finite element space
  using FiniteElementSpace =
      typename IntegrandTopoLinearElasticity<T, D>::FiniteElementSpace;

  // Mapping of the solution from the reference element to the physical element
  using SolutionMapping = InteriorMapping<T, dim>;

  IntegrandTopoBodyForce(T q, const T tx_[]) : q(q) {
    for (index_t i = 0; i < dim; i++) {
      tx[i] = tx_[i];
    }
  }

  KOKKOS_FUNCTION void weak(T wdetJ, const DataSpace& data,
                            const FiniteElementGeometry& geo,
                            const FiniteElementSpace& s,
                            FiniteElementSpace& coef) const {
    T rho = data[0];
    T penalty = (q + 1.0) * rho / (q * rho + 1.0);

    // Add body force components
    Vec<T, dim>& Ub = (coef.template get<0>()).get_value();
    for (index_t i = 0; i < dim; i++) {
      Ub(i) = wdetJ * penalty * tx[i];
    }
  }

  class AdjVecProduct {
   public:
    KOKKOS_FUNCTION AdjVecProduct(
        const IntegrandTopoBodyForce<T, dim>& integrand, T wdetJ,
        const DataSpace& data, const FiniteElementGeometry& geo,
        const FiniteElementSpace& s)
        : q(integrand.q), rho(data[0]), wdetJ(wdetJ) {
      for (index_t i = 0; i < dim; i++) {
        tx[i] = integrand.tx[i];
      }
    }

    KOKKOS_FUNCTION void operator()(const FiniteElementSpace& psi,
                                    DataSpace& dfdx) {
      const Vec<T, dim>& Uadj = (psi.template get<0>()).get_value();
      T dpdrho = (q + 1.0) / ((q * rho + 1.0) * (q * rho + 1.0));

      for (index_t i = 0; i < dim; i++) {
        dfdx[0] += wdetJ * dpdrho * Uadj(i) * tx[i];
      }
    }

   private:
    T q, rho, wdetJ;
    T tx[dim];
  };

 private:
  T q;        // RAMP parameter
  T tx[dim];  // body force values
};

/*
  Evalute the KS functional of the stress, given the constitutive class
*/
template <typename T, index_t D>
class IntegrandTopoVonMisesKS {
 public:
  IntegrandTopoVonMisesKS(T E, T nu, T q, T design_stress, T ks_penalty)
      : q(q), design_stress(design_stress), ks_penalty(ks_penalty) {
    mu = 0.5 * E / (1.0 + nu);
    lambda = E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu));

    max_failure_index = 1.0;
    failure_index_integral = 1.0;
  }

  // Number of dimensions
  static const index_t dim = D;

  // Number of data dimensions
  static const index_t data_dim = 1;

  // Space for the finite-element data
  using DataSpace = typename IntegrandTopoLinearElasticity<T, D>::DataSpace;

  // Space for the element geometry
  using FiniteElementGeometry =
      typename IntegrandTopoLinearElasticity<T, D>::FiniteElementGeometry;

  // Finite element space
  using FiniteElementSpace =
      typename IntegrandTopoLinearElasticity<T, D>::FiniteElementSpace;

  // Mapping of the solution from the reference element to the physical element
  using SolutionMapping =
      typename IntegrandTopoLinearElasticity<T, D>::SolutionMapping;

  // Material parameters
  T mu;
  T lambda;

  // The RAMP penalty parameter
  T q;

  // Design stress value used in the constraint
  T design_stress;

  // The KS penalty parameter
  T ks_penalty;

  // Offset value - should be the maximum failure index anywhere in the domain
  T max_failure_index;

  // Integral of e^{ks_penalty * (failure_index - offset)}
  T failure_index_integral;

  /**
   * @brief Set the maximum failure index value (or approximate value)
   *
   * @param max_failure_index_ Maximum value of the failure index anywhere in
   * the domain
   */
  void set_max_failure_index(T max_failure_index_) {
    max_failure_index = max_failure_index_;
  }

  /**
   * @brief Evaluate the functional value based on the max failure index in the
   * domain and the failure index integral
   *
   * @param failure_index_integral_ Integral of the failure index
   * @return T The failure value
   */
  T evaluate_functional(T failure_index_integral_) {
    failure_index_integral = failure_index_integral_;
    return max_failure_index + log(failure_index_integral) / ks_penalty;
  }

  /**
   * @brief Compute the failure index at a quadrature point
   *
   * @param data The data at the quadrature point
   * @param geo The geometry at the quadrature point
   * @param s The solution at the quadurature point
   * @return T The integrand contribution
   */
  T max(const DataSpace& data, const FiniteElementGeometry& geo,
        const FiniteElementSpace& s) const {
    const Mat<T, dim, dim>& Ux = (s.template get<0>()).get_grad();
    SymMat<T, dim> E, S;
    T trS, trSS;

    MatLinearGreenStrain(Ux, E);
    SymmIsotropicConstitutive(mu, lambda, E, S);
    SymmTrace(S, trS);
    SymmSymmMultTrace(S, S, trSS);

    // Extract the design density value
    T rho = data[0];

    // Compute the penalty = (q + 1) * rho/(q * rho + 1)
    T penalty = (q + 1.0) * rho / (q * rho + 1.0);

    // von Mises^2 = 1.5 * tr(S * S) - 0.5 * tr(S)**2;
    T vm = std::sqrt(1.5 * trSS - 0.5 * trS * trS) / design_stress;

    // Compute the failure index
    T failure_index = penalty * vm;

    return failure_index;
  }

  /**
   * @brief Compute the integrand for this functional
   *
   * @param wdetJ The determinant of the Jacobian times the quadrature weight
   * @param data The data at the quadrature point
   * @param geo The geometry at the quadrature point
   * @param s The solution at the quadurature point
   * @return T The integrand contribution
   */
  T integrand(T wdetJ, const DataSpace& data, const FiniteElementGeometry& geo,
              const FiniteElementSpace& s) const {
    const Mat<T, dim, dim>& Ux = (s.template get<0>()).get_grad();
    SymMat<T, dim> E, S;
    T trS, trSS;

    MatLinearGreenStrain(Ux, E);
    SymmIsotropicConstitutive(mu, lambda, E, S);
    SymmTrace(S, trS);
    SymmSymmMultTrace(S, S, trSS);

    // Extract the design density value
    T rho = data[0];

    // Compute the penalty = (q + 1) * rho/(q * rho + 1)
    T penalty = (q + 1.0) * rho / (q * rho + 1.0);

    // von Mises^2 = 1.5 * tr(S * S) - 0.5 * tr(S)**2;
    T vm = std::sqrt(1.5 * trSS - 0.5 * trS * trS) / design_stress;

    // Compute the failure index
    T failure_index = penalty * vm;

    return wdetJ * exp(ks_penalty * (failure_index - max_failure_index));
  }

  /**
   * @brief Evaluate the weak form coefficients for linear elasticity
   *
   * @param wdetJ The quadrature weight times determinant of the Jacobian
   * @param data The data at the quadrature point
   * @param geo The geometry at the quadrature point
   * @param s The trial solution
   * @param coef Output weak form coefficients of the test space
   */
  void weak(T wdetJ, const DataSpace& data, const FiniteElementGeometry& geo,
            const FiniteElementSpace& s, FiniteElementSpace& coef) const {
    Mat<T, dim, dim> Ux0 = (s.template get<0>()).get_grad();
    Mat<T, dim, dim>& Uxb = (coef.template get<0>()).get_grad();
    SymMat<T, dim> E0, Eb;
    SymMat<T, dim> S0, Sb;

    ADMat<Mat<T, dim, dim>> Ux(Ux0, Uxb);
    ADMat<SymMat<T, dim>> E(E0, Eb);
    ADMat<SymMat<T, dim>> S(S0, Sb);
    ADScalar<T> trS, trSS;

    auto strain = MatLinearGreenStrain(Ux, E);
    auto cons = SymmIsotropicConstitutive(mu, lambda, E, S);
    auto trace1 = SymmTrace(S, trS);
    auto trace2 = SymmSymmMultTrace(S, S, trSS);

    // Extract the design density value
    T rho = data[0];

    // Compute the penalty = (q + 1) * x/(q * x + 1)
    T penalty = (q + 1.0) * rho / (q * rho + 1.0);

    // von Mises = 1.5 * tr(S * S) - 0.5 * tr(S)**2;
    T vm = std::sqrt(1.5 * trSS.value - 0.5 * trS.value * trS.value) /
           design_stress;

    // Compute the failure index
    T failure_index = penalty * vm;

    // Compute the exponential contribution
    T ks_exp = exp(ks_penalty * (failure_index - max_failure_index));

    T scale = 0.5 * wdetJ * penalty * ks_exp /
              (vm * design_stress * design_stress * failure_index_integral);

    trSS.bvalue = 1.5 * scale;
    trS.bvalue = -trS.value * scale;

    trace2.reverse();
    trace1.reverse();
    cons.reverse();
    strain.reverse();
  }

  /**
   * @brief Derivative of the integrand with respect to the data
   *
   * @param wdetJ The determinant of the Jacobian times the quadrature weight
   * @param data The data at the quadrature point
   * @param geo The geometry at the quadrature point
   * @param s The solution at the quadurature point
   * @param dfdx The output derivative value
   */
  void data_derivative(T wdetJ, const DataSpace& data,
                       const FiniteElementGeometry& geo,
                       const FiniteElementSpace& s, DataSpace& dfdx) const {
    const Mat<T, dim, dim>& Ux = (s.template get<0>()).get_grad();
    SymMat<T, dim> E, S;
    T trS, trSS;

    MatLinearGreenStrain(Ux, E);
    SymmIsotropicConstitutive(mu, lambda, E, S);
    SymmTrace(S, trS);
    SymmSymmMultTrace(S, S, trSS);

    // Extract the design density value
    T rho = data[0];

    // Compute the penalty = (q + 1) * x/(q * x + 1)
    T penalty = (q + 1.0) * rho / (q * rho + 1.0);

    // Compute the penalty = (q + 1) * x/(q * x + 1)
    T denom = (q * rho + 1.0) * (q * rho + 1.0);
    T dpenalty = (q + 1.0) / denom;

    // von Mises = 1.5 * tr(S * S) - 0.5 * tr(S)**2;
    T vm = std::sqrt(1.5 * trSS - 0.5 * trS * trS) / design_stress;

    // Compute the failure index
    T failure_index = penalty * vm;

    // Compute the exponential contribution
    T ks_exp = exp(ks_penalty * (failure_index - max_failure_index));

    T scale = wdetJ * vm * ks_exp / (failure_index_integral);

    dfdx[0] = scale * dpenalty;
  }
};

/**
 * @brief Apply surface traction and/or surface torque.
 */
template <typename T, index_t D>
class IntegrandTopoSurfaceTraction {
 public:
  IntegrandTopoSurfaceTraction(const T tx_[] = nullptr,
                               const T torx_[] = nullptr,
                               const T x0_[] = nullptr) {
    has_traction = false;
    has_torque = false;

    if (tx_) {
      for (index_t i = 0; i < dim; i++) {
        tx[i] = tx_[i];
      }
      has_traction = true;
    }

    if (torx_ && x0_) {
      for (index_t i = 0; i < dim; i++) {
        x0[i] = x0_[i];
        if constexpr (dim == 3) {
          torx[i] = torx_[i];
        }
      }
      if constexpr (dim == 2) {
        torx[0] = torx_[0];
      }
      has_torque = true;
    }
  }

  // Number of dimensions
  static const index_t dim = D;

  // Number of data dimensions
  static const index_t data_dim = 1;

  // Space for the finite-element data
  using DataSpace = FESpace<T, dim>;

  // Space for the element geometry
  using FiniteElementGeometry = FESpace<T, dim, H1Space<T, dim, dim - 1>>;

  // Finite element space
  using FiniteElementSpace = FESpace<T, dim, H1Space<T, dim, dim - 1>>;

  // Mapping of the solution from the reference element to the physical element
  using SolutionMapping = SurfaceMapping<T, dim>;

  T tx[dim];  // surface traction vector
  T torx[conditional_value<index_t, dim == 3, 3, 1>::value];  // surface torque
                                                              // vector
  T x0[dim];                                                  // torque origin
  bool has_traction;
  bool has_torque;

  /**
   * @brief Evaluate the weak form coefficients for linear elasticity
   *
   * @param wdetJ The quadrature weight times determinant of the Jacobian
   * @param data The data at the quadrature point
   * @param geo The geometry at the quadrature point
   * @param s The trial solution
   * @param coef Output weak form coefficients of the test space
   */
  KOKKOS_FUNCTION void weak(T wdetJ, const DataSpace& data,
                            const FiniteElementGeometry& geo,
                            const FiniteElementSpace& s,
                            FiniteElementSpace& coef) const {
    // Extract the solution
    Vec<T, dim>& U = (coef.template get<0>()).get_value();
    for (index_t i = 0; i < dim; i++) {
      U(i) = 0.0;
    }

    if (has_traction) {
      for (index_t i = 0; i < dim; i++) {
        U(i) -= wdetJ * tx[i];
      }
    }

    if (has_torque) {
      // Extract location
      const Vec<T, dim>& x = (geo.template get<0>()).get_value();

      if constexpr (dim == 2) {
        // Force at this point is (x - x0) cross torque
        U(0) += wdetJ * torx[0] * (x(1) - x0[1]);
        U(1) += -wdetJ * torx[0] * (x(0) - x0[0]);
      } else {  // dim == 3
        // Force at this point is (x - x0) cross torque
        U(0) += wdetJ * ((x(1) - x0[1]) * torx[2] - (x(2) - x0[2]) * torx[1]);
        U(1) += wdetJ * ((x(2) - x0[2]) * torx[0] - (x(0) - x0[0]) * torx[2]);
        U(2) += wdetJ * ((x(0) - x0[0]) * torx[1] - (x(1) - x0[1]) * torx[0]);
      }
    }
  }
};

}  // namespace A2D

#endif  // A2D_ELASTICITY_H
