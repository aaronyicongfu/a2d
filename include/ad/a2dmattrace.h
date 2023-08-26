#ifndef A2D_MAT_TRACE_H
#define A2D_MAT_TRACE_H

#include <type_traits>

#include "a2denum.h"
#include "a2dmat.h"
#include "a2dobjs.h"
#include "a2dscalar.h"

namespace A2D {

template <typename T, int M>
KOKKOS_FUNCTION T MatTraceCore(const T* A) {
  T trace = T(0.0);
  for (int i = 0; i < M; i++) {
    trace += A[0];
    A += M + 1;
  }
  return trace;
}

template <typename T, int M>
KOKKOS_FUNCTION void MatAddDiagCore(const T diag, T* A) {
  for (int i = 0; i < M; i++) {
    A[0] += diag;
    A += M + 1;
  }
}

template <typename T, int M>
KOKKOS_FUNCTION void MatTrace(Mat<T, M, M>& A, T& trace) {
  trace = MatTraceCore<T, M>(get_data(A));
}

template <typename T, int M, ADorder order, ADiffType adA>
class MatTraceExpr {
 private:
  using Atype = ADMatType<adA, order, Mat<T, M, M>>;
  using ScalarType = ADScalarType<ADiffType::ACTIVE, order, T>;

 public:
  KOKKOS_FUNCTION MatTraceExpr(Atype& A, ScalarType& tr) : A(A), tr(tr) {
    get_data(tr) = MatTraceCore<T, M>(get_data(A));
  }

  template <ADorder forder>
  KOKKOS_FUNCTION void forward() {
    static_assert(
        !(order == ADorder::FIRST and forder == ADorder::SECOND),
        "Can't perform second order forward with first order objects");
    constexpr ADseed seed = conditional_value<ADseed, forder == ADorder::FIRST,
                                              ADseed::b, ADseed::p>::value;
    if constexpr (adA == ADiffType::ACTIVE) {
      GetSeed<seed>::get_data(tr) =
          MatTraceCore<T, M>(GetSeed<seed>::get_data(A));
    }
  }

  KOKKOS_FUNCTION void reverse() {
    if constexpr (adA == ADiffType::ACTIVE) {
      MatAddDiagCore<T, M>(GetSeed<ADseed::b>::get_data(tr),
                           GetSeed<ADseed::b>::get_data(A));
    }
  }

  KOKKOS_FUNCTION void hreverse() {
    if constexpr (adA == ADiffType::ACTIVE) {
      MatAddDiagCore<T, M>(GetSeed<ADseed::h>::get_data(tr),
                           GetSeed<ADseed::h>::get_data(A));
    }
  }

 private:
  Atype& A;
  ScalarType& tr;
};

template <typename T, int M>
KOKKOS_FUNCTION auto MatTrace(ADMat<Mat<T, M, M>>& A, ADScalar<T>& tr) {
  return MatTraceExpr<T, M, ADorder::FIRST, ADiffType::ACTIVE>(A, tr);
}
template <typename T, int M>
KOKKOS_FUNCTION auto MatTrace(A2DMat<Mat<T, M, M>>& A, A2DScalar<T>& tr) {
  return MatTraceExpr<T, M, ADorder::SECOND, ADiffType::ACTIVE>(A, tr);
}

namespace Test {

template <typename T, int N>
class MatTraceTest : public A2DTest<T, T, Mat<T, N, N>> {
 public:
  using Input = VarTuple<T, Mat<T, N, N>>;
  using Output = VarTuple<T, T>;

  // Assemble a string to describe the test
  std::string name() {
    std::stringstream s;
    s << "MatTrace<" << N << "," << N << ">";
    return s.str();
  }

  // Evaluate the matrix-matrix product
  Output eval(const Input& x) {
    T trace;
    Mat<T, N, N> A;
    x.get_values(A);
    MatTrace(A, trace);
    return MakeVarTuple<T>(trace);
  }

  // Compute the derivative
  void deriv(const Output& seed, const Input& x, Input& g) {
    ADScalar<T> trace;
    Mat<T, N, N> A0, Ab;
    ADMat<Mat<T, N, N>> A(A0, Ab);

    x.get_values(A0);
    auto op = MatTrace(A, trace);
    auto stack = MakeStack(op);
    seed.get_values(trace.bvalue);
    stack.reverse();
    g.set_values(Ab);
  }

  // Compute the second-derivative
  void hprod(const Output& seed, const Output& hval, const Input& x,
             const Input& p, Input& h) {
    A2DScalar<T> trace;
    A2DMat<Mat<T, N, N>> A;
    x.get_values(A.value());
    p.get_values(A.pvalue());

    auto op = MatTrace(A, trace);
    auto stack = MakeStack(op);

    seed.get_values(trace.bvalue);
    hval.get_values(trace.hvalue);
    stack.reverse();
    stack.hforward();
    stack.hreverse();
    h.set_values(A.hvalue());
  }
};

void MatTraceTestAll() {
  using Tc = std::complex<double>;
  MatTraceTest<Tc, 2> test1;
  Run(test1);
  MatTraceTest<Tc, 4> test2;
  Run(test2);
}

}  // namespace Test

}  // namespace A2D

#endif  // A2D_MAT_TRACE_H