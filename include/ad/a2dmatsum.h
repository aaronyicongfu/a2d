#ifndef A2D_MAT_ADD_H
#define A2D_MAT_ADD_H

#include <type_traits>

#include "a2denum.h"
#include "a2dmat.h"
#include "a2dobjs.h"
#include "a2dstack.h"
#include "a2dtest.h"
#include "ad/core/a2dveccore.h"

namespace A2D {

template <typename T, int N, int M>
KOKKOS_FUNCTION void MatSum(const Mat<T, N, M> &A, const Mat<T, N, M> &B,
                            Mat<T, N, M> &C) {
  VecSumCore<T, N * M>(get_data(A), get_data(B), get_data(C));
}

template <typename T, int N, int M>
KOKKOS_FUNCTION void MatSum(const T alpha, const Mat<T, N, M> &A, const T beta,
                            const Mat<T, N, M> &B, Mat<T, N, M> &C) {
  VecSumCore<T, N * M>(alpha, get_data(A), beta, get_data(B), get_data(C));
}

template <typename T, int N>
KOKKOS_FUNCTION void MatSum(const SymMat<T, N> &A, const SymMat<T, N> &B,
                            SymMat<T, N> &C) {
  VecSumCore<T, (N * (N + 1)) / 2>(get_data(A), get_data(B), get_data(C));
}

template <typename T, int N>
KOKKOS_FUNCTION void MatSum(const T alpha, const SymMat<T, N> &A, const T beta,
                            const SymMat<T, N> &B, SymMat<T, N> &C) {
  VecSumCore<T, (N * (N + 1)) / 2>(alpha, get_data(A), beta, get_data(B),
                                   get_data(C));
}

template <typename T, int N, int M, ADorder order, ADiffType adA, ADiffType adB,
          MatSymType sym>
class MatSumExpr {
 public:
  using Atype =
      ADMatType<adA, order,
                typename std::conditional<sym == MatSymType::NORMAL,
                                          Mat<T, N, M>, SymMat<T, N>>::type>;
  using Btype =
      ADMatType<adB, order,
                typename std::conditional<sym == MatSymType::NORMAL,
                                          Mat<T, N, M>, SymMat<T, N>>::type>;
  using Ctype =
      ADMatType<ADiffType::ACTIVE, order,
                typename std::conditional<sym == MatSymType::NORMAL,
                                          Mat<T, N, M>, SymMat<T, N>>::type>;
  static const int size = conditional_value<int, sym == MatSymType::NORMAL,
                                            N * M, (N * (N + 1) / 2)>::value;

  KOKKOS_FUNCTION
  MatSumExpr(Atype &A, Btype &B, Ctype &C) : A(A), B(B), C(C) {}

  KOKKOS_FUNCTION void eval() {
    VecSumCore<T, size>(get_data(A), get_data(B), get_data(C));
  }

  template <ADorder forder>
  KOKKOS_FUNCTION void forward() {
    constexpr ADseed seed = conditional_value<ADseed, forder == ADorder::FIRST,
                                              ADseed::b, ADseed::p>::value;
    if constexpr (adA == ADiffType::ACTIVE && adB == ADiffType::ACTIVE) {
      VecSumCore<T, size>(GetSeed<seed>::get_data(A),
                          GetSeed<seed>::get_data(B),
                          GetSeed<seed>::get_data(C));
    } else if constexpr (adA == ADiffType::ACTIVE) {
      VecCopyCore<T, size>(GetSeed<seed>::get_data(A),
                           GetSeed<seed>::get_data(C));
    } else if constexpr (adB == ADiffType::ACTIVE) {
      VecCopyCore<T, size>(GetSeed<seed>::get_data(B),
                           GetSeed<seed>::get_data(C));
    }
  }

  KOKKOS_FUNCTION void reverse() {
    constexpr ADseed seed = ADseed::b;
    if constexpr (adA == ADiffType::ACTIVE) {
      VecAddCore<T, size>(GetSeed<seed>::get_data(C),
                          GetSeed<seed>::get_data(A));
    }
    if constexpr (adB == ADiffType::ACTIVE) {
      VecAddCore<T, size>(GetSeed<seed>::get_data(C),
                          GetSeed<seed>::get_data(B));
    }
  }

  KOKKOS_FUNCTION void hreverse() {
    constexpr ADseed seed = ADseed::h;
    if constexpr (adA == ADiffType::ACTIVE) {
      VecAddCore<T, size>(GetSeed<seed>::get_data(C),
                          GetSeed<seed>::get_data(A));
    }
    if constexpr (adB == ADiffType::ACTIVE) {
      VecAddCore<T, size>(GetSeed<seed>::get_data(C),
                          GetSeed<seed>::get_data(B));
    }
  }

  Atype &A;
  Btype &B;
  Ctype &C;
};

// First-order AD
template <typename T, int N, int M>
KOKKOS_FUNCTION auto MatSum(ADMat<Mat<T, N, M>> &A, ADMat<Mat<T, N, M>> &B,
                            ADMat<Mat<T, N, M>> &C) {
  return MatSumExpr<T, N, M, ADorder::FIRST, ADiffType::ACTIVE,
                    ADiffType::ACTIVE, MatSymType::NORMAL>(A, B, C);
}

template <typename T, int N, int M>
KOKKOS_FUNCTION auto MatSum(Mat<T, N, M> &A, ADMat<Mat<T, N, M>> &B,
                            ADMat<Mat<T, N, M>> &C) {
  return MatSumExpr<T, N, M, ADorder::FIRST, ADiffType::PASSIVE,
                    ADiffType::ACTIVE, MatSymType::NORMAL>(A, B, C);
}

template <typename T, int N, int M>
KOKKOS_FUNCTION auto MatSum(ADMat<Mat<T, N, M>> &A, Mat<T, N, M> &B,
                            ADMat<Mat<T, N, M>> &C) {
  return MatSumExpr<T, N, M, ADorder::FIRST, ADiffType::ACTIVE,
                    ADiffType::PASSIVE, MatSymType::NORMAL>(A, B, C);
}

template <typename T, int N, int M>
KOKKOS_FUNCTION auto MatSum(ADMat<SymMat<T, N>> &A, ADMat<SymMat<T, N>> &B,
                            ADMat<SymMat<T, N>> &C) {
  return MatSumExpr<T, N, M, ADorder::FIRST, ADiffType::ACTIVE,
                    ADiffType::ACTIVE, MatSymType::SYMMETRIC>(A, B, C);
}

template <typename T, int N, int M>
KOKKOS_FUNCTION auto MatSum(SymMat<T, N> &A, ADMat<SymMat<T, N>> &B,
                            ADMat<SymMat<T, N>> &C) {
  return MatSumExpr<T, N, M, ADorder::FIRST, ADiffType::PASSIVE,
                    ADiffType::ACTIVE, MatSymType::SYMMETRIC>(A, B, C);
}

template <typename T, int N, int M>
KOKKOS_FUNCTION auto MatSum(ADMat<SymMat<T, N>> &A, SymMat<T, N> &B,
                            ADMat<SymMat<T, N>> &C) {
  return MatSumExpr<T, N, M, ADorder::FIRST, ADiffType::ACTIVE,
                    ADiffType::PASSIVE, MatSymType::SYMMETRIC>(A, B, C);
}

// Second-order AD
template <typename T, int N, int M>
KOKKOS_FUNCTION auto MatSum(A2DMat<Mat<T, N, M>> &A, A2DMat<Mat<T, N, M>> &B,
                            A2DMat<Mat<T, N, M>> &C) {
  return MatSumExpr<T, N, M, ADorder::SECOND, ADiffType::ACTIVE,
                    ADiffType::ACTIVE, MatSymType::NORMAL>(A, B, C);
}

template <typename T, int N, int M>
KOKKOS_FUNCTION auto MatSum(Mat<T, N, M> &A, A2DMat<Mat<T, N, M>> &B,
                            A2DMat<Mat<T, N, M>> &C) {
  return MatSumExpr<T, N, M, ADorder::SECOND, ADiffType::PASSIVE,
                    ADiffType::ACTIVE, MatSymType::NORMAL>(A, B, C);
}

template <typename T, int N, int M>
KOKKOS_FUNCTION auto MatSum(A2DMat<Mat<T, N, M>> &A, Mat<T, N, M> &B,
                            A2DMat<Mat<T, N, M>> &C) {
  return MatSumExpr<T, N, M, ADorder::SECOND, ADiffType::ACTIVE,
                    ADiffType::PASSIVE, MatSymType::NORMAL>(A, B, C);
}

template <typename T, int N, int M>
KOKKOS_FUNCTION auto MatSum(A2DMat<SymMat<T, N>> &A, A2DMat<SymMat<T, N>> &B,
                            A2DMat<SymMat<T, N>> &C) {
  return MatSumExpr<T, N, M, ADorder::SECOND, ADiffType::ACTIVE,
                    ADiffType::ACTIVE, MatSymType::SYMMETRIC>(A, B, C);
}

template <typename T, int N, int M>
KOKKOS_FUNCTION auto MatSum(SymMat<T, N> &A, A2DMat<SymMat<T, N>> &B,
                            A2DMat<SymMat<T, N>> &C) {
  return MatSumExpr<T, N, M, ADorder::SECOND, ADiffType::PASSIVE,
                    ADiffType::ACTIVE, MatSymType::SYMMETRIC>(A, B, C);
}

template <typename T, int N, int M>
KOKKOS_FUNCTION auto MatSum(A2DMat<SymMat<T, N>> &A, SymMat<T, N> &B,
                            A2DMat<SymMat<T, N>> &C) {
  return MatSumExpr<T, N, M, ADorder::SECOND, ADiffType::ACTIVE,
                    ADiffType::PASSIVE, MatSymType::SYMMETRIC>(A, B, C);
}

template <typename T, int N, int M, ADorder order, ADiffType ada, ADiffType adA,
          ADiffType adb, ADiffType adB, MatSymType sym>
class MatSumScaleExpr {
 public:
  using atype = ADScalarInputType<ada, order, T>;
  using Atype =
      ADMatType<adA, order,
                typename std::conditional<sym == MatSymType::NORMAL,
                                          Mat<T, N, M>, SymMat<T, N>>::type>;
  using btype = ADScalarInputType<adb, order, T>;
  using Btype =
      ADMatType<adB, order,
                typename std::conditional<sym == MatSymType::NORMAL,
                                          Mat<T, N, M>, SymMat<T, N>>::type>;

  using Ctype =
      ADMatType<ADiffType::ACTIVE, order,
                typename std::conditional<sym == MatSymType::NORMAL,
                                          Mat<T, N, M>, SymMat<T, N>>::type>;

  static const int size = conditional_value<int, sym == MatSymType::NORMAL,
                                            N * M, (N * (N + 1) / 2)>::value;

  KOKKOS_FUNCTION
  MatSumScaleExpr(atype alpha, Atype &A, btype beta, Btype &B, Ctype &C)
      : alpha(alpha), A(A), beta(beta), B(B), C(C) {}

  KOKKOS_FUNCTION void eval() {
    VecSumCore<T, size>(get_data(alpha), get_data(A), get_data(beta),
                        get_data(B), get_data(C));
  }

  template <ADorder forder>
  KOKKOS_FUNCTION void forward() {
    constexpr ADseed seed = conditional_value<ADseed, forder == ADorder::FIRST,
                                              ADseed::b, ADseed::p>::value;

    if constexpr (ada == ADiffType::ACTIVE && adA == ADiffType::ACTIVE &&
                  adb == ADiffType::ACTIVE && adB == ADiffType::ACTIVE) {
      VecSumCore<T, size>(get_data(alpha), GetSeed<seed>::get_data(A),
                          get_data(beta), GetSeed<seed>::get_data(B),
                          GetSeed<seed>::get_data(C));
      VecSumCoreAdd<T, size>(GetSeed<seed>::get_data(alpha), get_data(A),
                             GetSeed<seed>::get_data(beta), get_data(B),
                             GetSeed<seed>::get_data(C));
    } else {
      VecZeroCore<T, size>(GetSeed<seed>::get_data(C));
      if constexpr (adA == ADiffType::ACTIVE) {
        VecAddCore<T, size>(get_data(alpha), GetSeed<seed>::get_data(A),
                            GetSeed<seed>::get_data(C));
      }
      if constexpr (adB == ADiffType::ACTIVE) {
        VecAddCore<T, size>(get_data(beta), GetSeed<seed>::get_data(B),
                            GetSeed<seed>::get_data(C));
      }
      if constexpr (ada == ADiffType::ACTIVE) {
        VecAddCore<T, size>(GetSeed<seed>::get_data(alpha), get_data(A),
                            GetSeed<seed>::get_data(C));
      }
      if constexpr (adb == ADiffType::ACTIVE) {
        VecAddCore<T, size>(GetSeed<seed>::get_data(beta), get_data(B),
                            GetSeed<seed>::get_data(C));
      }
    }
  }

  KOKKOS_FUNCTION void reverse() {
    constexpr ADseed seed = ADseed::b;
    if constexpr (adA == ADiffType::ACTIVE) {
      VecAddCore<T, size>(get_data(alpha), GetSeed<seed>::get_data(C),
                          GetSeed<seed>::get_data(A));
    }
    if constexpr (adB == ADiffType::ACTIVE) {
      VecAddCore<T, size>(get_data(beta), GetSeed<seed>::get_data(C),
                          GetSeed<seed>::get_data(B));
    }
    if constexpr (ada == ADiffType::ACTIVE) {
      GetSeed<seed>::get_data(alpha) +=
          VecDotCore<T, size>(GetSeed<seed>::get_data(C), get_data(A));
    }
    if constexpr (adb == ADiffType::ACTIVE) {
      GetSeed<seed>::get_data(beta) +=
          VecDotCore<T, size>(GetSeed<seed>::get_data(C), get_data(B));
    }
  }

  KOKKOS_FUNCTION void hreverse() {
    constexpr ADseed seed = ADseed::h;
    if constexpr (adA == ADiffType::ACTIVE) {
      VecAddCore<T, size>(get_data(alpha), GetSeed<seed>::get_data(C),
                          GetSeed<seed>::get_data(A));
    }
    if constexpr (adB == ADiffType::ACTIVE) {
      VecAddCore<T, size>(get_data(beta), GetSeed<seed>::get_data(C),
                          GetSeed<seed>::get_data(B));
    }
    if constexpr (ada == ADiffType::ACTIVE) {
      GetSeed<seed>::get_data(alpha) +=
          VecDotCore<T, size>(GetSeed<seed>::get_data(C), get_data(A));
    }
    if constexpr (adb == ADiffType::ACTIVE) {
      GetSeed<seed>::get_data(beta) +=
          VecDotCore<T, size>(GetSeed<seed>::get_data(C), get_data(B));
    }
    if constexpr (adA == ADiffType::ACTIVE && ada == ADiffType::ACTIVE) {
      GetSeed<seed>::get_data(alpha) += VecDotCore<T, size>(
          GetSeed<ADseed::b>::get_data(C), GetSeed<ADseed::p>::get_data(A));
      VecAddCore<T, size>(GetSeed<ADseed::p>::get_data(alpha),
                          GetSeed<ADseed::b>::get_data(C),
                          GetSeed<seed>::get_data(A));
    }
    if constexpr (adB == ADiffType::ACTIVE && adb == ADiffType::ACTIVE) {
      GetSeed<seed>::get_data(beta) += VecDotCore<T, size>(
          GetSeed<ADseed::b>::get_data(C), GetSeed<ADseed::p>::get_data(B));
      VecAddCore<T, size>(GetSeed<ADseed::p>::get_data(beta),
                          GetSeed<ADseed::b>::get_data(C),
                          GetSeed<seed>::get_data(B));
    }
  }

  atype alpha;
  Atype &A;
  btype beta;
  Btype &B;
  Ctype &C;
};

// First-order AD
template <typename T, int N, int M>
KOKKOS_FUNCTION auto MatSum(ADScalar<T> &alpha, ADMat<Mat<T, N, M>> &A,
                            ADScalar<T> &beta, ADMat<Mat<T, N, M>> &B,
                            ADMat<Mat<T, N, M>> &C) {
  return MatSumScaleExpr<T, N, M, ADorder::FIRST, ADiffType::ACTIVE,
                         ADiffType::ACTIVE, ADiffType::ACTIVE,
                         ADiffType::ACTIVE, MatSymType::NORMAL>(alpha, A, beta,
                                                                B, C);
}

template <typename T, int N, int M>
KOKKOS_FUNCTION auto MatSum(const T alpha, ADMat<Mat<T, N, M>> &A, const T beta,
                            ADMat<Mat<T, N, M>> &B, ADMat<Mat<T, N, M>> &C) {
  return MatSumScaleExpr<T, N, M, ADorder::FIRST, ADiffType::PASSIVE,
                         ADiffType::ACTIVE, ADiffType::PASSIVE,
                         ADiffType::ACTIVE, MatSymType::NORMAL>(alpha, A, beta,
                                                                B, C);
}

template <typename T, int N>
KOKKOS_FUNCTION auto MatSum(ADScalar<T> &alpha, ADMat<SymMat<T, N>> &A,
                            ADScalar<T> &beta, ADMat<SymMat<T, N>> &B,
                            ADMat<SymMat<T, N>> &C) {
  return MatSumScaleExpr<T, N, N, ADorder::FIRST, ADiffType::ACTIVE,
                         ADiffType::ACTIVE, ADiffType::ACTIVE,
                         ADiffType::ACTIVE, MatSymType::SYMMETRIC>(alpha, A,
                                                                   beta, B, C);
}

template <typename T, int N>
KOKKOS_FUNCTION auto MatSum(const T alpha, ADMat<SymMat<T, N>> &A, const T beta,
                            ADMat<SymMat<T, N>> &B, ADMat<SymMat<T, N>> &C) {
  return MatSumScaleExpr<T, N, N, ADorder::FIRST, ADiffType::PASSIVE,
                         ADiffType::ACTIVE, ADiffType::PASSIVE,
                         ADiffType::ACTIVE, MatSymType::SYMMETRIC>(alpha, A,
                                                                   beta, B, C);
}

// Second-order AD
template <typename T, int N, int M>
KOKKOS_FUNCTION auto MatSum(A2DScalar<T> &alpha, A2DMat<Mat<T, N, M>> &A,
                            A2DScalar<T> &beta, A2DMat<Mat<T, N, M>> &B,
                            A2DMat<Mat<T, N, M>> &C) {
  return MatSumScaleExpr<T, N, M, ADorder::SECOND, ADiffType::ACTIVE,
                         ADiffType::ACTIVE, ADiffType::ACTIVE,
                         ADiffType::ACTIVE, MatSymType::NORMAL>(alpha, A, beta,
                                                                B, C);
}

template <typename T, int N, int M>
KOKKOS_FUNCTION auto MatSum(const T alpha, A2DMat<Mat<T, N, M>> &A,
                            const T beta, A2DMat<Mat<T, N, M>> &B,
                            A2DMat<Mat<T, N, M>> &C) {
  return MatSumScaleExpr<T, N, M, ADorder::SECOND, ADiffType::PASSIVE,
                         ADiffType::ACTIVE, ADiffType::PASSIVE,
                         ADiffType::ACTIVE, MatSymType::NORMAL>(alpha, A, beta,
                                                                B, C);
}

template <typename T, int N>
KOKKOS_FUNCTION auto MatSum(A2DScalar<T> &alpha, A2DMat<SymMat<T, N>> &A,
                            A2DScalar<T> &beta, A2DMat<SymMat<T, N>> &B,
                            A2DMat<SymMat<T, N>> &C) {
  return MatSumScaleExpr<T, N, N, ADorder::SECOND, ADiffType::ACTIVE,
                         ADiffType::ACTIVE, ADiffType::ACTIVE,
                         ADiffType::ACTIVE, MatSymType::SYMMETRIC>(alpha, A,
                                                                   beta, B, C);
}

template <typename T, int N>
KOKKOS_FUNCTION auto MatSum(const T alpha, A2DMat<SymMat<T, N>> &A,
                            const T beta, A2DMat<SymMat<T, N>> &B,
                            A2DMat<SymMat<T, N>> &C) {
  return MatSumScaleExpr<T, N, N, ADorder::SECOND, ADiffType::PASSIVE,
                         ADiffType::ACTIVE, ADiffType::PASSIVE,
                         ADiffType::ACTIVE, MatSymType::SYMMETRIC>(alpha, A,
                                                                   beta, B, C);
}

namespace Test {

template <typename T, int N, int M>
class MatSumTest : public A2DTest<T, Mat<T, N, M>, Mat<T, N, M>, Mat<T, N, M>> {
 public:
  using Input = VarTuple<T, Mat<T, N, M>, Mat<T, N, M>>;
  using Output = VarTuple<T, Mat<T, N, M>>;

  // Assemble a string to describe the test
  std::string name() {
    std::stringstream s;
    s << "MatSum<" << N << "," << M << ">";
    return s.str();
  }

  // Evaluate the matrix-matrix product
  Output eval(const Input &x) {
    Mat<T, N, M> A, B, C;
    x.get_values(A, B);
    MatSum(A, B, C);
    return MakeVarTuple<T>(C);
  }

  // Compute the derivative
  void deriv(const Output &seed, const Input &x, Input &g) {
    Mat<T, N, M> A0, Ab, B0, Bb, C0, Cb;
    ADMat<Mat<T, N, M>> A(A0, Ab), B(B0, Bb), C(C0, Cb);
    x.get_values(A0, B0);
    auto op = MatSum(A, B, C);
    auto stack = MakeStack(op);
    seed.get_values(Cb);
    stack.reverse();
    g.set_values(Ab, Bb);
  }

  // Compute the second-derivative
  void hprod(const Output &seed, const Output &hval, const Input &x,
             const Input &p, Input &h) {
    A2DMat<Mat<T, N, M>> A, B, C;
    x.get_values(A.value(), B.value());
    p.get_values(A.pvalue(), B.pvalue());
    auto op = MatSum(A, B, C);
    auto stack = MakeStack(op);
    seed.get_values(C.bvalue());
    hval.get_values(C.hvalue());
    stack.reverse();
    stack.hforward();
    stack.hreverse();
    h.set_values(A.hvalue(), B.hvalue());
  }
};

template <typename T, int N, int M>
class MatSumScaleTest
    : public A2DTest<T, Mat<T, N, M>, T, Mat<T, N, M>, T, Mat<T, N, M>> {
 public:
  using Input = VarTuple<T, T, Mat<T, N, M>, T, Mat<T, N, M>>;
  using Output = VarTuple<T, Mat<T, N, M>>;

  // Assemble a string to describe the test
  std::string name() {
    std::stringstream s;
    s << "MatSum<" << N << "," << M << ">";
    return s.str();
  }

  // Evaluate the matrix-matrix product
  Output eval(const Input &x) {
    T alpha, beta;
    Mat<T, N, M> A, B, C;
    x.get_values(alpha, A, beta, B);
    MatSum(alpha, A, beta, B, C);
    return MakeVarTuple<T>(C);
  }

  // Compute the derivative
  void deriv(const Output &seed, const Input &x, Input &g) {
    ADScalar<T> alpha, beta;
    Mat<T, N, M> A0, Ab, B0, Bb, C0, Cb;
    ADMat<Mat<T, N, M>> A(A0, Ab), B(B0, Bb), C(C0, Cb);
    x.get_values(alpha.value, A0, beta.value, B0);
    auto op = MatSum(alpha, A, beta, B, C);
    auto stack = MakeStack(op);
    seed.get_values(Cb);
    stack.reverse();
    g.set_values(alpha.bvalue, Ab, beta.bvalue, Bb);
  }

  // Compute the second-derivative
  void hprod(const Output &seed, const Output &hval, const Input &x,
             const Input &p, Input &h) {
    A2DScalar<T> alpha, beta;
    A2DMat<Mat<T, N, M>> A, B, C;
    x.get_values(alpha.value, A.value(), beta.value, B.value());
    p.get_values(alpha.pvalue, A.pvalue(), beta.pvalue, B.pvalue());
    auto op = MatSum(alpha, A, beta, B, C);
    auto stack = MakeStack(op);
    seed.get_values(C.bvalue());
    hval.get_values(C.hvalue());
    stack.reverse();
    stack.hforward();
    stack.hreverse();
    h.set_values(alpha.hvalue, A.hvalue(), beta.hvalue, B.hvalue());
  }
};

bool MatSumTestAll(bool component = false, bool write_output = true) {
  using Tc = std::complex<double>;

  bool passed = true;
  MatSumTest<Tc, 3, 4> test1;
  passed = passed && Run(test1, component, write_output);
  MatSumTest<Tc, 5, 3> test2;
  passed = passed && Run(test2, component, write_output);

  MatSumScaleTest<Tc, 3, 4> test3;
  passed = passed && Run(test3, component, write_output);

  return passed;
}

}  // namespace Test

}  // namespace A2D

#endif  // A2D_MAT_ADD_H