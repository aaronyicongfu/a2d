#ifndef A2D_MAT_OPS_2D_H
#define A2D_MAT_OPS_2D_H

#include <stdlib.h>

#include "a2dmatcore2d.h"
#include "a2dobjs.h"
#include "a2dtypes.h"
#include "ad/a2dmat.h"
#include "ad/a2dscalar.h"
#include "ad/a2dvec.h"

namespace A2D {

// MatMatMult
template <typename ScalarType, bool AT = false, bool BT = false>
KOKKOS_FUNCTION void MatMatMult(const Mat<ScalarType, 2, 2>& A,
                                const Mat<ScalarType, 2, 2>& B,
                                Mat<ScalarType, 2, 2>& C) {
  if (AT && BT) {
    MatTrans2x2MatTransMultCore(A, B, C);
  } else if (AT) {
    MatTrans2x2MatMultCore(A, B, C);
  } else if (BT) {
    Mat2x2MatTransMultCore(A, B, C);
  } else {
    Mat2x2MatMultCore(A, B, C);
  }
}

template <typename ScalarType, bool AT = false, bool BT = false>
class ADMat2x2MatMultExpr
    : public ADExpression<ADMat2x2MatMultExpr<ScalarType, AT, BT>> {
 public:
  KOKKOS_FUNCTION ADMat2x2MatMultExpr(ADMat<Mat<ScalarType, 2, 2>>& AObj,
                                      ADMat<Mat<ScalarType, 2, 2>>& BObj,
                                      ADMat<Mat<ScalarType, 2, 2>>& CObj)
      : AObj(AObj), BObj(BObj), CObj(CObj) {
    const Mat<ScalarType, 2, 2>& A = AObj.value();
    const Mat<ScalarType, 2, 2>& B = BObj.value();
    Mat<ScalarType, 2, 2>& C = CObj.value();

    if (AT && BT) {
      MatTrans2x2MatTransMultCore(A, B, C);
    } else if (AT) {
      MatTrans2x2MatMultCore(A, B, C);
    } else if (BT) {
      Mat2x2MatTransMultCore(A, B, C);
    } else {
      Mat2x2MatMultCore(A, B, C);
    }
  }

  KOKKOS_FUNCTION void forward() {
    const Mat<ScalarType, 2, 2>& A = AObj.value();
    const Mat<ScalarType, 2, 2>& Ab = AObj.bvalue();
    const Mat<ScalarType, 2, 2>& B = BObj.value();
    const Mat<ScalarType, 2, 2>& Bb = BObj.bvalue();
    Mat<ScalarType, 2, 2>& Cb = CObj.bvalue();

    if (AT && BT) {
      MatTrans2x2MatTransMultCore(Ab, B, Cb);
      MatTrans2x2MatTransMultAddCore(A, Bb, Cb);
    } else if (AT) {
      MatTrans2x2MatMultCore(Ab, B, Cb);
      MatTrans2x2MatMultAddCore(A, Bb, Cb);
    } else if (BT) {
      Mat2x2MatTransMultCore(Ab, B, Cb);
      Mat2x2MatTransMultAddCore(A, Bb, Cb);
    } else {
      Mat2x2MatMultCore(Ab, B, Cb);
      Mat2x2MatMultAddCore(A, Bb, Cb);
    }
  }

  KOKKOS_FUNCTION void reverse() {
    const Mat<ScalarType, 2, 2>& A = AObj.value();
    Mat<ScalarType, 2, 2>& Ab = AObj.bvalue();
    const Mat<ScalarType, 2, 2>& B = BObj.value();
    Mat<ScalarType, 2, 2>& Bb = BObj.bvalue();
    const Mat<ScalarType, 2, 2>& Cb = CObj.bvalue();

    if (AT && BT) {
      MatTrans2x2MatTransMultAddCore(B, Cb, Ab);
      MatTrans2x2MatTransMultAddCore(Cb, A, Bb);
    } else if (AT) {
      Mat2x2MatTransMultAddCore(B, Cb, Ab);
      Mat2x2MatMultAddCore(A, Cb, Bb);
    } else if (BT) {
      Mat2x2MatMultAddCore(Cb, B, Ab);
      MatTrans2x2MatMultAddCore(Cb, A, Bb);
    } else {
      Mat2x2MatTransMultAddCore(Cb, B, Ab);
      MatTrans2x2MatMultAddCore(A, Cb, Bb);
    }
  }

  ADMat<Mat<ScalarType, 2, 2>>& AObj;
  ADMat<Mat<ScalarType, 2, 2>>& BObj;
  ADMat<Mat<ScalarType, 2, 2>>& CObj;
};

template <typename ScalarType, bool AT = false, bool BT = false>
KOKKOS_FUNCTION ADMat2x2MatMultExpr<ScalarType, AT, BT> MatMatMult(
    ADMat<Mat<ScalarType, 2, 2>>& AObj, ADMat<Mat<ScalarType, 2, 2>>& BObj,
    ADMat<Mat<ScalarType, 2, 2>>& CObj) {
  return ADMat2x2MatMultExpr<ScalarType, AT, BT>(AObj, BObj, CObj);
}

template <typename ScalarType, bool AT = false, bool BT = false>
class A2DMat2x2MatMultExpr
    : public A2DExpression<A2DMat2x2MatMultExpr<ScalarType, AT, BT>> {
 public:
  KOKKOS_FUNCTION A2DMat2x2MatMultExpr(A2DMat<Mat<ScalarType, 2, 2>>& AObj,
                                       A2DMat<Mat<ScalarType, 2, 2>>& BObj,
                                       A2DMat<Mat<ScalarType, 2, 2>>& CObj)
      : AObj(AObj), BObj(BObj), CObj(CObj) {
    const Mat<ScalarType, 2, 2>& A = AObj.value();
    const Mat<ScalarType, 2, 2>& B = BObj.value();
    Mat<ScalarType, 2, 2>& C = CObj.value();

    if (AT && BT) {
      MatTrans2x2MatTransMultCore(A, B, C);
    } else if (AT) {
      MatTrans2x2MatMultCore(A, B, C);
    } else if (BT) {
      Mat2x2MatTransMultCore(A, B, C);
    } else {
      Mat2x2MatMultCore(A, B, C);
    }
  }

  KOKKOS_FUNCTION void reverse() {
    const Mat<ScalarType, 2, 2>& A = AObj.value();
    Mat<ScalarType, 2, 2>& Ab = AObj.bvalue();
    const Mat<ScalarType, 2, 2>& B = BObj.value();
    Mat<ScalarType, 2, 2>& Bb = BObj.bvalue();
    const Mat<ScalarType, 2, 2>& Cb = CObj.bvalue();

    if (AT && BT) {
      MatTrans2x2MatTransMultAddCore(B, Cb, Ab);
      MatTrans2x2MatTransMultAddCore(Cb, A, Bb);
    } else if (AT) {
      Mat2x2MatTransMultAddCore(B, Cb, Ab);
      Mat2x2MatMultAddCore(A, Cb, Bb);
    } else if (BT) {
      Mat2x2MatMultAddCore(Cb, B, Ab);
      MatTrans2x2MatMultAddCore(Cb, A, Bb);
    } else {
      Mat2x2MatTransMultAddCore(Cb, B, Ab);
      MatTrans2x2MatMultAddCore(A, Cb, Bb);
    }
  }

  KOKKOS_FUNCTION void hforward() {
    const Mat<ScalarType, 2, 2>& A = AObj.value();
    const Mat<ScalarType, 2, 2>& B = BObj.value();

    const Mat<ScalarType, 2, 2>& Ap = AObj.pvalue();
    const Mat<ScalarType, 2, 2>& Bp = BObj.pvalue();
    Mat<ScalarType, 2, 2>& Cp = CObj.pvalue();

    if (AT && BT) {
      MatTrans2x2MatTransMultCore(Ap, B, Cp);
      MatTrans2x2MatTransMultAddCore(A, Bp, Cp);
    } else if (AT) {
      MatTrans2x2MatMultCore(Ap, B, Cp);
      MatTrans2x2MatMultAddCore(A, Bp, Cp);
    } else if (BT) {
      Mat2x2MatTransMultCore(Ap, B, Cp);
      Mat2x2MatTransMultAddCore(A, Bp, Cp);
    } else {
      Mat2x2MatMultCore(Ap, B, Cp);
      Mat2x2MatMultAddCore(A, Bp, Cp);
    }
  }

  KOKKOS_FUNCTION void hreverse() {
    const Mat<ScalarType, 2, 2>& A = AObj.value();
    const Mat<ScalarType, 2, 2>& B = BObj.value();

    Mat<ScalarType, 2, 2>& Ah = AObj.hvalue();
    Mat<ScalarType, 2, 2>& Bh = BObj.hvalue();
    const Mat<ScalarType, 2, 2>& Ch = CObj.hvalue();
    const Mat<ScalarType, 2, 2>& Cb = CObj.bvalue();
    const Mat<ScalarType, 2, 2>& Ap = AObj.pvalue();
    const Mat<ScalarType, 2, 2>& Bp = BObj.pvalue();

    if (AT && BT) {
      MatTrans2x2MatTransMultAddCore(B, Ch, Ah);
      MatTrans2x2MatTransMultAddCore(Ch, A, Bh);

      for (int ii = 0; ii < 2; ii++) {
        for (int jj = 0; jj < 2; jj++) {
          for (int kk = 0; kk < 2; kk++) {
            Ah(jj, ii) += Cb(ii, kk) * Bp(kk, jj);
            Bh(jj, ii) += Cb(kk, jj) * Ap(ii, kk);
          }
        }
      }
    } else if (AT) {
      Mat2x2MatTransMultAddCore(B, Ch, Ah);
      Mat2x2MatMultAddCore(A, Ch, Bh);

      for (int ii = 0; ii < 2; ii++) {
        for (int jj = 0; jj < 2; jj++) {
          for (int kk = 0; kk < 2; kk++) {
            Ah(jj, ii) += Cb(ii, kk) * Bp(jj, kk);
            Bh(ii, jj) += Cb(kk, jj) * Ap(ii, kk);
          }
        }
      }
    } else if (BT) {
      Mat2x2MatMultAddCore(Ch, B, Ah);
      MatTrans2x2MatMultAddCore(Ch, A, Bh);

      for (int ii = 0; ii < 2; ii++) {
        for (int jj = 0; jj < 2; jj++) {
          for (int kk = 0; kk < 2; kk++) {
            Ah(ii, jj) += Cb(ii, kk) * Bp(kk, jj);
            Bh(jj, ii) += Cb(kk, jj) * Ap(kk, ii);
          }
        }
      }
    } else {
      Mat2x2MatTransMultAddCore(Ch, B, Ah);
      MatTrans2x2MatMultAddCore(A, Ch, Bh);

      for (int ii = 0; ii < 2; ii++) {
        for (int jj = 0; jj < 2; jj++) {
          for (int kk = 0; kk < 2; kk++) {
            Ah(ii, jj) += Cb(ii, kk) * Bp(jj, kk);
            Bh(ii, jj) += Cb(kk, jj) * Ap(kk, ii);
          }
        }
      }
    }
  }

  A2DMat<Mat<ScalarType, 2, 2>>& AObj;
  A2DMat<Mat<ScalarType, 2, 2>>& BObj;
  A2DMat<Mat<ScalarType, 2, 2>>& CObj;
};

template <typename ScalarType, bool AT = false, bool BT = false>
KOKKOS_FUNCTION A2DMat2x2MatMultExpr<ScalarType, AT, BT> MatMatMult(
    A2DMat<Mat<ScalarType, 2, 2>>& AObj, A2DMat<Mat<ScalarType, 2, 2>>& BObj,
    A2DMat<Mat<ScalarType, 2, 2>>& CObj) {
  return A2DMat2x2MatMultExpr<ScalarType, AT, BT>(AObj, BObj, CObj);
}

template <typename ScalarType, bool AT = false, bool BT = false>
class ADpMat2x2MatMultExpr
    : public ADExpression<ADpMat2x2MatMultExpr<ScalarType, AT, BT>> {
 public:
  typedef Mat<ScalarType, 2, 2> Mat2x2;

  KOKKOS_FUNCTION ADpMat2x2MatMultExpr(Mat2x2& A,
                                       ADMat<Mat<ScalarType, 2, 2>>& BObj,
                                       ADMat<Mat<ScalarType, 2, 2>>& CObj)
      : A(A), BObj(BObj), CObj(CObj) {
    const Mat<ScalarType, 2, 2>& B = BObj.value();
    Mat<ScalarType, 2, 2>& C = CObj.value();

    if (AT && BT) {
      MatTrans2x2MatTransMultCore(A, B, C);
    } else if (AT) {
      MatTrans2x2MatMultCore(A, B, C);
    } else if (BT) {
      Mat2x2MatTransMultCore(A, B, C);
    } else {
      Mat2x2MatMultCore(A, B, C);
    }
  }

  KOKKOS_FUNCTION void forward() {
    const Mat<ScalarType, 2, 2>& Bb = BObj.bvalue();
    Mat<ScalarType, 2, 2>& Cb = CObj.bvalue();

    if (AT && BT) {
      MatTrans2x2MatTransMultCore(A, Bb, Cb);
    } else if (AT) {
      MatTrans2x2MatMultCore(A, Bb, Cb);
    } else if (BT) {
      Mat2x2MatTransMultCore(A, Bb, Cb);
    } else {
      Mat2x2MatMultCore(A, Bb, Cb);
    }
  }

  KOKKOS_FUNCTION void reverse() {
    Mat<ScalarType, 2, 2>& Bb = BObj.bvalue();
    const Mat<ScalarType, 2, 2>& Cb = CObj.bvalue();

    if (AT && BT) {
      MatTrans2x2MatTransMultAddCore(Cb, A, Bb);
    } else if (AT) {
      Mat2x2MatMultAddCore(A, Cb, Bb);
    } else if (BT) {
      MatTrans2x2MatMultAddCore(Cb, A, Bb);
    } else {
      MatTrans2x2MatMultAddCore(A, Cb, Bb);
    }
  }

  Mat2x2& A;
  ADMat<Mat<ScalarType, 2, 2>>& BObj;
  ADMat<Mat<ScalarType, 2, 2>>& CObj;
};

template <class ScalarType, bool AT = false, bool BT = false>
KOKKOS_FUNCTION ADpMat2x2MatMultExpr<ScalarType, AT, BT> MatMatMult(
    Mat<ScalarType, 2, 2>& A, ADMat<Mat<ScalarType, 2, 2>>& BObj,
    ADMat<Mat<ScalarType, 2, 2>>& CObj) {
  return ADpMat2x2MatMultExpr<ScalarType, AT, BT>(A, BObj, CObj);
}

template <typename ScalarType, bool AT = false, bool BT = false>
class ADMat2x2pMatMultExpr
    : public ADExpression<ADMat2x2pMatMultExpr<ScalarType, AT, BT>> {
 public:
  typedef Mat<ScalarType, 2, 2> Mat2x2;

  KOKKOS_FUNCTION ADMat2x2pMatMultExpr(ADMat<Mat<ScalarType, 2, 2>>& AObj,
                                       Mat2x2& B,
                                       ADMat<Mat<ScalarType, 2, 2>>& CObj)
      : AObj(AObj), B(B), CObj(CObj) {
    const Mat<ScalarType, 2, 2>& A = AObj.value();
    Mat<ScalarType, 2, 2>& C = CObj.value();

    if (AT && BT) {
      MatTrans2x2MatTransMultCore(A, B, C);
    } else if (AT) {
      MatTrans2x2MatMultCore(A, B, C);
    } else if (BT) {
      Mat2x2MatTransMultCore(A, B, C);
    } else {
      Mat2x2MatMultCore(A, B, C);
    }
  }

  KOKKOS_FUNCTION void forward() {
    const Mat<ScalarType, 2, 2>& Ab = AObj.bvalue();
    Mat<ScalarType, 2, 2>& Cb = CObj.bvalue();

    if (AT && BT) {
      MatTrans2x2MatTransMultCore(Ab, B, Cb);
    } else if (AT) {
      MatTrans2x2MatMultCore(Ab, B, Cb);
    } else if (BT) {
      Mat2x2MatTransMultCore(Ab, B, Cb);
    } else {
      Mat2x2MatMultCore(Ab, B, Cb);
    }
  }

  KOKKOS_FUNCTION void reverse() {
    Mat<ScalarType, 2, 2>& Ab = AObj.bvalue();
    const Mat<ScalarType, 2, 2>& Cb = CObj.bvalue();

    if (AT && BT) {
      MatTrans2x2MatTransMultAddCore(B, Cb, Ab);
    } else if (AT) {
      Mat2x2MatTransMultAddCore(B, Cb, Ab);
    } else if (BT) {
      Mat2x2MatMultAddCore(Cb, B, Ab);
    } else {
      Mat2x2MatTransMultAddCore(Cb, B, Ab);
    }
  }

  ADMat<Mat<ScalarType, 2, 2>>& AObj;
  Mat2x2& B;
  ADMat<Mat<ScalarType, 2, 2>>& CObj;
};

template <typename ScalarType, bool AT = false, bool BT = false>
KOKKOS_FUNCTION ADMat2x2pMatMultExpr<ScalarType, AT, BT> MatMatMult(
    ADMat<Mat<ScalarType, 2, 2>>& AObj, Mat<ScalarType, 2, 2>& B,
    ADMat<Mat<ScalarType, 2, 2>>& CObj) {
  return ADMat2x2pMatMultExpr<ScalarType, AT, BT>(AObj, B, CObj);
}

template <typename ScalarType, bool AT = false, bool BT = false>
class A2DpMat2x2MatMultExpr
    : public A2DExpression<A2DpMat2x2MatMultExpr<ScalarType, AT, BT>> {
 public:
  typedef Mat<ScalarType, 2, 2> Mat2x2;

  KOKKOS_FUNCTION A2DpMat2x2MatMultExpr(Mat2x2& A,
                                        A2DMat<Mat<ScalarType, 2, 2>>& BObj,
                                        A2DMat<Mat<ScalarType, 2, 2>>& CObj)
      : A(A), BObj(BObj), CObj(CObj) {
    const Mat<ScalarType, 2, 2>& B = BObj.value();
    Mat<ScalarType, 2, 2>& C = CObj.value();

    if (AT && BT) {
      MatTrans2x2MatTransMultCore(A, B, C);
    } else if (AT) {
      MatTrans2x2MatMultCore(A, B, C);
    } else if (BT) {
      Mat2x2MatTransMultCore(A, B, C);
    } else {
      Mat2x2MatMultCore(A, B, C);
    }
  }

  KOKKOS_FUNCTION void reverse() {
    Mat<ScalarType, 2, 2>& Bb = BObj.bvalue();
    const Mat<ScalarType, 2, 2>& Cb = CObj.bvalue();

    if (AT && BT) {
      MatTrans2x2MatTransMultAddCore(Cb, A, Bb);
    } else if (AT) {
      Mat2x2MatMultAddCore(A, Cb, Bb);
    } else if (BT) {
      MatTrans2x2MatMultAddCore(Cb, A, Bb);
    } else {
      MatTrans2x2MatMultAddCore(A, Cb, Bb);
    }
  }

  KOKKOS_FUNCTION void hforward() {
    const Mat<ScalarType, 2, 2>& Bp = BObj.pvalue();
    Mat<ScalarType, 2, 2>& Cp = CObj.pvalue();

    if (AT && BT) {
      MatTrans2x2MatTransMultCore(A, Bp, Cp);
    } else if (AT) {
      MatTrans2x2MatMultCore(A, Bp, Cp);
    } else if (BT) {
      Mat2x2MatTransMultCore(A, Bp, Cp);
    } else {
      Mat2x2MatMultCore(A, Bp, Cp);
    }
  }

  KOKKOS_FUNCTION void hreverse() {
    Mat<ScalarType, 2, 2>& Bh = BObj.hvalue();
    const Mat<ScalarType, 2, 2>& Ch = CObj.hvalue();

    if (AT && BT) {
      MatTrans2x2MatTransMultAddCore(Ch, A, Bh);
    } else if (AT) {
      Mat2x2MatMultAddCore(A, Ch, Bh);
    } else if (BT) {
      MatTrans2x2MatMultAddCore(Ch, A, Bh);
    } else {
      MatTrans2x2MatMultAddCore(A, Ch, Bh);
    }
  }

  Mat2x2& A;
  A2DMat<Mat<ScalarType, 2, 2>>& BObj;
  A2DMat<Mat<ScalarType, 2, 2>>& CObj;
};

template <class ScalarType, bool AT = false, bool BT = false>
KOKKOS_FUNCTION A2DpMat2x2MatMultExpr<ScalarType, AT, BT> MatMatMult(
    Mat<ScalarType, 2, 2>& A, A2DMat<Mat<ScalarType, 2, 2>>& BObj,
    A2DMat<Mat<ScalarType, 2, 2>>& CObj) {
  return A2DpMat2x2MatMultExpr<ScalarType, AT, BT>(A, BObj, CObj);
}

template <typename ScalarType, bool AT = false, bool BT = false>
class A2DMat2x2pMatMultExpr
    : public ADExpression<A2DMat2x2pMatMultExpr<ScalarType, AT, BT>> {
 public:
  typedef Mat<ScalarType, 2, 2> Mat2x2;

  KOKKOS_FUNCTION A2DMat2x2pMatMultExpr(A2DMat<Mat<ScalarType, 2, 2>>& AObj,
                                        Mat2x2& B,
                                        A2DMat<Mat<ScalarType, 2, 2>>& CObj)
      : AObj(AObj), B(B), CObj(CObj) {
    const Mat<ScalarType, 2, 2>& A = AObj.value();
    Mat<ScalarType, 2, 2>& C = CObj.value();

    if (AT && BT) {
      MatTrans2x2MatTransMultCore(A, B, C);
    } else if (AT) {
      MatTrans2x2MatMultCore(A, B, C);
    } else if (BT) {
      Mat2x2MatTransMultCore(A, B, C);
    } else {
      Mat2x2MatMultCore(A, B, C);
    }
  }

  KOKKOS_FUNCTION void reverse() {
    Mat<ScalarType, 2, 2>& Ab = AObj.bvalue();
    const Mat<ScalarType, 2, 2>& Cb = CObj.bvalue();

    if (AT && BT) {
      MatTrans2x2MatTransMultAddCore(B, Cb, Ab);
    } else if (AT) {
      Mat2x2MatTransMultAddCore(B, Cb, Ab);
    } else if (BT) {
      Mat2x2MatMultAddCore(Cb, B, Ab);
    } else {
      Mat2x2MatTransMultAddCore(Cb, B, Ab);
    }
  }

  KOKKOS_FUNCTION void hforward() {
    const Mat<ScalarType, 2, 2>& Ap = AObj.pvalue();
    Mat<ScalarType, 2, 2>& Cp = CObj.pvalue();

    if (AT && BT) {
      MatTrans2x2MatTransMultCore(Ap, B, Cp);
    } else if (AT) {
      MatTrans2x2MatMultCore(Ap, B, Cp);
    } else if (BT) {
      Mat2x2MatTransMultCore(Ap, B, Cp);
    } else {
      Mat2x2MatMultCore(Ap, B, Cp);
    }
  }

  KOKKOS_FUNCTION void hreverse() {
    Mat<ScalarType, 2, 2>& Ah = AObj.hvalue();
    const Mat<ScalarType, 2, 2>& Ch = CObj.hvalue();

    if (AT && BT) {
      MatTrans2x2MatTransMultAddCore(B, Ch, Ah);
    } else if (AT) {
      Mat2x2MatTransMultAddCore(B, Ch, Ah);
    } else if (BT) {
      Mat2x2MatMultAddCore(Ch, B, Ah);
    } else {
      Mat2x2MatTransMultAddCore(Ch, B, Ah);
    }
  }

  A2DMat<Mat<ScalarType, 2, 2>>& AObj;
  Mat2x2& B;
  A2DMat<Mat<ScalarType, 2, 2>>& CObj;
};

template <typename ScalarType, bool AT = false, bool BT = false>
KOKKOS_FUNCTION A2DMat2x2pMatMultExpr<ScalarType, AT, BT> MatMatMult(
    A2DMat<Mat<ScalarType, 2, 2>>& AObj, Mat<ScalarType, 2, 2>& B,
    A2DMat<Mat<ScalarType, 2, 2>>& CObj) {
  return A2DMat2x2pMatMultExpr<ScalarType, AT, BT>(AObj, B, CObj);
}

// Mat2x2Det
template <typename ScalarType>
KOKKOS_FUNCTION void MatDet(const Mat<ScalarType, 2, 2>& A, ScalarType& det) {
  det = A(1, 1) * A(0, 0) - A(0, 1) * A(1, 0);
}

template <class ScalarType>
class ADMat2x2DetExpr : public ADExpression<ADMat2x2DetExpr<ScalarType>> {
 public:
  KOKKOS_FUNCTION ADMat2x2DetExpr(ADMat<Mat<ScalarType, 2, 2>>& AObj,
                                  ADScalar<ScalarType>& detObj)
      : AObj(AObj), detObj(detObj) {
    const Mat<ScalarType, 2, 2>& A = AObj.value();

    detObj.value = (A(0, 0) * A(1, 1) - A(0, 1) * A(1, 0));
  }

  KOKKOS_FUNCTION void forward() {
    const Mat<ScalarType, 2, 2>& A = AObj.value();
    const Mat<ScalarType, 2, 2>& Ad = AObj.bvalue();

    detObj.bvalue = (Ad(0, 0) * A(1, 1) + Ad(1, 1) * A(0, 0) -
                     Ad(0, 1) * A(1, 0) - Ad(1, 0) * A(0, 1));
  }

  KOKKOS_FUNCTION void reverse() {
    const ScalarType& bdet = detObj.bvalue;
    const Mat<ScalarType, 2, 2>& A = AObj.value();
    Mat<ScalarType, 2, 2>& Ab = AObj.bvalue();

    Ab(0, 0) += A(1, 1) * bdet;
    Ab(1, 1) += A(0, 0) * bdet;
    Ab(0, 1) -= A(1, 0) * bdet;
    Ab(1, 0) -= A(0, 1) * bdet;
  }

  ADMat<Mat<ScalarType, 2, 2>>& AObj;
  ADScalar<ScalarType>& detObj;
};

template <typename ScalarType>
KOKKOS_FUNCTION ADMat2x2DetExpr<ScalarType> MatDet(
    ADMat<Mat<ScalarType, 2, 2>>& A, ADScalar<ScalarType>& det) {
  return ADMat2x2DetExpr<ScalarType>(A, det);
}

template <class ScalarType>
class A2DMat2x2DetExpr : public ADExpression<A2DMat2x2DetExpr<ScalarType>> {
 public:
  KOKKOS_FUNCTION A2DMat2x2DetExpr(A2DMat<Mat<ScalarType, 2, 2>>& AObj,
                                   A2DScalar<ScalarType>& detObj)
      : AObj(AObj), detObj(detObj) {
    const Mat<ScalarType, 2, 2>& A = AObj.value();

    detObj.value = (A(0, 0) * A(1, 1) - A(0, 1) * A(1, 0));
  }

  KOKKOS_FUNCTION void reverse() {
    const ScalarType& bdet = detObj.bvalue;
    const Mat<ScalarType, 2, 2>& A = AObj.value();
    Mat<ScalarType, 2, 2>& Ab = AObj.bvalue();

    Ab(0, 0) += A(1, 1) * bdet;
    Ab(1, 1) += A(0, 0) * bdet;
    Ab(0, 1) -= A(1, 0) * bdet;
    Ab(1, 0) -= A(0, 1) * bdet;
  }

  KOKKOS_FUNCTION void hforward() {
    const Mat<ScalarType, 2, 2>& A = AObj.value();
    const Mat<ScalarType, 2, 2>& Ap = AObj.pvalue();
    detObj.pvalue = (Ap(0, 0) * A(1, 1) + Ap(1, 1) * A(0, 0) -
                     Ap(0, 1) * A(1, 0) - Ap(1, 0) * A(0, 1));
  }

  KOKKOS_FUNCTION void hreverse() {
    const ScalarType bdet = detObj.bvalue;
    const Mat<ScalarType, 2, 2>& A = AObj.value();
    const ScalarType hdet = detObj.hvalue;
    const Mat<ScalarType, 2, 2>& Ap = AObj.pvalue();
    Mat<ScalarType, 2, 2>& Ah = AObj.hvalue();

    Ah(0, 0) += Ap(1, 1) * bdet;
    Ah(0, 1) -= Ap(1, 0) * bdet;
    Ah(1, 0) -= Ap(0, 1) * bdet;
    Ah(1, 1) += Ap(0, 0) * bdet;

    Ah(0, 0) += A(1, 1) * hdet;
    Ah(1, 1) += A(0, 0) * hdet;
    Ah(0, 1) -= A(1, 0) * hdet;
    Ah(1, 0) -= A(0, 1) * hdet;
  }

  A2DMat<Mat<ScalarType, 2, 2>>& AObj;
  A2DScalar<ScalarType>& detObj;
};

template <typename ScalarType>
KOKKOS_FUNCTION A2DMat2x2DetExpr<ScalarType> MatDet(
    A2DMat<Mat<ScalarType, 2, 2>>& A, A2DScalar<ScalarType>& det) {
  return A2DMat2x2DetExpr<ScalarType>(A, det);
}

// Mat2x2Inverse
template <typename ScalarType>
KOKKOS_FUNCTION void MatInverse(const Mat<ScalarType, 2, 2>& A,
                                Mat<ScalarType, 2, 2>& Ainv) {
  ScalarType det = A(0, 0) * A(1, 1) - A(1, 0) * A(0, 1);
  ScalarType detinv = 1.0 / det;

  Ainv(0, 0) = A(1, 1) * detinv;
  Ainv(0, 1) = -A(0, 1) * detinv;
  Ainv(1, 0) = -A(1, 0) * detinv;
  Ainv(1, 1) = A(0, 0) * detinv;
}

template <typename ScalarType>
class ADMat2x2InverseExpr
    : public ADExpression<ADMat2x2InverseExpr<ScalarType>> {
 public:
  KOKKOS_FUNCTION ADMat2x2InverseExpr(ADMat<Mat<ScalarType, 2, 2>>& AObj,
                                      ADMat<Mat<ScalarType, 2, 2>>& AinvObj)
      : AObj(AObj), AinvObj(AinvObj) {
    const Mat<ScalarType, 2, 2>& A = AObj.value();
    Mat<ScalarType, 2, 2>& Ainv = AinvObj.value();

    ScalarType det = A(0, 0) * A(1, 1) - A(1, 0) * A(0, 1);
    ScalarType detinv = 1.0 / det;

    Ainv(0, 0) = A(1, 1) * detinv;
    Ainv(0, 1) = -A(0, 1) * detinv;
    Ainv(1, 0) = -A(1, 0) * detinv;
    Ainv(1, 1) = A(0, 0) * detinv;
  }

  KOKKOS_FUNCTION void forward() {
    const Mat<ScalarType, 2, 2>& Ainv = AinvObj.value();
    const Mat<ScalarType, 2, 2>& Ad = AObj.bvalue();
    Mat<ScalarType, 2, 2>& Ainvd = AinvObj.bvalue();

    Mat<ScalarType, 2, 2> tmp;
    Mat2x2MatMultCore(Ainv, Ad, tmp);
    Mat2x2MatMultScaleCore(ScalarType(-1.0), tmp, Ainv, Ainvd);
  }

  KOKKOS_FUNCTION void reverse() {
    const Mat<ScalarType, 2, 2>& Ainv = AinvObj.value();
    const Mat<ScalarType, 2, 2>& Ainvb = AinvObj.bvalue();
    Mat<ScalarType, 2, 2>& Ab = AObj.bvalue();

    Mat<ScalarType, 2, 2> tmp;
    MatTrans2x2MatMultCore(Ainv, Ainvb, tmp);
    Mat2x2MatTransMultAddScaleCore(ScalarType(-1.0), tmp, Ainv, Ab);
  }

  ADMat<Mat<ScalarType, 2, 2>>& AObj;
  ADMat<Mat<ScalarType, 2, 2>>& AinvObj;
};

template <typename ScalarType>
KOKKOS_FUNCTION ADMat2x2InverseExpr<ScalarType> MatInverse(
    ADMat<Mat<ScalarType, 2, 2>>& AObj, ADMat<Mat<ScalarType, 2, 2>>& AinvObj) {
  return ADMat2x2InverseExpr<ScalarType>(AObj, AinvObj);
}

template <typename ScalarType>
class A2DMat2x2InverseExpr
    : public A2DExpression<A2DMat2x2InverseExpr<ScalarType>> {
 public:
  KOKKOS_FUNCTION A2DMat2x2InverseExpr(A2DMat<Mat<ScalarType, 2, 2>>& AObj,
                                       A2DMat<Mat<ScalarType, 2, 2>>& AinvObj)
      : AObj(AObj), AinvObj(AinvObj) {
    const Mat<ScalarType, 2, 2>& A = AObj.value();
    Mat<ScalarType, 2, 2>& Ainv = AinvObj.value();

    ScalarType det = A(0, 0) * A(1, 1) - A(1, 0) * A(0, 1);
    ScalarType detinv = 1.0 / det;

    Ainv(0, 0) = A(1, 1) * detinv;
    Ainv(0, 1) = -A(0, 1) * detinv;
    Ainv(1, 0) = -A(1, 0) * detinv;
    Ainv(1, 1) = A(0, 0) * detinv;
  }

  KOKKOS_FUNCTION void reverse() {
    const Mat<ScalarType, 2, 2>& Ainv = AinvObj.value();
    const Mat<ScalarType, 2, 2>& Ainvd = AinvObj.bvalue();
    Mat<ScalarType, 2, 2>& Ad = AObj.bvalue();

    Mat<ScalarType, 2, 2> tmp;
    MatTrans2x2MatMultCore(Ainv, Ainvd, tmp);
    Mat2x2MatTransMultAddScaleCore(ScalarType(-1.0), tmp, Ainv, Ad);
  }

  KOKKOS_FUNCTION void hforward() {
    Mat<ScalarType, 2, 2> tmp;
    const Mat<ScalarType, 2, 2>& Ainv = AinvObj.value();

    const Mat<ScalarType, 2, 2>& Ap = AObj.pvalue();
    Mat<ScalarType, 2, 2>& Ainvp = AinvObj.pvalue();

    Mat2x2MatMultCore(Ainv, Ap, tmp);
    Mat2x2MatMultScaleCore(ScalarType(-1.0), tmp, Ainv, Ainvp);
  }

  // hA = A^{-T} * Ap^{T} * A^{-T} * Ainvb * A^{-T} +
  //      A^{-T} * Ainvb * A^{-T} * Ap^{T} * A^{-T} =
  //    = - (A^{-T} * Ap^{T} * Ab + Ab * Ap^{T} * A^{-T})

  KOKKOS_FUNCTION void hreverse() {
    // Temporary matrix
    Mat<ScalarType, 2, 2> tmp, tmp2;

    const Mat<ScalarType, 2, 2>& Ainv = AinvObj.value();
    const Mat<ScalarType, 2, 2>& Ab = AObj.bvalue();

    const Mat<ScalarType, 2, 2>& Ap = AObj.pvalue();
    const Mat<ScalarType, 2, 2>& Ainvh = AinvObj.hvalue();
    Mat<ScalarType, 2, 2>& Ah = AObj.hvalue();

    // Ainv^{T} * Ap^{T} * Ab
    MatTrans2x2MatTransMultCore(Ainv, Ap, tmp);
    Mat2x2MatMultAddScaleCore(ScalarType(-1.0), tmp, Ab, Ah);

    // Ab * Ap^{T} * A^{-T}
    Mat2x2MatTransMultCore(Ab, Ap, tmp);
    Mat2x2MatTransMultAddScaleCore(ScalarType(-1.0), tmp, Ainv, Ah);

    MatTrans2x2MatMultCore(Ainv, Ainvh, tmp);
    Mat2x2MatTransMultAddScaleCore(ScalarType(-1.0), tmp, Ainv, Ah);
  }

  A2DMat<Mat<ScalarType, 2, 2>>& AObj;
  A2DMat<Mat<ScalarType, 2, 2>>& AinvObj;
};

template <typename ScalarType>
KOKKOS_FUNCTION A2DMat2x2InverseExpr<ScalarType> MatInverse(
    A2DMat<Mat<ScalarType, 2, 2>>& AObj,
    A2DMat<Mat<ScalarType, 2, 2>>& AinvObj) {
  return A2DMat2x2InverseExpr<ScalarType>(AObj, AinvObj);
}

// SymmTrace
template <typename ScalarType>
KOKKOS_FUNCTION void SymmTrace(const SymMat<ScalarType, 2>& S,
                               ScalarType& trace) {
  trace = S(0, 0) + S(1, 1);
}

template <class ScalarType>
class ADSymm2x2TraceExpr : public ADExpression<ADSymm2x2TraceExpr<ScalarType>> {
 public:
  KOKKOS_FUNCTION ADSymm2x2TraceExpr(ADMat<SymMat<ScalarType, 2>>& SObj,
                                     ADScalar<ScalarType>& output)
      : SObj(SObj), output(output) {
    const SymMat<ScalarType, 2>& S = SObj.value();
    output.value = S(0, 0) + S(1, 1);
  }

  KOKKOS_FUNCTION void forward() {
    const SymMat<ScalarType, 2>& Sd = SObj.bvalue();
    output.bvalue = Sd(0, 0) + Sd(1, 1);
  }

  KOKKOS_FUNCTION void reverse() {
    SymMat<ScalarType, 2>& Sb = SObj.bvalue();

    Sb(0, 0) += output.bvalue;
    Sb(1, 1) += output.bvalue;
  }

  ADMat<SymMat<ScalarType, 2>>& SObj;
  ADScalar<ScalarType>& output;
};

template <class ScalarType>
KOKKOS_FUNCTION ADSymm2x2TraceExpr<ScalarType> SymmTrace(
    ADMat<SymMat<ScalarType, 2>>& S, ADScalar<ScalarType>& trace) {
  return ADSymm2x2TraceExpr<ScalarType>(S, trace);
}

template <class ScalarType>
class A2DSymm2x2TraceExpr
    : public A2DExpression<A2DSymm2x2TraceExpr<ScalarType>> {
 public:
  KOKKOS_FUNCTION A2DSymm2x2TraceExpr(A2DMat<SymMat<ScalarType, 2>>& SObj,
                                      A2DScalar<ScalarType>& output)
      : SObj(SObj), output(output) {
    const SymMat<ScalarType, 2>& S = SObj.value();
    output.value = S(0, 0) + S(1, 1);
  }

  KOKKOS_FUNCTION void reverse() {
    SymMat<ScalarType, 2>& Sb = SObj.bvalue();

    Sb(0, 0) += output.bvalue;
    Sb(1, 1) += output.bvalue;
  }

  // Compute E.pvalue() = J * Ux.pvalue()
  KOKKOS_FUNCTION void hforward() {
    const SymMat<ScalarType, 2>& Sp = SObj.pvalue();
    output.pvalue = Sp(0, 0) + Sp(1, 1);
  }

  KOKKOS_FUNCTION void hreverse() {
    SymMat<ScalarType, 2>& Sh = SObj.hvalue();

    Sh(0, 0) += output.hvalue;
    Sh(1, 1) += output.hvalue;
  }

  A2DMat<SymMat<ScalarType, 2>>& SObj;
  A2DScalar<ScalarType>& output;
};

template <class ScalarType>
KOKKOS_FUNCTION A2DSymm2x2TraceExpr<ScalarType> SymmTrace(
    A2DMat<SymMat<ScalarType, 2>>& S, A2DScalar<ScalarType>& trace) {
  return A2DSymm2x2TraceExpr<ScalarType>(S, trace);
}

// Symm2x2SymmMultTrace
template <typename ScalarType>
KOKKOS_FUNCTION void SymmSymmMultTrace(const SymMat<ScalarType, 2>& S,
                                       const SymMat<ScalarType, 2>& E,
                                       ScalarType& trace) {
  trace = S(0, 0) * E(0, 0) + S(1, 1) * E(1, 1) + 2.0 * (S(0, 1) * E(0, 1));
}

template <class ScalarType>
class ADSymm2x2SymmMultTraceExpr
    : public ADExpression<ADSymm2x2SymmMultTraceExpr<ScalarType>> {
 public:
  KOKKOS_FUNCTION ADSymm2x2SymmMultTraceExpr(ADMat<SymMat<ScalarType, 2>>& SObj,
                                             ADMat<SymMat<ScalarType, 2>>& EObj,
                                             ADScalar<ScalarType>& output)
      : SObj(SObj), EObj(EObj), output(output) {
    const SymMat<ScalarType, 2>& S = SObj.value();
    const SymMat<ScalarType, 2>& E = EObj.value();

    output.value =
        S(0, 0) * E(0, 0) + S(1, 1) * E(1, 1) + 2.0 * (S(0, 1) * E(0, 1));
  }

  KOKKOS_FUNCTION void forward() {
    const SymMat<ScalarType, 2>& E = EObj.value();
    const SymMat<ScalarType, 2>& Ed = EObj.bvalue();
    const SymMat<ScalarType, 2>& S = SObj.value();
    const SymMat<ScalarType, 2>& Sd = SObj.bvalue();

    output.bvalue = S(0, 0) * Ed(0, 0) + S(1, 1) * Ed(1, 1) +
                    2.0 * (S(0, 1) * Ed(0, 1)) + Sd(0, 0) * E(0, 0) +
                    Sd(1, 1) * E(1, 1) + 2.0 * (Sd(0, 1) * E(0, 1));
  }

  KOKKOS_FUNCTION void reverse() {
    const SymMat<ScalarType, 2>& E = EObj.value();
    SymMat<ScalarType, 2>& Eb = EObj.bvalue();
    const SymMat<ScalarType, 2>& S = SObj.value();
    SymMat<ScalarType, 2>& Sb = SObj.bvalue();

    Eb(0, 0) += output.bvalue * S(0, 0);
    Eb(1, 1) += output.bvalue * S(1, 1);
    Eb(0, 1) += 2.0 * output.bvalue * S(0, 1);

    Sb(0, 0) += output.bvalue * E(0, 0);
    Sb(1, 1) += output.bvalue * E(1, 1);
    Sb(0, 1) += 2.0 * output.bvalue * E(0, 1);
  }

  ADMat<SymMat<ScalarType, 2>>& SObj;
  ADMat<SymMat<ScalarType, 2>>& EObj;
  ADScalar<ScalarType>& output;
};

template <class ScalarType>
KOKKOS_FUNCTION ADSymm2x2SymmMultTraceExpr<ScalarType> SymmSymmMultTrace(
    ADMat<SymMat<ScalarType, 2>>& S, ADMat<SymMat<ScalarType, 2>>& E,
    ADScalar<ScalarType>& trace) {
  return ADSymm2x2SymmMultTraceExpr<ScalarType>(S, E, trace);
}

template <class ScalarType>
class A2DSymm2x2SymmMultTraceExpr
    : public A2DExpression<A2DSymm2x2SymmMultTraceExpr<ScalarType>> {
 public:
  KOKKOS_FUNCTION A2DSymm2x2SymmMultTraceExpr(
      A2DMat<SymMat<ScalarType, 2>>& SObj, A2DMat<SymMat<ScalarType, 2>>& EObj,
      A2DScalar<ScalarType>& output)
      : SObj(SObj), EObj(EObj), output(output) {
    const SymMat<ScalarType, 2>& S = SObj.value();
    const SymMat<ScalarType, 2>& E = EObj.value();

    output.value =
        S(0, 0) * E(0, 0) + S(1, 1) * E(1, 1) + 2.0 * (S(0, 1) * E(0, 1));
  }

  KOKKOS_FUNCTION void reverse() {
    const SymMat<ScalarType, 2>& E = EObj.value();
    SymMat<ScalarType, 2>& Eb = EObj.bvalue();
    const SymMat<ScalarType, 2>& S = SObj.value();
    SymMat<ScalarType, 2>& Sb = SObj.bvalue();

    Eb(0, 0) += output.bvalue * S(0, 0);
    Eb(1, 1) += output.bvalue * S(1, 1);
    Eb(0, 1) += 2.0 * output.bvalue * S(0, 1);

    Sb(0, 0) += output.bvalue * E(0, 0);
    Sb(1, 1) += output.bvalue * E(1, 1);
    Sb(0, 1) += 2.0 * output.bvalue * E(0, 1);
  }

  // Compute E.pvalue() = J * Ux.pvalue()
  KOKKOS_FUNCTION void hforward() {
    const SymMat<ScalarType, 2>& E = EObj.value();
    const SymMat<ScalarType, 2>& S = SObj.value();

    const SymMat<ScalarType, 2>& Ep = EObj.pvalue();
    const SymMat<ScalarType, 2>& Sp = SObj.pvalue();

    output.pvalue = S(0, 0) * Ep(0, 0) + S(1, 1) * Ep(1, 1) +
                    2.0 * S(0, 1) * Ep(0, 1) + Sp(0, 0) * E(0, 0) +
                    Sp(1, 1) * E(1, 1) + 2.0 * Sp(0, 1) * E(0, 1);
  }

  KOKKOS_FUNCTION void hreverse() {
    const SymMat<ScalarType, 2>& E = EObj.value();
    const SymMat<ScalarType, 2>& S = SObj.value();

    const SymMat<ScalarType, 2>& Ep = EObj.pvalue();
    const SymMat<ScalarType, 2>& Sp = SObj.pvalue();
    SymMat<ScalarType, 2>& Eh = EObj.hvalue();
    SymMat<ScalarType, 2>& Sh = SObj.hvalue();

    Eh(0, 0) += output.bvalue * Sp(0, 0);
    Eh(1, 1) += output.bvalue * Sp(1, 1);
    Eh(0, 1) += 2.0 * output.bvalue * Sp(0, 1);

    Sh(0, 0) += output.bvalue * Ep(0, 0);
    Sh(1, 1) += output.bvalue * Ep(1, 1);
    Sh(0, 1) += 2.0 * output.bvalue * Ep(0, 1);

    Eh(0, 0) += output.hvalue * S(0, 0);
    Eh(1, 1) += output.hvalue * S(1, 1);
    Eh(0, 1) += 2.0 * output.hvalue * S(0, 1);

    Sh(0, 0) += output.hvalue * E(0, 0);
    Sh(1, 1) += output.hvalue * E(1, 1);
    Sh(0, 1) += 2.0 * output.hvalue * E(0, 1);
  }

  A2DMat<SymMat<ScalarType, 2>>& SObj;
  A2DMat<SymMat<ScalarType, 2>>& EObj;
  A2DScalar<ScalarType>& output;
};

template <class ScalarType>
KOKKOS_FUNCTION A2DSymm2x2SymmMultTraceExpr<ScalarType> SymmSymmMultTrace(
    A2DMat<SymMat<ScalarType, 2>>& S, A2DMat<SymMat<ScalarType, 2>>& E,
    A2DScalar<ScalarType>& trace) {
  return A2DSymm2x2SymmMultTraceExpr<ScalarType>(S, E, trace);
}

template <class ScalarType>
KOKKOS_FUNCTION void SymmIsotropicConstitutive(const ScalarType& mu,
                                               const ScalarType& lambda,
                                               const SymMat<ScalarType, 2>& E,
                                               SymMat<ScalarType, 2>& S) {
  ScalarType tr = lambda * (E(0, 0) + E(1, 1));
  ScalarType mu2 = 2.0 * mu;
  S(0, 0) = mu2 * E(0, 0) + tr;
  S(0, 1) = mu2 * E(0, 1);
  S(1, 1) = mu2 * E(1, 1) + tr;
}

template <class ScalarType>
class ADSymm2x2IsotropicConstitutiveExpr
    : public ADExpression<ADSymm2x2IsotropicConstitutiveExpr<ScalarType>> {
 public:
  KOKKOS_FUNCTION ADSymm2x2IsotropicConstitutiveExpr(
      const ScalarType& mu, const ScalarType& lambda,
      ADMat<SymMat<ScalarType, 2>>& EObj, ADMat<SymMat<ScalarType, 2>>& SObj)
      : mu(mu), lambda(lambda), EObj(EObj), SObj(SObj) {
    const SymMat<ScalarType, 2>& E = EObj.value();
    SymMat<ScalarType, 2>& S = SObj.value();
    ScalarType tr = lambda * (E(0, 0) + E(1, 1));
    ScalarType mu2 = 2.0 * mu;
    S(0, 0) = mu2 * E(0, 0) + tr;
    S(0, 1) = mu2 * E(0, 1);
    S(1, 1) = mu2 * E(1, 1) + tr;
  }

  KOKKOS_FUNCTION void forward() {
    const SymMat<ScalarType, 2>& Ed = EObj.bvalue();
    SymMat<ScalarType, 2>& Sd = SObj.bvalue();

    ScalarType tr = lambda * (Ed(0, 0) + Ed(1, 1));
    ScalarType mu2 = 2.0 * mu;
    Sd(0, 0) = mu2 * Ed(0, 0) + tr;
    Sd(0, 1) = mu2 * Ed(0, 1);
    Sd(1, 1) = mu2 * Ed(1, 1) + tr;
  }

  KOKKOS_FUNCTION void reverse() {
    const SymMat<ScalarType, 2>& Sb = SObj.bvalue();
    SymMat<ScalarType, 2>& Eb = EObj.bvalue();

    ScalarType tr = lambda * (Sb(0, 0) + Sb(1, 1));
    ScalarType mu2 = 2.0 * mu;
    Eb(0, 0) += mu2 * Sb(0, 0) + tr;
    Eb(0, 1) += mu2 * Sb(0, 1);
    Eb(1, 1) += mu2 * Sb(1, 1) + tr;
  }

  const ScalarType& mu;
  const ScalarType& lambda;
  ADMat<SymMat<ScalarType, 2>>& EObj;
  ADMat<SymMat<ScalarType, 2>>& SObj;
};

template <class ScalarType>
KOKKOS_FUNCTION ADSymm2x2IsotropicConstitutiveExpr<ScalarType>
SymmIsotropicConstitutive(const ScalarType& mu, const ScalarType& lambda,
                          ADMat<SymMat<ScalarType, 2>>& E,
                          ADMat<SymMat<ScalarType, 2>>& S) {
  return ADSymm2x2IsotropicConstitutiveExpr<ScalarType>(mu, lambda, E, S);
}

template <class ScalarType>
class A2DSymm2x2IsotropicConstitutiveExpr
    : public A2DExpression<A2DSymm2x2IsotropicConstitutiveExpr<ScalarType>> {
 public:
  KOKKOS_FUNCTION A2DSymm2x2IsotropicConstitutiveExpr(
      const ScalarType& mu, const ScalarType& lambda,
      A2DMat<SymMat<ScalarType, 2>>& EObj, A2DMat<SymMat<ScalarType, 2>>& SObj)
      : mu(mu), lambda(lambda), EObj(EObj), SObj(SObj) {
    const SymMat<ScalarType, 2>& E = EObj.value();
    SymMat<ScalarType, 2>& S = SObj.value();
    ScalarType tr = lambda * (E(0, 0) + E(1, 1));
    ScalarType mu2 = 2.0 * mu;
    S(0, 0) = mu2 * E(0, 0) + tr;
    S(0, 1) = mu2 * E(0, 1);
    S(1, 1) = mu2 * E(1, 1) + tr;
  }

  KOKKOS_FUNCTION void reverse() {
    const SymMat<ScalarType, 2>& Sb = SObj.bvalue();
    SymMat<ScalarType, 2>& Eb = EObj.bvalue();

    ScalarType tr = lambda * (Sb(0, 0) + Sb(1, 1));
    ScalarType mu2 = 2.0 * mu;
    Eb(0, 0) += mu2 * Sb(0, 0) + tr;
    Eb(0, 1) += mu2 * Sb(0, 1);
    Eb(1, 1) += mu2 * Sb(1, 1) + tr;
  }

  KOKKOS_FUNCTION void hforward() {
    const SymMat<ScalarType, 2>& Ep = EObj.pvalue();
    SymMat<ScalarType, 2>& Sp = SObj.pvalue();

    ScalarType tr = lambda * (Ep(0, 0) + Ep(1, 1));
    ScalarType mu2 = 2.0 * mu;
    Sp(0, 0) = mu2 * Ep(0, 0) + tr;
    Sp(0, 1) = mu2 * Ep(0, 1);
    Sp(1, 1) = mu2 * Ep(1, 1) + tr;
  }

  KOKKOS_FUNCTION void hreverse() {
    const SymMat<ScalarType, 2>& Sh = SObj.hvalue();
    SymMat<ScalarType, 2>& Eh = EObj.hvalue();

    ScalarType tr = lambda * (Sh(0, 0) + Sh(1, 1));
    ScalarType mu2 = 2.0 * mu;
    Eh(0, 0) += mu2 * Sh(0, 0) + tr;
    Eh(0, 1) += mu2 * Sh(0, 1);
    Eh(1, 1) += mu2 * Sh(1, 1) + tr;
  }

  const ScalarType& mu;
  const ScalarType& lambda;
  A2DMat<SymMat<ScalarType, 2>>& EObj;
  A2DMat<SymMat<ScalarType, 2>>& SObj;
};

template <class ScalarType>
KOKKOS_FUNCTION A2DSymm2x2IsotropicConstitutiveExpr<ScalarType>
SymmIsotropicConstitutive(const ScalarType& mu, const ScalarType& lambda,
                          A2DMat<SymMat<ScalarType, 2>>& E,
                          A2DMat<SymMat<ScalarType, 2>>& S) {
  return A2DSymm2x2IsotropicConstitutiveExpr<ScalarType>(mu, lambda, E, S);
}

template <class ScalarType>
KOKKOS_FUNCTION void SymmIsotropicEnergy(const ScalarType& mu,
                                         const ScalarType& lambda,
                                         const SymMat<ScalarType, 2>& E,
                                         ScalarType& output) {
  ScalarType tr = E(0, 0) + E(1, 1);
  ScalarType trE =
      E(0, 0) * E(0, 0) + E(1, 1) * E(1, 1) + 2.0 * E(0, 1) * E(0, 1);

  output = mu * trE + 0.5 * lambda * tr * tr;
}

template <class ScalarType>
class ADSymm2x2IsotropicEnergyExpr
    : public ADExpression<ADSymm2x2IsotropicEnergyExpr<ScalarType>> {
 public:
  KOKKOS_FUNCTION ADSymm2x2IsotropicEnergyExpr(
      const ScalarType& mu, const ScalarType& lambda,
      ADMat<SymMat<ScalarType, 2>>& EObj, ADScalar<ScalarType>& output)
      : mu(mu), lambda(lambda), EObj(EObj), output(output) {
    const SymMat<ScalarType, 2>& E = EObj.value();
    ScalarType tr = E(0, 0) + E(1, 1);
    ScalarType trE =
        E(0, 0) * E(0, 0) + E(1, 1) * E(1, 1) + 2.0 * E(0, 1) * E(0, 1);

    output.value = mu * trE + 0.5 * lambda * tr * tr;
  }

  KOKKOS_FUNCTION void forward() {
    const SymMat<ScalarType, 2>& E = EObj.value();
    const SymMat<ScalarType, 2>& Ed = EObj.bvalue();
    ScalarType tr = (E(0, 0) + E(1, 1));
    ScalarType trd = (Ed(0, 0) + Ed(1, 1));
    ScalarType trEd = 2.0 * (E(0, 0) * Ed(0, 0) + E(1, 1) * Ed(1, 1) +
                             2.0 * E(0, 1) * Ed(0, 1));

    output.bvalue = mu * trEd + lambda * tr * trd;
  }

  KOKKOS_FUNCTION void reverse() {
    const SymMat<ScalarType, 2>& E = EObj.value();
    SymMat<ScalarType, 2>& Eb = EObj.bvalue();
    const ScalarType mu2 = 2.0 * mu;

    ScalarType tr = (E(0, 0) + E(1, 1));
    Eb(0, 0) += (mu2 * E(0, 0) + lambda * tr) * output.bvalue;
    Eb(1, 1) += (mu2 * E(1, 1) + lambda * tr) * output.bvalue;

    Eb(0, 1) += 2.0 * mu2 * E(0, 1) * output.bvalue;
  }

  const ScalarType& mu;
  const ScalarType& lambda;
  ADMat<SymMat<ScalarType, 2>>& EObj;
  ADScalar<ScalarType>& output;
};

template <class ScalarType>
KOKKOS_FUNCTION ADSymm2x2IsotropicEnergyExpr<ScalarType> SymmIsotropicEnergy(
    const ScalarType& mu, const ScalarType& lambda,
    ADMat<SymMat<ScalarType, 2>>& E, ADScalar<ScalarType>& output) {
  return ADSymm2x2IsotropicEnergyExpr<ScalarType>(mu, lambda, E, output);
}

template <class ScalarType>
class ADSymm2x2ADIsotropicEnergyExpr
    : public ADExpression<ADSymm2x2ADIsotropicEnergyExpr<ScalarType>> {
 public:
  KOKKOS_FUNCTION ADSymm2x2ADIsotropicEnergyExpr(
      ADScalar<ScalarType>& mu, ADScalar<ScalarType>& lambda,
      ADMat<SymMat<ScalarType, 2>>& EObj, ADScalar<ScalarType>& output)
      : mu(mu), lambda(lambda), EObj(EObj), output(output) {
    const SymMat<ScalarType, 2>& E = EObj.value();
    ScalarType tr = E(0, 0) + E(1, 1);
    ScalarType trE =
        E(0, 0) * E(0, 0) + E(1, 1) * E(1, 1) + 2.0 * E(0, 1) * E(0, 1);

    output.value = mu.value * trE + 0.5 * lambda.value * tr * tr;
  }

  KOKKOS_FUNCTION void forward() {
    const SymMat<ScalarType, 2>& E = EObj.value();
    const SymMat<ScalarType, 2>& Ed = EObj.bvalue();
    ScalarType tr = E(0, 0) + E(1, 1);
    ScalarType trd = Ed(0, 0) + Ed(1, 1);
    ScalarType trE =
        E(0, 0) * E(0, 0) + E(1, 1) * E(1, 1) + 2.0 * E(0, 1) * E(0, 1);
    ScalarType trEd = 2.0 * (E(0, 0) * Ed(0, 0) + E(1, 1) * Ed(1, 1) +
                             2.0 * E(0, 1) * Ed(0, 1));

    output.bvalue = mu.value * trEd + lambda.value * tr * trd +
                    mu.bvalue * trE + 0.5 * lambda.bvalue * tr * tr;
  }

  KOKKOS_FUNCTION void reverse() {
    const SymMat<ScalarType, 2>& E = EObj.value();
    SymMat<ScalarType, 2>& Eb = EObj.bvalue();
    const ScalarType mu2 = 2.0 * mu.value;

    ScalarType tr = E(0, 0) + E(1, 1);
    ScalarType trE =
        E(0, 0) * E(0, 0) + E(1, 1) * E(1, 1) + 2.0 * E(0, 1) * E(0, 1);

    Eb(0, 0) += (mu2 * E(0, 0) + lambda.value * tr) * output.bvalue;
    Eb(1, 1) += (mu2 * E(1, 1) + lambda.value * tr) * output.bvalue;

    Eb(0, 1) += 2.0 * mu2 * E(0, 1) * output.bvalue;

    mu.bvalue = trE * output.bvalue;
    lambda.bvalue = 0.5 * tr * tr * output.bvalue;
  }

  ADScalar<ScalarType>& mu;
  ADScalar<ScalarType>& lambda;
  ADMat<SymMat<ScalarType, 2>>& EObj;
  ADScalar<ScalarType>& output;
};

template <class ScalarType>
KOKKOS_FUNCTION ADSymm2x2ADIsotropicEnergyExpr<ScalarType> SymmIsotropicEnergy(
    ADScalar<ScalarType>& mu, ADScalar<ScalarType>& lambda,
    ADMat<SymMat<ScalarType, 2>>& E, ADScalar<ScalarType>& output) {
  return ADSymm2x2ADIsotropicEnergyExpr<ScalarType>(mu, lambda, E, output);
}

template <class ScalarType>
class A2DSymm2x2IsotropicEnergyExpr
    : public A2DExpression<A2DSymm2x2IsotropicEnergyExpr<ScalarType>> {
 public:
  KOKKOS_FUNCTION A2DSymm2x2IsotropicEnergyExpr(
      const ScalarType& mu, const ScalarType& lambda,
      A2DMat<SymMat<ScalarType, 2>>& EObj, A2DScalar<ScalarType>& output)
      : mu(mu), lambda(lambda), EObj(EObj), output(output) {
    const SymMat<ScalarType, 2>& E = EObj.value();
    ScalarType tr = E(0, 0) + E(1, 1);
    ScalarType trE =
        E(0, 0) * E(0, 0) + E(1, 1) * E(1, 1) + 2.0 * E(0, 1) * E(0, 1);

    output.value = mu * trE + 0.5 * lambda * tr * tr;
  }

  KOKKOS_FUNCTION void reverse() {
    const SymMat<ScalarType, 2>& E = EObj.value();
    SymMat<ScalarType, 2>& Eb = EObj.bvalue();
    const ScalarType mu2 = 2.0 * mu;

    ScalarType tr = E(0, 0) + E(1, 1);
    Eb(0, 0) += (mu2 * E(0, 0) + lambda * tr) * output.bvalue;
    Eb(1, 1) += (mu2 * E(1, 1) + lambda * tr) * output.bvalue;

    Eb(0, 1) += 2.0 * mu2 * E(0, 1) * output.bvalue;
  }

  KOKKOS_FUNCTION void hforward() {
    const SymMat<ScalarType, 2>& E = EObj.value();
    ScalarType tr = E(0, 0) + E(1, 1);

    const SymMat<ScalarType, 2>& Ep = EObj.pvalue();
    ScalarType trd = Ep(0, 0) + Ep(1, 1);
    ScalarType trEd = 2.0 * (E(0, 0) * Ep(0, 0) + E(1, 1) * Ep(1, 1) +
                             2.0 * E(0, 1) * Ep(0, 1));

    output.pvalue = mu * trEd + lambda * tr * trd;
  }

  KOKKOS_FUNCTION void hreverse() {
    const SymMat<ScalarType, 2>& E = EObj.value();
    const ScalarType mu2 = 2.0 * mu;
    ScalarType tr = E(0, 0) + E(1, 1);

    const SymMat<ScalarType, 2>& Ep = EObj.pvalue();
    SymMat<ScalarType, 2>& Eh = EObj.hvalue();

    // by * (d^2y/dx^2 * px)
    ScalarType trp = Ep(0, 0) + Ep(1, 1);
    Eh(0, 0) += (mu2 * Ep(0, 0) + lambda * trp) * output.bvalue;
    Eh(1, 1) += (mu2 * Ep(1, 1) + lambda * trp) * output.bvalue;

    Eh(0, 1) += 2.0 * mu2 * Ep(0, 1) * output.bvalue;

    // hy * (dy/dx)
    Eh(0, 0) += (mu2 * E(0, 0) + lambda * tr) * output.hvalue;
    Eh(1, 1) += (mu2 * E(1, 1) + lambda * tr) * output.hvalue;

    Eh(0, 1) += 2.0 * mu2 * E(0, 1) * output.hvalue;
  }

  const ScalarType& mu;
  const ScalarType& lambda;
  A2DMat<SymMat<ScalarType, 2>>& EObj;
  A2DScalar<ScalarType>& output;
};

template <class ScalarType>
KOKKOS_FUNCTION A2DSymm2x2IsotropicEnergyExpr<ScalarType> SymmIsotropicEnergy(
    const ScalarType& mu, const ScalarType& lambda,
    A2DMat<SymMat<ScalarType, 2>>& E, A2DScalar<ScalarType>& output) {
  return A2DSymm2x2IsotropicEnergyExpr<ScalarType>(mu, lambda, E, output);
}

template <class ScalarType>
class A2DSymm2x2A2DIsotropicEnergyExpr
    : public A2DExpression<A2DSymm2x2A2DIsotropicEnergyExpr<ScalarType>> {
 public:
  KOKKOS_FUNCTION A2DSymm2x2A2DIsotropicEnergyExpr(
      A2DScalar<ScalarType>& mu, A2DScalar<ScalarType>& lambda,
      A2DMat<SymMat<ScalarType, 2>>& EObj, A2DScalar<ScalarType>& output)
      : mu(mu), lambda(lambda), EObj(EObj), output(output) {
    const SymMat<ScalarType, 2>& E = EObj.value();
    ScalarType tr = E(0, 0) + E(1, 1);
    ScalarType trE =
        E(0, 0) * E(0, 0) + E(1, 1) * E(1, 1) + 2.0 * E(0, 1) * E(0, 1);

    output.value = mu.value * trE + 0.5 * lambda.value * tr * tr;
  }

  KOKKOS_FUNCTION void reverse() {
    const SymMat<ScalarType, 2>& E = EObj.value();
    SymMat<ScalarType, 2>& Eb = EObj.bvalue();
    const ScalarType mu2 = 2.0 * mu.value;

    ScalarType tr = E(0, 0) + E(1, 1);
    ScalarType trE =
        E(0, 0) * E(0, 0) + E(1, 1) * E(1, 1) + 2.0 * E(0, 1) * E(0, 1);

    Eb(0, 0) += (mu2 * E(0, 0) + lambda.value * tr) * output.bvalue;
    Eb(1, 1) += (mu2 * E(1, 1) + lambda.value * tr) * output.bvalue;

    Eb(0, 1) += 2.0 * mu2 * E(0, 1) * output.bvalue;

    mu.bvalue = trE * output.bvalue;
    lambda.bvalue = 0.5 * tr * tr * output.bvalue;
  }

  KOKKOS_FUNCTION void hforward() {
    const SymMat<ScalarType, 2>& E = EObj.value();
    ScalarType tr = E(0, 0) + E(1, 1);
    ScalarType trE =
        E(0, 0) * E(0, 0) + E(1, 1) * E(1, 1) + 2.0 * E(0, 1) * E(0, 1);

    const SymMat<ScalarType, 2>& Ep = EObj.pvalue();
    ScalarType trd = Ep(0, 0) + Ep(1, 1);
    ScalarType trEd = 2.0 * (E(0, 0) * Ep(0, 0) + E(1, 1) * Ep(1, 1) +
                             +2.0 * E(0, 1) * Ep(0, 1));

    output.pvalue = mu.value * trEd + lambda.value * tr * trd +
                    mu.pvalue * trE + 0.5 * lambda.pvalue * tr * tr;
  }

  KOKKOS_FUNCTION void hreverse() {
    const SymMat<ScalarType, 2>& E = EObj.value();
    const ScalarType mu2 = 2.0 * mu.value;
    ScalarType tr = E(0, 0) + E(1, 1);
    ScalarType trE =
        E(0, 0) * E(0, 0) + E(1, 1) * E(1, 1) + 2.0 * E(0, 1) * E(0, 1);

    const SymMat<ScalarType, 2>& Ep = EObj.pvalue();
    SymMat<ScalarType, 2>& Eh = EObj.hvalue();

    // by * (d^2y/dx^2 * px)
    ScalarType trp = Ep(0, 0) + Ep(1, 1);
    Eh(0, 0) += (mu2 * Ep(0, 0) + lambda.value * trp) * output.bvalue;
    Eh(1, 1) += (mu2 * Ep(1, 1) + lambda.value * trp) * output.bvalue;

    Eh(0, 1) += 2.0 * mu2 * Ep(0, 1) * output.bvalue;

    // hy * (dy/dx)
    Eh(0, 0) += (mu2 * E(0, 0) + lambda.value * tr) * output.hvalue;
    Eh(1, 1) += (mu2 * E(1, 1) + lambda.value * tr) * output.hvalue;

    Eh(0, 1) += 2.0 * mu2 * E(0, 1) * output.hvalue;

    mu.hvalue += trE * output.hvalue;
    lambda.hvalue += 0.5 * tr * tr * output.hvalue;

    // account for Hessian blocks w.r.t. E and mu or lambda
    Eh(0, 0) +=
        (2.0 * E(0, 0) * mu.pvalue + tr * lambda.pvalue) * output.bvalue;
    Eh(1, 1) +=
        (2.0 * E(1, 1) * mu.pvalue + tr * lambda.pvalue) * output.bvalue;
    Eh(0, 1) += 4.0 * E(0, 1) * mu.pvalue * output.bvalue;

    mu.hvalue += (2.0 * (E(0, 0) * Ep(0, 0) + E(1, 1) * Ep(1, 1)) +
                  4.0 * E(0, 1) * Ep(0, 1)) *
                 output.bvalue;
    lambda.hvalue += tr * trp * output.bvalue;
  }

  A2DScalar<ScalarType>& mu;
  A2DScalar<ScalarType>& lambda;
  A2DMat<SymMat<ScalarType, 2>>& EObj;
  A2DScalar<ScalarType>& output;
};

template <class ScalarType>
KOKKOS_FUNCTION A2DSymm2x2A2DIsotropicEnergyExpr<ScalarType>
SymmIsotropicEnergy(A2DScalar<ScalarType>& mu, A2DScalar<ScalarType>& lambda,
                    A2DMat<SymMat<ScalarType, 2>>& E,
                    A2DScalar<ScalarType>& output) {
  return A2DSymm2x2A2DIsotropicEnergyExpr<ScalarType>(mu, lambda, E, output);
}

template <class ScalarType>
KOKKOS_FUNCTION void MatGreenStrain(const Mat<ScalarType, 2, 2>& Ux,
                                    SymMat<ScalarType, 2>& E) {
  E(0, 0) = Ux(0, 0) + 0.5 * (Ux(0, 0) * Ux(0, 0) + Ux(1, 0) * Ux(1, 0));
  E(1, 1) = Ux(1, 1) + 0.5 * (Ux(0, 1) * Ux(0, 1) + Ux(1, 1) * Ux(1, 1));

  E(0, 1) =
      0.5 * (Ux(0, 1) + Ux(1, 0) + Ux(0, 0) * Ux(0, 1) + Ux(1, 0) * Ux(1, 1));
}

template <typename ScalarType>
class ADMat2x2GreenStrainExpr
    : public ADExpression<ADMat2x2GreenStrainExpr<ScalarType>> {
 public:
  KOKKOS_FUNCTION ADMat2x2GreenStrainExpr(ADMat<Mat<ScalarType, 2, 2>>& UxObj,
                                          ADMat<SymMat<ScalarType, 2>>& EObj)
      : UxObj(UxObj), EObj(EObj) {
    const Mat<ScalarType, 2, 2>& Ux = UxObj.value();
    SymMat<ScalarType, 2>& E = EObj.value();
    E(0, 0) = Ux(0, 0) + 0.5 * (Ux(0, 0) * Ux(0, 0) + Ux(1, 0) * Ux(1, 0));
    E(1, 1) = Ux(1, 1) + 0.5 * (Ux(0, 1) * Ux(0, 1) + Ux(1, 1) * Ux(1, 1));

    E(0, 1) =
        0.5 * (Ux(0, 1) + Ux(1, 0) + Ux(0, 0) * Ux(0, 1) + Ux(1, 0) * Ux(1, 1));
  }

  KOKKOS_FUNCTION void forward() {
    const Mat<ScalarType, 2, 2>& Ux = UxObj.value();
    const Mat<ScalarType, 2, 2>& Uxd = UxObj.bvalue();
    SymMat<ScalarType, 2>& Ed = EObj.bvalue();

    Ed(0, 0) = Uxd(0, 0) + Ux(0, 0) * Uxd(0, 0) + Ux(1, 0) * Uxd(1, 0);
    Ed(1, 1) = Uxd(1, 1) + Ux(0, 1) * Uxd(0, 1) + Ux(1, 1) * Uxd(1, 1);

    Ed(0, 1) = 0.5 * (Uxd(0, 1) + Uxd(1, 0) + Ux(0, 0) * Uxd(0, 1) +
                      Ux(1, 0) * Uxd(1, 1) + Uxd(0, 0) * Ux(0, 1) +
                      Uxd(1, 0) * Ux(1, 1));
  }

  KOKKOS_FUNCTION void reverse() {
    const Mat<ScalarType, 2, 2>& Ux = UxObj.value();
    const SymMat<ScalarType, 2>& Eb = EObj.bvalue();
    Mat<ScalarType, 2, 2>& Uxb = UxObj.bvalue();

    // Uxb = (I + Ux) * Eb
    Uxb(0, 0) += (Ux(0, 0) + 1.0) * Eb(0, 0) + 0.5 * Ux(0, 1) * Eb(0, 1);
    Uxb(0, 1) += 0.5 * (Ux(0, 0) + 1.0) * Eb(0, 1) + Ux(0, 1) * Eb(1, 1);

    Uxb(1, 0) += Ux(1, 0) * Eb(0, 0) + 0.5 * (Ux(1, 1) + 1.0) * Eb(0, 1);
    Uxb(1, 1) += 0.5 * Ux(1, 0) * Eb(0, 1) + (Ux(1, 1) + 1.0) * Eb(1, 1);
  }

  ADMat<Mat<ScalarType, 2, 2>>& UxObj;
  ADMat<SymMat<ScalarType, 2>>& EObj;
};

template <typename ScalarType>
KOKKOS_FUNCTION ADMat2x2GreenStrainExpr<ScalarType> MatGreenStrain(
    ADMat<Mat<ScalarType, 2, 2>>& Ux, ADMat<SymMat<ScalarType, 2>>& E) {
  return ADMat2x2GreenStrainExpr<ScalarType>(Ux, E);
}

template <typename ScalarType>
class A2DMat2x2GreenStrainExpr
    : public A2DExpression<A2DMat2x2GreenStrainExpr<ScalarType>> {
 public:
  KOKKOS_FUNCTION A2DMat2x2GreenStrainExpr(A2DMat<Mat<ScalarType, 2, 2>>& UxObj,
                                           A2DMat<SymMat<ScalarType, 2>>& EObj)
      : UxObj(UxObj), EObj(EObj) {
    const Mat<ScalarType, 2, 2>& Ux = UxObj.value();
    SymMat<ScalarType, 2>& E = EObj.value();
    E(0, 0) = Ux(0, 0) + 0.5 * (Ux(0, 0) * Ux(0, 0) + Ux(1, 0) * Ux(1, 0));
    E(1, 1) = Ux(1, 1) + 0.5 * (Ux(0, 1) * Ux(0, 1) + Ux(1, 1) * Ux(1, 1));

    E(0, 1) =
        0.5 * (Ux(0, 1) + Ux(1, 0) + Ux(0, 0) * Ux(0, 1) + Ux(1, 0) * Ux(1, 1));
  }

  KOKKOS_FUNCTION void reverse() {
    const Mat<ScalarType, 2, 2>& Ux = UxObj.value();
    const SymMat<ScalarType, 2>& Eb = EObj.bvalue();
    Mat<ScalarType, 2, 2>& Uxb = UxObj.bvalue();

    // Uxb = (I + Ux) * Eb
    Uxb(0, 0) += (Ux(0, 0) + 1.0) * Eb(0, 0) + 0.5 * Ux(0, 1) * Eb(0, 1);
    Uxb(0, 1) += 0.5 * (Ux(0, 0) + 1.0) * Eb(0, 1) + Ux(0, 1) * Eb(1, 1);

    Uxb(1, 0) += Ux(1, 0) * Eb(0, 0) + 0.5 * (Ux(1, 1) + 1.0) * Eb(0, 1);
    Uxb(1, 1) += 0.5 * Ux(1, 0) * Eb(0, 1) + (Ux(1, 1) + 1.0) * Eb(1, 1);
  }

  KOKKOS_FUNCTION void hforward() {
    const Mat<ScalarType, 2, 2>& Ux = UxObj.value();

    const Mat<ScalarType, 2, 2>& Uxp = UxObj.pvalue();
    SymMat<ScalarType, 2>& Ep = EObj.pvalue();
    Ep(0, 0) = Uxp(0, 0) + Ux(0, 0) * Uxp(0, 0) + Ux(1, 0) * Uxp(1, 0);
    Ep(1, 1) = Uxp(1, 1) + Ux(0, 1) * Uxp(0, 1) + Ux(1, 1) * Uxp(1, 1);

    Ep(0, 1) = 0.5 * (Uxp(0, 1) + Uxp(1, 0) + Ux(0, 0) * Uxp(0, 1) +
                      Ux(1, 0) * Uxp(1, 1) + Uxp(0, 0) * Ux(0, 1) +
                      Uxp(1, 0) * Ux(1, 1));
  }

  KOKKOS_FUNCTION void hreverse() {
    const Mat<ScalarType, 2, 2>& Eb = EObj.bvalue();
    const Mat<ScalarType, 2, 2>& Ux = UxObj.value();

    const Mat<ScalarType, 2, 2>& Uxp = UxObj.pvalue();
    const SymMat<ScalarType, 2>& Eh = EObj.hvalue();
    Mat<ScalarType, 2, 2>& Uxh = UxObj.hvalue();

    Uxh(0, 0) += Uxp(0, 0) * Eb(0, 0) + 0.5 * Uxp(0, 1) * Eb(0, 1);
    Uxh(0, 1) += 0.5 * Uxp(0, 0) * Eb(0, 1) + Uxp(0, 1) * Eb(1, 1);

    Uxh(1, 0) += Uxp(1, 0) * Eb(0, 0) + 0.5 * Uxp(1, 1) * Eb(0, 1);
    Uxh(1, 1) += 0.5 * Uxp(1, 0) * Eb(0, 1) + Uxp(1, 1) * Eb(1, 1);

    Uxh(0, 0) += (Ux(0, 0) + 1.0) * Eh(0, 0) + 0.5 * Ux(0, 1) * Eh(0, 1);
    Uxh(0, 1) += 0.5 * (Ux(0, 0) + 1.0) * Eh(0, 1) + Ux(0, 1) * Eh(1, 1);

    Uxh(1, 0) += Ux(1, 0) * Eh(0, 0) + 0.5 * (Ux(1, 1) + 1.0) * Eh(0, 1);
    Uxh(1, 1) += 0.5 * Ux(1, 0) * Eh(0, 1) + (Ux(1, 1) + 1.0) * Eh(1, 1);
  }

  A2DMat<Mat<ScalarType, 2, 2>>& UxObj;
  A2DMat<SymMat<ScalarType, 2>>& EObj;
};

template <typename ScalarType>
KOKKOS_FUNCTION A2DMat2x2GreenStrainExpr<ScalarType> MatGreenStrain(
    A2DMat<Mat<ScalarType, 2, 2>>& Ux, A2DMat<SymMat<ScalarType, 2>>& E) {
  return A2DMat2x2GreenStrainExpr<ScalarType>(Ux, E);
}

template <class ScalarType>
KOKKOS_FUNCTION void MatLinearGreenStrain(const Mat<ScalarType, 2, 2>& Ux,
                                          SymMat<ScalarType, 2>& E) {
  E(0, 0) = Ux(0, 0);
  E(1, 1) = Ux(1, 1);

  E(0, 1) = 0.5 * (Ux(0, 1) + Ux(1, 0));
}

template <typename ScalarType>
class ADMat2x2LinearGreenStrainExpr
    : public ADExpression<ADMat2x2LinearGreenStrainExpr<ScalarType>> {
 public:
  KOKKOS_FUNCTION ADMat2x2LinearGreenStrainExpr(
      ADMat<Mat<ScalarType, 2, 2>>& UxObj, ADMat<SymMat<ScalarType, 2>>& EObj)
      : UxObj(UxObj), EObj(EObj) {
    const Mat<ScalarType, 2, 2>& Ux = UxObj.value();
    SymMat<ScalarType, 2>& E = EObj.value();
    E(0, 0) = Ux(0, 0);
    E(1, 1) = Ux(1, 1);

    E(0, 1) = 0.5 * (Ux(0, 1) + Ux(1, 0));
  }

  KOKKOS_FUNCTION void forward() {
    const Mat<ScalarType, 2, 2>& Uxb = UxObj.bvalue();
    SymMat<ScalarType, 2>& Eb = EObj.bvalue();

    Eb(0, 0) = Uxb(0, 0);
    Eb(1, 1) = Uxb(1, 1);

    Eb(0, 1) = 0.5 * (Uxb(0, 1) + Uxb(1, 0));
  }

  KOKKOS_FUNCTION void reverse() {
    const SymMat<ScalarType, 2>& Eb = EObj.bvalue();
    Mat<ScalarType, 2, 2>& Uxb = UxObj.bvalue();

    // Uxb = (I + Ux) * Eb
    Uxb(0, 0) += Eb(0, 0);
    Uxb(0, 1) += 0.5 * Eb(0, 1);

    Uxb(1, 0) += 0.5 * Eb(0, 1);
    Uxb(1, 1) += Eb(1, 1);
  }

  ADMat<Mat<ScalarType, 2, 2>>& UxObj;
  ADMat<SymMat<ScalarType, 2>>& EObj;
};

template <typename ScalarType>
KOKKOS_FUNCTION ADMat2x2LinearGreenStrainExpr<ScalarType> MatLinearGreenStrain(
    ADMat<Mat<ScalarType, 2, 2>>& Ux, ADMat<SymMat<ScalarType, 2>>& E) {
  return ADMat2x2LinearGreenStrainExpr<ScalarType>(Ux, E);
}

template <typename ScalarType>
class A2DMat2x2LinearGreenStrainExpr
    : public A2DExpression<A2DMat2x2LinearGreenStrainExpr<ScalarType>> {
 public:
  KOKKOS_FUNCTION A2DMat2x2LinearGreenStrainExpr(
      A2DMat<Mat<ScalarType, 2, 2>>& UxObj, A2DMat<SymMat<ScalarType, 2>>& EObj)
      : UxObj(UxObj), EObj(EObj) {
    const Mat<ScalarType, 2, 2>& Ux = UxObj.value();
    SymMat<ScalarType, 2>& E = EObj.value();
    E(0, 0) = Ux(0, 0);
    E(1, 1) = Ux(1, 1);

    E(0, 1) = 0.5 * (Ux(0, 1) + Ux(1, 0));
  }

  KOKKOS_FUNCTION void reverse() {
    const SymMat<ScalarType, 2>& Eb = EObj.bvalue();
    Mat<ScalarType, 2, 2>& Uxb = UxObj.bvalue();

    // Uxb = Eb
    Uxb(0, 0) += Eb(0, 0);
    Uxb(0, 1) += 0.5 * Eb(0, 1);

    Uxb(1, 0) += 0.5 * Eb(0, 1);
    Uxb(1, 1) += Eb(1, 1);
  }

  KOKKOS_FUNCTION void hforward() {
    const Mat<ScalarType, 2, 2>& Uxp = UxObj.pvalue();
    SymMat<ScalarType, 2>& Ep = EObj.pvalue();

    Ep(0, 0) = Uxp(0, 0);
    Ep(1, 1) = Uxp(1, 1);

    Ep(0, 1) = 0.5 * (Uxp(0, 1) + Uxp(1, 0));
  }

  KOKKOS_FUNCTION void hreverse() {
    const SymMat<ScalarType, 2>& Eh = EObj.hvalue();
    Mat<ScalarType, 2, 2>& Uxh = UxObj.hvalue();

    Uxh(0, 0) += Eh(0, 0);
    Uxh(0, 1) += 0.5 * Eh(0, 1);

    Uxh(1, 0) += 0.5 * Eh(0, 1);
    Uxh(1, 1) += Eh(1, 1);
  }

  A2DMat<Mat<ScalarType, 2, 2>>& UxObj;
  A2DMat<SymMat<ScalarType, 2>>& EObj;
};

template <typename ScalarType>
KOKKOS_FUNCTION A2DMat2x2LinearGreenStrainExpr<ScalarType> MatLinearGreenStrain(
    A2DMat<Mat<ScalarType, 2, 2>>& Ux, A2DMat<SymMat<ScalarType, 2>>& E) {
  return A2DMat2x2LinearGreenStrainExpr<ScalarType>(Ux, E);
}

}  // namespace A2D

#endif  // A2D_MAT_OPS_2D_H
