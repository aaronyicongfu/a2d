#include <iostream>
#include <memory>

#include "block_numeric.h"
#include "multiarray.h"

using namespace A2D;
using namespace std;

int main(int argc, char* argv[]) {
  using T = std::complex<double>;
  using A_Layout = FLayout<3>;
  using Ainv_Layout = CLayout<3>;

  auto A = MultiArray<T, A_Layout>(A_Layout(3));
  auto Ainv = MultiArray<T, Ainv_Layout>(Ainv_Layout(3));
  A(0, 0) = T(0.95213344, 0.67950085);
  A(0, 1) = T(0.55520316, 0.12368274);
  A(0, 2) = T(0.49634073, 0.38212409);
  A(1, 0) = T(0.25232854, 0.3560515);
  A(1, 1) = T(0.49572308, 0.17893203);
  A(1, 2) = T(0.9012511, 0.37591567);
  A(2, 0) = T(0.74097586, 0.82404429);
  A(2, 1) = T(0.59794044, 0.73940968);
  A(2, 2) = T(0.13228835, 0.80532406);

  printf("A:\n");
  for (int i = 0; i != 3; i++) {
    for (int j = 0; j != 3; j++) {
      printf("%8f+%8fj   ", A(i, j).real(), A(i, j).imag());
    }
    printf("\n");
  }

  blockPseudoInverse<T, 3>(A, Ainv);

  printf("inv(A):\n");
  for (int i = 0; i != 3; i++) {
    for (int j = 0; j != 3; j++) {
      printf("%8f+%8fj   ", Ainv(i, j).real(), Ainv(i, j).imag());
    }
    printf("\n");
  }
}