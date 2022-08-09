#ifndef A2D_MULTI_ARRAY_H
#define A2D_MULTI_ARRAY_H

#include <cstddef>
#include <cstdlib>
#include <numeric>
#include <type_traits>

#include "a2dobjs.h"

namespace A2D {

/*
  Unpack the argument list to get the extent of indices
*/
template <int r, index_t... dims>
struct ___get_extent;

template <index_t dim0, index_t... dims>
struct ___get_extent<0, dim0, dims...> {
  static const index_t extent = dim0;
};

template <int r, index_t dim0, index_t... dims>
struct ___get_extent<r, dim0, dims...> {
  static const index_t extent = ___get_extent<r - 1, dims...>::extent;
};

template <index_t... dims>
struct __get_size;

template <index_t dim0, index_t... dims>
struct __get_size<dim0, dims...> {
  static const index_t size = dim0 * __get_size<dims...>::size;
};

template <index_t dim0>
struct __get_size<dim0> {
  static const index_t size = dim0;
};

/*
  Fortran ordering

  Given an entry (i, j, k) in a multi-dimensional array, the fortran ordering
  stores entry A(i, j, k) in the location i + dim1 * (j + dim2 * k).
*/
template <index_t... dims>
class FLayout {
 public:
  FLayout() {}
  FLayout(const index_t dim0) : dim0(dim0), static_extents{dims...} {
    size = dim0;
    for (index_t i = 1; i < get_rank(); i++) {
      size *= get_extent(i);
    }
  }
  FLayout(const FLayout<dims...>& src)
      : dim0(src.dim0), static_extents{dims...} {
    size = dim0;
    for (index_t i = 1; i < get_rank(); i++) {
      size *= get_extent(i);
    }
  }
  const index_t dim0;
  static const index_t rank = sizeof...(dims) + 1;
  static const index_t get_rank() { return rank; }

  // Get the size of the array required given the first dimension
  static index_t get_size(index_t dim) {
    return dim * __get_size<dims...>::size;
  }

  const index_t get_extent(index_t index) const {
    if (index == 0) {
      return dim0;
    }
    return static_extents[index - 1];
  }

  template <class Idx, class... IdxType>
  const index_t compute_index(Idx i1, IdxType... idx) const {
    return i1 + dim0 * __compute_index<0>(idx...);
  }

  const index_t get_size() { return size; }

 private:
  index_t size;
  index_t static_extents[rank - 1];

  template <int r, class Idx, class... IdxType>
  static const index_t __compute_index(Idx i, IdxType... idx) {
    return i +
           ___get_extent<r, dims...>::extent * __compute_index<r + 1>(idx...);
  }

  template <int r, class Idx>
  static const index_t __compute_index(Idx i) {
    return i;
  }
};

/*
  C ordering

  Given an entry (i, j, k) in a multi-dimensional array, the c ordering
  stores entry A(i, j, k) in the location (i * dim0 + j) * dim1 + k.
*/
template <index_t... dims>
class CLayout {
 public:
  CLayout() {}
  CLayout(const index_t dim0) : dim0(dim0), static_extents{dims...} {
    size = dim0;
    for (index_t i = 1; i < get_rank(); i++) {
      size *= get_extent(i);
    }
  }
  CLayout(const CLayout<dims...>& src)
      : dim0(src.dim0), static_extents{dims...} {
    size = dim0;
    for (index_t i = 1; i < get_rank(); i++) {
      size *= get_extent(i);
    }
  }

  index_t dim0;
  static const index_t rank = sizeof...(dims) + 1;  // total number of
                                                    // dimensions (fixed and
                                                    // variable)
  static const index_t get_rank() { return rank; }

  // Get the size of the array required given the first dimension
  static index_t get_size(index_t dim) {
    return dim * __get_size<dims...>::size;
  }

  const index_t get_extent(index_t index) const {
    if (index == 0) {
      return dim0;
    }
    return static_extents[index - 1];
  }

  template <class Idx, class... IdxType>
  const index_t compute_index(Idx i1, IdxType... idx) const {
    return __compute_index<0>(i1, idx...);
  }

  const index_t get_size() { return size; }

 private:
  index_t size;
  index_t static_extents[rank - 1];

  template <int r, class Idx, class... IdxType>
  static const index_t __compute_index(const index_t index, Idx i,
                                       IdxType... idx) {
    return __compute_index<r + 1>(index * ___get_extent<r, dims...>::extent + i,
                                  idx...);
  }

  template <int r, class Idx>
  static const index_t __compute_index(const index_t index, Idx i) {
    return index * ___get_extent<r, dims...>::extent + i;
  }
};

/*
  A multi-dimensional array
*/
template <typename T, class Layout>
class MultiArray {
 public:
  typedef T type;

  MultiArray() { data = nullptr; }

  MultiArray(Layout layout, T* data_ = NULL) : layout(layout), data(data_) {
    if (data) {
      data_owner = false;
    } else {
      data_owner = true;
      data = new T[layout.get_size()];
      zero();
    }
  }
  MultiArray(const MultiArray<T, Layout>& src)
      : layout(src.layout), data(src.data) {
    data_owner = false;
  }
  ~MultiArray() {
    if (data_owner) {
      delete[] data;
    }
  }

  Layout layout;
  T* data;

  T* get_data() { return data; }

  /*
    Non-constant access to array elements
  */
  template <class... IdxType>
  T& operator()(IdxType... idx) {
    return data[layout.compute_index(idx...)];
  }

  /*
    Constant access to array elements
  */
  template <class... IdxType>
  const T& operator()(IdxType... idx) const {
    return data[layout.compute_index(idx...)];
  }

  /*
    Get the rank of the the multi-dimensional array
  */
  const index_t get_rank() const { return layout.get_rank(); }

  /*
    Get the extent of one of the dimensions
  */
  const index_t extent(index_t index) const { return layout.get_extent(index); }

  /*
    Zero all elements in the array
  */
  void zero() {
    const index_t len = layout.get_size();
    std::fill(data, data + len, T(0));
  }

  /*
    Fill all the values in the array with the specified value
  */
  void fill(T value) {
    const index_t len = layout.get_size();
    std::fill(data, data + len, value);
  }

  /*
    Copy elements from the source to this vector
  */
  void copy(MultiArray<T, Layout>& src) {
    const index_t len = layout.get_size();
    std::copy(src.data, src.data + len, data);
  }

  /*
    Copy elements from the source to this vector
  */
  void scale(T alpha) {
    const index_t len = layout.get_size();
    for (index_t i = 0; i < len; i++) {
      data[i] *= alpha;
    }
  }

  /*
    Duplicate the array
  */
  MultiArray<T, Layout>* duplicate() {
    const index_t len = layout.get_size();
    MultiArray<T, Layout>* array = new MultiArray<T, Layout>(layout);
    std::copy(data, data + len, array->data);
    return array;
  }

  /*
    Set a random seed for the data
  */
  void random(T lower = -1.0, T upper = 1.0) {
    const index_t len = layout.get_size();
    for (index_t i = 0; i < len; i++) {
      data[i] =
          lower + ((upper - lower) * (1.0 * std::rand())) / (1.0 * RAND_MAX);
    }
  }

  /*
    Take the dot product with the source vector data
  */
  T dot(MultiArray<T, Layout>& src) {
    const index_t len = layout.get_size();
    return std::inner_product(data, data + len, src.data, T(0));
  }

  /*
    Norm of the array
  */
  T norm() {
    const index_t len = layout.get_size();
    return std::sqrt(std::inner_product(data, data + len, data, T(0)));
  }

  /*
    Axpy: this = alpha * x + this
  */
  void axpy(T alpha, MultiArray<T, Layout>& x) {
    const index_t len = layout.get_size();
    for (index_t i = 0; i < len; i++) {
      data[i] += alpha * x.data[i];
    }
  }

  /*
    Axpby: this = alpha * x + beta * this
  */
  void axpby(T alpha, T beta, MultiArray<T, Layout>& x) {
    const index_t len = layout.get_size();
    for (index_t i = 0; i < len; i++) {
      data[i] = alpha * x.data[i] + beta * data[i];
    }
  }

 private:
  bool data_owner;
};

/**
 * @brief Note: MultiArraySlice only works for CLayout
 */
template <typename T, index_t... dims>
class MultiArraySlice {
 public:
  template <class IdxType>
  MultiArraySlice(MultiArray<T, CLayout<dims...>>& array, const IdxType idx) {
    data = &array.data[array.layout.get_size(idx)];
  }

  /*
    Non-constant access to array elements
  */
  template <class... IdxType>
  T& operator()(IdxType... idx) {
    return data[__compute_index<0, IdxType...>(0, idx...)];
  }

  /*
    Constant access to array elements
  */
  template <class... IdxType>
  const T& operator()(IdxType... idx) const {
    return data[__compute_index<0, IdxType...>(0, idx...)];
  }

  /*
    Zero the entries of the slice
  */
  void zero() {
    for (index_t i = 0; i < size; i++) {
      data[i] = 0.0;
    }
  }

  /**
   * Get the underlying data pointer
   */
  T* get_data() { return data; }

 private:
  T* data;

  static const index_t size = __get_size<dims...>::size;

  template <int r, class Idx, class... IdxType>
  static const index_t __compute_index(const index_t index, Idx i,
                                       IdxType... idx) {
    return __compute_index<r + 1>(index * ___get_extent<r, dims...>::extent + i,
                                  idx...);
  }

  template <int r, class Idx>
  static const index_t __compute_index(const index_t index, Idx i) {
    return index * ___get_extent<r, dims...>::extent + i;
  }
};

template <typename T, typename IdxType, index_t... ldims>
MultiArraySlice<T, ldims...> MakeSlice(MultiArray<T, CLayout<ldims...>>& array,
                                       IdxType idx) {
  return MultiArraySlice<T, ldims...>(array, idx);
}

template <typename T, typename IdxType, index_t dim0, index_t... ldims>
MultiArraySlice<T, ldims...> MakeSlice(
    MultiArray<T, CLayout<dim0, ldims...>>& array, IdxType idx0, IdxType idx1) {
  return MultiArraySlice<T, ldims...>(array, idx0, idx1);
}

}  // namespace A2D

#endif  // A2D_MULTI_ARRAY_H
