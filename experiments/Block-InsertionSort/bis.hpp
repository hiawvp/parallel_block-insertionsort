#ifndef BIS
#define BIS

#include "common/utils.hpp"
#include <cstddef>
#include <cstring>
#include <iostream>
#include <math.h>
#include <type_traits>

// TODO: Aceptar una funcion comp
// implementaciones personales de bis
namespace bis {
template <typename T> void insertion_sort(T *A, size_t l, size_t r);
const int DEFAULT_BLOCK_SIZE = 16;
const int SEQUENTIAL_THRESHOLD = 2 << 15;

namespace internal {
template <typename T>
void insert_second_block(size_t l, size_t r, size_t current_block_size, T *S,
                         T *E) {
  size_t dt = current_block_size;
  long j = r - current_block_size - 1;
  T key;
  for (long i = current_block_size - 1; i >= 0; i--) {
    key = E[i];
    while (j >= (long)l && key < S[j]) {
      S[j + dt] = S[j];
      j--;
    }
    S[j + dt] = key;
    dt--;
  }
  return;
}

template <typename T>
void insert_all_but_second_block(size_t l, size_t r, size_t block_size, T *S) {
  // inicio del tercer bloque
  size_t ini = l + 2 * block_size;
  // fin del tercer bloque
  size_t end = ini + block_size;
  long j = l + block_size - 1;
  size_t dt;
  T key;
  while (ini < r) {
    if (end > r) {
      end = r;
      dt = end - ini;
    } else {
      dt = block_size;
    }
    for (long i = end - 1; i >= (long)ini; i--) {
      key = S[i];
      while (j >= (long)l && key < S[j]) {
        S[j + dt] = S[j];
        j--;
      }
      S[j + dt] = key;
      dt--;
    }
    j = ini - 1;
    ini = end;
    end += block_size;
  }
  return;
}

template <typename T>
void rec_bis_clean(T *S, T *E, size_t l, size_t r, size_t e_start, int &k) {
  if (l >= r) {
    return;
  }
  size_t n = r - l;
  if (n <= (size_t)k) {
    return insertion_sort(S, l, r);
  }
  size_t kpow = ceil(log(n) / log(k));
  size_t current_block_size = pow(k, kpow - 1);
  size_t next_block_size = pow(k, kpow - 2);
  for (size_t i = l, e_idx = e_start; i < r;
       i += current_block_size, e_idx += next_block_size) {
    size_t ir = std::min(i + current_block_size, r);
    internal::rec_bis_clean(S, E, i, ir, e_idx, k);
  }
  size_t second_block_start = l + current_block_size;
  size_t second_block_end = second_block_start + current_block_size;
  size_t copied_block_size = current_block_size;
  if (second_block_end > r) {
    second_block_end = r;
    copied_block_size = second_block_end - second_block_start;
  }
  if (second_block_end >= second_block_start) {
    memcpy(E + e_start, S + second_block_start, copied_block_size * sizeof(T));
  }
  insert_all_but_second_block(l, r, current_block_size, S);
  insert_second_block(l, r, copied_block_size, S, E + e_start);
}

} // namespace internal

template <typename T>
void block_insertion_sort(T *S, size_t n, size_t k = DEFAULT_BLOCK_SIZE) {
  if (n <= k) {
    return insertion_sort(S, 0, n);
  }
  for (size_t l = 0; l < n; l += k) {
    size_t r = std::min(l + k, n);
    insertion_sort(S, l, r);
  }
  size_t current_block_size = k;
  size_t starting_point, ending_point;
  size_t sorted_block_start_idx, sorted_block_end_idx;
  size_t delta;
  long j;

  size_t dt = std::log(n) / std::log(k);
  dt = std::pow(k, dt);
  size_t potBIS = log2(k);

  if (dt == n) {
    dt >>= potBIS;
  }
  T *E = new T[dt];
  T key;
  while (current_block_size < n) {
    starting_point = 0;
    while (starting_point < n) {
      // pasamos de un bloque de tamanio k^i a uno de tamanio k^(i+1)
      ending_point = starting_point + (k * current_block_size) - 1;
      // si la entrada no es potencia de k
      if (ending_point >= n) {
        ending_point = n - 1;
      }

      // inicio del segundo bloque
      sorted_block_start_idx = starting_point + current_block_size;
      // inicio del tercer bloque
      sorted_block_end_idx = sorted_block_start_idx + current_block_size;
      size_t copied_block_size = current_block_size;
      if (sorted_block_end_idx > ending_point + 1) {
        sorted_block_end_idx = ending_point + 1;
        copied_block_size = sorted_block_end_idx - sorted_block_start_idx;
      }
      //
      if (sorted_block_end_idx >= sorted_block_start_idx) {
        memcpy(E, S + sorted_block_start_idx, copied_block_size * sizeof(T));
      }

      // partimos insertando el tercer bloque...
      sorted_block_start_idx = starting_point + 2 * current_block_size;
      sorted_block_end_idx = sorted_block_start_idx + current_block_size;
      j = starting_point + current_block_size - 1;

      // save the second block in E[]

      while (sorted_block_start_idx <= ending_point) {
        if (sorted_block_end_idx > ending_point + 1) {
          sorted_block_end_idx = ending_point + 1;
          delta = sorted_block_end_idx - sorted_block_start_idx;
        } else {
          delta = current_block_size;
        }

        for (long i = sorted_block_end_idx - 1;
             i >= (long)sorted_block_start_idx; i--) {
          key = S[i];
          while (j >= (long)starting_point && key < S[j]) {
            S[j + delta] = S[j];
            j--;
          }
          S[j + delta] = key;
          delta--;
        }
        j = sorted_block_start_idx - 1;
        sorted_block_start_idx = sorted_block_end_idx;
        sorted_block_end_idx += current_block_size;
      }

      delta = copied_block_size;
      j = ending_point - copied_block_size;

      // aqui se inserta el bloque guardado en E
      for (int i = copied_block_size - 1; i >= 0; i--) {
        key = E[i];
        while (j >= (long)starting_point && key < S[j]) {
          S[j + delta] = S[j];
          j--;
        }
        S[j + delta] = key;
        delta--;
      }
      starting_point = ending_point + 1;
    }
    current_block_size *= k;
  }
  delete[] E;
  return;
}

template <typename T>
void rec_block_insertion_sort(T *S, size_t l, size_t r, int k) {
  size_t n = r - l;
  if (l >= r) {
    return;
  }
  if (r - l <= (size_t)k) {
    return insertion_sort(S, l, r);
  }
  // todo: compute required aux mem size;
  size_t dt = std::log(n) / std::log(k);
  dt = std::pow(k, dt);
  size_t potBIS = log2(k);
  if (dt == n) {
    dt >>= potBIS;
  }
  T *E = new T[dt];
  internal::rec_bis_clean(S, E, l, r, 0, k);
  delete[] E;
}

// template <typename T> void insertion_sort(T *A, size_t l, size_t r) {
//   T key;
//   long j;
//   for (size_t i = l + 1; i < r; i++) {
//     key = A[i];
//     j = i - 1;
//     while (j >= l && key < A[j]) {
//       A[j + 1] = A[j];
//       j = j - 1;
//     }
//     A[j + 1] = key;
//   }
// }
template <typename T> void insertion_sort(T *A, size_t l, size_t r) {
  if (r <= l + 1) {
    return;
  }

  for (size_t i = l + 1; i < r; i++) {
    T key = A[i];
    long j = i - 1;

    while (j >= static_cast<long>(l) && key < A[static_cast<size_t>(j)]) {
      A[static_cast<size_t>(j + 1)] = A[static_cast<size_t>(j)];
      j--;
    }
    A[static_cast<size_t>(j + 1)] = key;
  }
}

} // namespace bis
#endif /* BASIC_DRF_H_ */
