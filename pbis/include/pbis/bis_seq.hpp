#ifndef BIS
#define BIS
#include <cmath>
#include <cstddef>
#include <cstring>
#include <functional>
#include <iterator>

namespace bis {
namespace __internal {
template <typename T, typename Compare>
void insertion_sort(T *A, size_t l, size_t r, Compare comp);
}

const int SEQUENTIAL_THRESHOLD = 2 << 15;
const int SEQ_DEFAULT_BLOCK_SIZE = 16;
const int PAR_DEFAULT_BLOCK_SIZE = 8;
const int MWM_SWITCH_STRAT = 0;
namespace sequential_sorter {

template <typename T, typename Compare>
void block_insertion_sort(T *S, size_t n, Compare comp,
                          size_t k = SEQ_DEFAULT_BLOCK_SIZE) {
  if (n <= k) {
    return bis::__internal::insertion_sort(S, 0, n, comp);
  }
  for (size_t l = 0; l < n; l += k) {
    size_t r = std::min(l + k, n);
    bis::__internal::insertion_sort(S, l, r, comp);
  }
  size_t current_block_size = k;
  size_t starting_point, ending_point;
  size_t sorted_block_start_idx, sorted_block_end_idx;
  size_t delta;
  long j;

  size_t dt = std::log(n) / std::log(k);
  dt = std::pow(k, dt);
  size_t potBIS = std::log2(k);

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
        std::memcpy(E, S + sorted_block_start_idx,
                    copied_block_size * sizeof(T));
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
          while (j >= (long)starting_point && comp(key, S[j])) {
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
        while (j >= (long)starting_point && comp(key, S[j])) {
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
} // namespace sequential_sorter

//
namespace __internal {
template <typename T, typename Compare>
void insertion_sort(T *A, size_t l, size_t r, Compare comp) {
  if (r <= l + 1) {
    return;
  }

  for (size_t i = l + 1; i < r; i++) {
    T key = A[i];
    long j = i - 1;

    while (j >= static_cast<long>(l) && comp(key, A[static_cast<size_t>(j)])) {
      A[static_cast<size_t>(j + 1)] = A[static_cast<size_t>(j)];
      j--;
    }
    A[static_cast<size_t>(j + 1)] = key;
  }
}
} // namespace __internal
} // namespace bis

#endif /* BASIC_DRF_H_ */
