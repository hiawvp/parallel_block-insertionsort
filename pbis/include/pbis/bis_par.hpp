
#ifndef PBIS
#define PBIS
#include "./bis_seq.hpp"
#include "execution"
#include "omp.h"
#include <cstdio>
#include <parallel/algorithm>

namespace bis {

namespace parallel_sorter {

namespace {

template <typename T, typename Compare>
void multiway_merge_output(T *arr, T *output, size_t start, size_t end,
                           size_t step, Compare comp) {
  // #if defined(_USETBB)
  size_t n = end - start;
  int nseqs = (n + step - 1) / step;

  // omp_set_num_threads(p);
  // tbb::global_control c(tbb::global_control::max_allowed_parallelism, p);
  std::vector<std::pair<T *, T *>> seqs;

  size_t r, l;
  for (int i = 0; i < nseqs; ++i) {
    l = start + (i * step);
    r = std::min(l + step, end);
    std::pair<T *, T *> pair = std::make_pair(arr + l, arr + r);
    seqs.push_back(pair);
  }
  __gnu_parallel::multiway_merge(seqs.begin(), seqs.end(), output, n, comp);

#if defined(_USETBB)
  std::copy(std::execution::par, output, output + n, arr + start);
#else
  std::copy(output, output + n, arr + start);
#endif // _USETBB
  return;
}

template <typename T, typename Compare>
void insert_second_block(size_t l, size_t r, size_t current_block_size, T *S,
                         T *E, Compare comp) {
  size_t dt = current_block_size;
  long j = r - current_block_size - 1;
  T key;
  for (long i = current_block_size - 1; i >= 0; i--) {
    key = E[i];
    while (j >= (long)l && comp(key, S[j])) {
      S[j + dt] = S[j];
      j--;
    }
    S[j + dt] = key;
    dt--;
  }
  return;
}

template <typename T, typename Compare>
void insert_all_but_second_block(size_t l, size_t r, size_t block_size, T *S,
                                 Compare comp) {
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
      while (j >= (long)l && comp(key, S[j])) {
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

template <typename T, typename Compare>
void __run_bis_iterative(T *S, size_t n, int p, int k, int min_segments,
                         Compare comp) {
  T *E = new T[n];
  omp_set_num_threads(p);

#pragma omp parallel for
  for (size_t l = 0; l < n; l += k) {
    size_t r = std::min(l + k, n);
    bis::__internal::insertion_sort(S, l, r, comp);
  }
  // current_block_size (k'): k, k^2, k^3... k^t-1 ?
  for (size_t current_block_size = k; current_block_size < n;
       current_block_size *= k) {
    size_t segment_size = current_block_size * k;
    int segment_count = (n + segment_size - 1) / segment_size;
    if (segment_count < min_segments) {
      // #if defined(_USETBB)
      multiway_merge_output(S, E, 0, n, current_block_size, comp);
      delete[] E;
      return;
      // #endif // _USETBB
    }
#pragma omp parallel for
    for (int segment_idx = 0; segment_idx < segment_count; segment_idx++) {
      size_t starting_point = segment_idx * segment_size;
      size_t ending_point = starting_point + segment_size;
      if (ending_point > n) {
        ending_point = n;
      }
      // starting point;
      size_t E_offset = segment_idx * segment_size;
      T *my_E = E + E_offset;
      size_t sorted_block_start_idx = starting_point + current_block_size;
      size_t sorted_block_end_idx = sorted_block_start_idx + current_block_size;
      size_t copied_block_size = current_block_size;
      if (sorted_block_end_idx > ending_point) {
        sorted_block_end_idx = ending_point;
        copied_block_size = sorted_block_end_idx - sorted_block_start_idx;
      }

      if (sorted_block_end_idx >= sorted_block_start_idx) {
        std::memcpy(my_E, S + sorted_block_start_idx,
                    copied_block_size * sizeof(T));
      }
      insert_all_but_second_block(starting_point, ending_point,
                                  current_block_size, S, comp);
      insert_second_block(starting_point, ending_point, copied_block_size, S,
                          my_E, comp);
    }
  }
  delete[] E;
}
} // namespace

template <typename T, typename Compare>
void block_insertion_sort(T *S, size_t n, Compare comp,
                          size_t k = PAR_DEFAULT_BLOCK_SIZE,
                          int ss = MWM_SWITCH_STRAT) {
#if !defined(_OPENMP)
  return bis::sequential_sorter::block_insertion_sort(S, n, comp);
#endif
  omp_set_max_active_levels(3);
  int p = omp_get_max_threads();
  if (p == 1) {
    return bis::sequential_sorter::block_insertion_sort(S, n, comp);
  }
  int min_segments = 0;
  switch (ss) {
  case 0: {
    min_segments = p / 2;
    break;
  };
  case 1: {
    min_segments = p;
    break;
  };
  default: {
    break;
  };
  }
  min_segments = std::max(1, min_segments);
  __run_bis_iterative(S, n, p, k, min_segments, comp);
}
} // namespace parallel_sorter

} // namespace bis

#endif /* BASIC_DRF_H_ */
