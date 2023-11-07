#ifndef PBIS_HELPERS
#define PBIS_HELPERS

#include "common/utils.hpp"
#include "execution"
#include "experimental_pbis.hpp"
#include <algorithm>
#include <cstddef>
#include <omp.h>
#include <parallel/algorithm>
#include <tbb/global_control.h>

namespace bis {
namespace helpers {

template <typename T>
void multiway_merge_output(T *arr, T *output, size_t start, size_t end,
                           size_t step, int p) {
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
  __gnu_parallel::multiway_merge(seqs.begin(), seqs.end(), output, n,
                                 std::less<T>());

  // std::copy(std::execution::par, output, output + n, arr + start);
  std::copy(output, output + n, arr + start);
  return;
}

template <typename T>
void multiway_merge(T *arr, size_t start, size_t end, size_t step, int p) {
  size_t n = end - start;
  T *output = new T[n];
  multiway_merge_output(arr, output, start, end, step, p);
  delete[] output;
}

template <typename T>
void mergeAdjacentSequences(T *arr, size_t start, size_t end, size_t step,
                            experimental::MergeStrat strat, int p) {
  if (step >= end - start) {
    return;
  }
  int nblocks = (end - start) / (2 * step);
  if (strat == experimental::gnu_multiway) {
    return multiway_merge(arr, start, end, step, p);
  }
  omp_set_num_threads(p);
  // tbb::global_control c(tbb::global_control::max_allowed_parallelism, p);
  experimental::MergeStrat next_strat = strat;

  switch (strat) {
  case (experimental::inplace_pfor): {
#pragma omp parallel for
    for (size_t i = start; i < end; i += 2 * step) {
      int mid = i + step - 1;
      int rightEnd = std::min(i + 2 * step - 1, end - 1);
      if (rightEnd > mid) {
        std::inplace_merge(arr + i, arr + mid + 1, arr + rightEnd + 1);
      }
    }
    break;
  }
  case (experimental::inplace_exec_policy): {
    for (size_t i = start; i < end; i += 2 * step) {
      int mid = i + step - 1;
      int rightEnd = std::min(i + 2 * step - 1, end - 1);
      if (rightEnd > mid) {
        std::inplace_merge(std::execution::par, arr + i, arr + mid + 1,
                           arr + rightEnd + 1);
      }
    }
    break;
  }
  case (experimental::inplace_mixed): {
#pragma omp parallel for
    for (size_t i = start; i < end; i += 2 * step) {
      int mid = i + step - 1;
      int rightEnd = std::min(i + 2 * step - 1, end - 1);
      if (rightEnd > mid) {
        std::inplace_merge(arr + i, arr + mid + 1, arr + rightEnd + 1);
      }
    }
    if (nblocks <= p) {
      next_strat = experimental::MergeStrat::inplace_exec_policy;
    }
    break;
  }
  default: {
    return;
  }
  }
  mergeAdjacentSequences(arr, start, end, 2 * step, next_strat, p);
}

// codigo de original_bis/sorting.hpp
template <typename T> void merge(T *A, size_t l, size_t m, size_t r) {
  size_t i, j, k;
  size_t n1 = m - l + 1;
  size_t n2 = r - m;
  T *L = new int[n1];
  T *R = new int[n2];

  for (i = 0; i < n1; i++)
    L[i] = A[l + i];
  for (i = 0; i < n2; i++)
    R[i] = A[m + 1 + i];

  i = j = 0;
  k = l;
  while (i < n1 && j < n2) {
    if (L[i] <= R[j]) {
      A[k] = L[i];
      i++;
    } else {
      A[k] = R[j];
      j++;
    }
    k++;
  }

  // Copy the remaining elements of L[], if there are any
  while (i < n1) {
    A[k] = L[i];
    i++;
    k++;
  }

  // Copy the remaining elements of R[], if there are any
  while (j < n2) {
    A[k] = R[j];
    j++;
    k++;
  }

  delete[] L;
  delete[] R;
}

} // namespace helpers
} // namespace bis

#endif // !PBIS_HELPERS
