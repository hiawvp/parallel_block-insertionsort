// https://stackoverflow.com/a/47495419
// https://www.cs.kent.edu/~jbaker/23Workshop/Chesebrough/mergesort/mergesortOMP.cpp
// https://cw.fel.cvut.cz/old/_media/courses/b4m35pag/lab6_slides_advanced_openmp.pdf
#include "omp.h"
#include <algorithm>
#include <cstddef>
#include <cstring>
#include <vector>

namespace collected {

void merge(int *X, int n, int *tmp) {
  int i = 0;
  int j = n / 2;
  int ti = 0;

  while (i < n / 2 && j < n) {
    if (X[i] < X[j]) {
      tmp[ti] = X[i];
      ti++;
      i++;
    } else {
      tmp[ti] = X[j];
      ti++;
      j++;
    }
  }
  while (i < n / 2) { /* finish up lower half */
    tmp[ti] = X[i];
    ti++;
    i++;
  }
  while (j < n) { /* finish up upper half */
    tmp[ti] = X[j];
    ti++;
    j++;
  }
  memcpy(X, tmp, n * sizeof(int));

} // end of merge()

void mergesort(int *X, int n, int *tmp) {
  if (n < 2)
    return;

#pragma omp task firstprivate(X, n, tmp)
  mergesort(X, n / 2, tmp);

#pragma omp task firstprivate(X, n, tmp)
  mergesort(X + (n / 2), n - (n / 2), tmp);

#pragma omp taskwait

  /* merge sorted halves into sorted list */
  merge(X, n, tmp);
}

// span parallel region outside once outside
void mergesort_omp(int *S, size_t n, int p) {
  omp_set_num_threads(p);
  int *tmp = new int[n];
#pragma omp parallel
  {
#pragma omp single
    mergesort(S, n, tmp);
  }
}

template <typename T>
void mergeSortRecursiveNaive(T *A, size_t left, size_t right) {
  if (left < right) {
    if (right - left >= 32) {
      unsigned long mid = (left + right) / 2;
#pragma omp taskgroup
      {
#pragma omp task shared(A)
        mergeSortRecursiveNaive(A, left, mid);
        mergeSortRecursiveNaive(A, mid + 1, right);
      }
      std::inplace_merge(A + left, A + mid + 1, A + right + 1);
    } else {
      std::sort(A + left, A + right + 1);
    }
  }
}

template <typename T> void mergeSortRecursive(T *A, size_t left, size_t right) {
  if (left < right) {
    if (right - left >= 32) {
      unsigned long mid = (left + right) / 2;
#pragma omp taskgroup
      {
#pragma omp task shared(A) untied if (right - left >= (1 << 14))
        mergeSortRecursive(A, left, mid);
#pragma omp task shared(A) untied if (right - left >= (1 << 14))
        mergeSortRecursive(A, mid + 1, right);
#pragma omp taskyield
      }
      std::inplace_merge(A + left, A + mid + 1, A + right + 1);
    } else {
      std::sort(A + left, A + right + 1);
    }
  }
}

template <typename T> void mergeSortImproved(T *A, size_t n, int p) {
  omp_set_num_threads(p);
#pragma omp parallel
#pragma omp single
  mergeSortRecursive(A, 0, n - 1);
}

template <typename T> void mergeSort(T *A, size_t n, int p) {
  omp_set_num_threads(p);
#pragma omp parallel
#pragma omp single
  mergeSortRecursiveNaive(A, 0, n - 1);
}

template void mergeSort<int>(int *A, size_t n, int p);
template void mergeSort<float>(float *A, size_t n, int p);
template void mergeSort<long>(long *A, size_t n, int p);
template void mergeSort<double>(double *A, size_t n, int p);

template void mergeSortImproved<int>(int *A, size_t n, int p);
template void mergeSortImproved<float>(float *A, size_t n, int p);
template void mergeSortImproved<long>(long *A, size_t n, int p);
template void mergeSortImproved<double>(double *A, size_t n, int p);
template void mergeSortImproved<size_t>(size_t *A, size_t n, int p);
} // namespace collected
