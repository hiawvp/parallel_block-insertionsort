#include "pbis.hpp"
#include "bis.hpp"
#include "common/utils.hpp"
#include <algorithm>
#include <chrono>
#include <cstddef>
#include <cstring>
#include <iostream>
#include <math.h>
#include <omp.h>

namespace bis {

namespace internal {

template <typename T>
void blockInsertionSort(T *A, T *E, size_t l, size_t r, size_t e_start, int k) {
  size_t pos, ini, end, dt, bLev, sp, ep, lenB, n = r - l + 1, potBIS = log2(k);
  long int i, j;
  T key;

  // std::cout << "potBIS = " << potBIS << std::endl;

  // 1- sort Blk_1, Blk_2, ..
  ini = l;
  end = ini + k;
  while (ini <= r) {
    if (end > r)
      end = r + 1;
    for (pos = ini + 1; pos < end; pos++) {
      key = A[pos];
      j = pos - 1;
      while (j >= (long int)ini && key < A[j]) {
        A[j + 1] = A[j];
        j--;
      }
      A[j + 1] = key;
    }
    ini = end;
    end += k;
  };
  if (n <= (size_t)k)
    return;

  dt = log(n) / log(k);
  dt = pow(k, dt);
  if (dt == n)
    dt >>= potBIS;
  T *ArrBI = E + e_start;
  // std::cout << "len ArrBI = " << dt << std::endl;

  bLev = k;
  while (bLev < n) {
    sp = l;
    while (sp <= r) {
      ep = sp + (bLev << potBIS) - 1;
      if (ep > r)
        ep = r;
      ini = sp + bLev;
      end = ini + bLev;
      if (end > ep + 1) {
        end = ep + 1;
        lenB = end - ini;
      } else
        lenB = bLev;

      // 2- copy Blk_1 of A[sp..ep] in ArrBI
      for (pos = ini; pos < end; pos++)
        ArrBI[pos - ini] = A[pos];

      // 3- insert all sorted blocks from blk_2 to ep in of A[sp..ep]
      ini = sp + (bLev << 1);
      end = ini + bLev;
      j = sp + bLev - 1;
      while (ini <= ep) {
        if (end > ep + 1) {
          end = ep + 1;
          dt = end - ini;
        } else
          dt = bLev;
        for (i = end - 1; i >= (long int)ini; i--) {
          key = A[i];
          while (j >= (long int)sp && key < A[j]) {
            A[j + dt] = A[j];
            j--;
          }
          A[j + dt] = key;
          dt--;
        }
        j = ini - 1;
        ini = end;
        end += bLev;
      }

      // 4- Insert ArrBI (blk_2) in A
      dt = lenB;
      j = ep - lenB;
      for (i = lenB - 1; i >= 0; i--) {
        key = ArrBI[i];
        while (j >= (long int)sp && key < A[j]) {
          A[j + dt] = A[j];
          j--;
        }
        A[j + dt] = key;
        dt--;
      }
      // A[sp..ep] is sorted
      sp = ep + 1;
    }
    // to multiply by k for the nest level
    bLev <<= potBIS;
  }
}

template <typename T>
void ptask_rec_bis(T *S, T *E, size_t l, size_t r, size_t e_start, int &k) {
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
  // after this, we should have kpow-1 sorted blocks of size current_block_size
  // TODO: handle task or thread limit, see if any var need to be explicit
  // private or shared
  // #pragma omp task firstprivate(i, e_start) if (current_block_size >
  // 16)
#pragma omp parallel
  {
#pragma omp single
    {
      for (size_t i = l, e_idx = e_start; i < r;
           i += current_block_size, e_idx += next_block_size) {
#pragma omp task shared(S, E) firstprivate(i, e_start)
        {
          size_t ir = std::min(i + current_block_size, r);
          ptask_rec_bis(S, E, i, ir, e_idx, k);
          // blockInsertionSort(S, E, i, ir - 1, e_idx, k);
        }
      }
#pragma omp taskwait
    }
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
} // namespace internal

} // namespace internal

namespace parallel {

template <typename T> void block_insertion_sort(T *S, size_t n, int p, int k) {
  if (p == 1) {
    return block_insertion_sort<T>(S, n);
  }
  if (n <= (size_t)k) {
    return insertion_sort(S, 0, n);
  }
  omp_set_num_threads(p);
#pragma omp parallel for
  for (size_t l = 0; l < n; l += k) {
    size_t r = std::min(l + k, n);
    insertion_sort(S, l, r);
  }
  size_t current_block_size = k;

  size_t dt = std::log(n) / std::log(k);
  dt = std::pow(k, dt);
  size_t potBIS = log2(k);

  // si la entrada es de tamanio potencia de k, el espacio adicional necesario =
  // n/k
  if (dt == n) {
    dt >>= potBIS;
  }
  T *E = new T[dt];
  T key;
  // current_block_size (k'): k, k^2, k^3... k^t-1 ?
  while (current_block_size < n) {
    int max_threads = std::min(n / (current_block_size * k), (size_t)p);
    omp_set_num_threads(max_threads);
#pragma omp parallel for
    for (size_t starting_point = 0; starting_point < n;
         starting_point += k * current_block_size) {
      int ID = omp_get_thread_num();
      int thread_offset = ID * current_block_size;
      size_t ending_point = starting_point + (k * current_block_size) - 1;

      if (ending_point >= n) {
        ending_point = n - 1;
      }

      size_t sorted_block_start_idx = starting_point + current_block_size;
      size_t sorted_block_end_idx = sorted_block_start_idx + current_block_size;
      size_t copied_block_size = current_block_size;
      if (sorted_block_end_idx > ending_point + 1) {
        sorted_block_end_idx = ending_point + 1;
        copied_block_size = sorted_block_end_idx - sorted_block_start_idx;
      }

      if (sorted_block_end_idx >= sorted_block_start_idx) {
        memcpy(E + thread_offset, S + sorted_block_start_idx,
               copied_block_size * sizeof(T));
      }

      sorted_block_start_idx = starting_point + 2 * current_block_size;
      sorted_block_end_idx = sorted_block_start_idx + current_block_size;
      long j = starting_point + current_block_size - 1;

      while (sorted_block_start_idx <= ending_point) {
        size_t delta = current_block_size;
        if (sorted_block_end_idx > ending_point + 1) {
          sorted_block_end_idx = ending_point + 1;
          delta = sorted_block_end_idx - sorted_block_start_idx;
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

      size_t delta = copied_block_size;
      j = ending_point - copied_block_size;

      // aqui se copia el bloque guardado en E
      for (int i = copied_block_size - 1; i >= 0; i--) {
        key = E[i + thread_offset];
        while (j >= (long)starting_point && key < S[j]) {
          S[j + delta] = S[j];
          j--;
        }
        S[j + delta] = key;
        delta--;
      }
      // starting_point = ending_point + 1;
    }
    current_block_size *= k;
  }
  delete[] E;
  return;
}

// rec_block_insertion_sort suugma S input array, l, r limits, k = initial block
// size
template <typename T>
void rec_block_insertion_sort(T *S, size_t l, size_t r, int p, int k) {
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
  if (p > 1) {
    omp_set_num_threads(p);
    internal::ptask_rec_bis(S, E, l, r, 0, k);

  } else {
    internal::rec_bis_clean(S, E, l, r, 0, k);
  }
  delete[] E;
}

template void block_insertion_sort<int>(int *S, size_t n, int p, int k);

template void block_insertion_sort<long>(long *S, size_t n, int p, int k);
template void block_insertion_sort<float>(float *S, size_t n, int p, int k);
template void block_insertion_sort<double>(double *S, size_t n, int p, int k);
template void block_insertion_sort<size_t>(size_t *S, size_t n, int p, int k);

template void rec_block_insertion_sort<int>(int *S, size_t l, size_t r, int p,
                                            int k);

template void rec_block_insertion_sort<float>(float *S, size_t l, size_t r,
                                              int p, int k);

template void rec_block_insertion_sort<long>(long *S, size_t l, size_t r, int p,
                                             int k);
template void rec_block_insertion_sort<double>(double *S, size_t l, size_t r,
                                               int p, int k);
template void rec_block_insertion_sort<size_t>(size_t *S, size_t l, size_t r,
                                               int p, int k);

} // namespace parallel
} // namespace bis
