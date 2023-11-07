#include "experimental_pbis.hpp"
#include "bis.hpp"
#include "common/utils.hpp"
#include "helpers.hpp"
#include <cstddef>
#include <cstdint>
#include <tbb/global_control.h>

// #include "pstl/"
#include "execution"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <omp.h>
// #include <parallel/algorithm>
// #include <parallel/multiway_merge.h>
#include <stdio.h>
#include <utility>

namespace bis {

namespace {
// r debe ser n -1, E debe venir con el offset aplicado
template <typename T>
void block_insertion_sort(T *S, T *E, size_t l, size_t r, int p, int k) {
  size_t n = r - l + 1;
  if (n <= (size_t)k) {
    return insertion_sort(S, l, r + 1);
  }
  omp_set_num_threads(p);
#pragma omp parallel for
  for (size_t i = l; i < r; i += k) {
    size_t end = std::min(i + k, r + 1);
    insertion_sort(S, i, end);
  }
  size_t current_block_size = k;

  size_t dt = std::log(n) / std::log(k);
  dt = std::pow(k, dt);
  size_t potBIS = log2(k);

  // si la entrada es de tamanio potencia de k, el espacio adicional necesario
  // = n/k
  if (dt == n) {
    dt >>= potBIS;
  }
  T key;
  // report_num_threads(3);
  // current_block_size (k'): k, k^2, k^3... k^t-1 ?
  while (current_block_size < n) {
    int max_threads = std::min(n / (current_block_size * k), (size_t)p);
    omp_set_num_threads(max_threads);
#pragma omp parallel for
    for (size_t starting_point = l; starting_point < r + 1;
         starting_point += k * current_block_size) {
      int ID = omp_get_thread_num();
      int thread_offset = ID * current_block_size;
      size_t ending_point = starting_point + (k * current_block_size) - 1;

      if (ending_point > r) {
        ending_point = r;
      }

      size_t sorted_block_start_idx = starting_point + current_block_size;
      size_t sorted_block_end_idx = sorted_block_start_idx + current_block_size;
      size_t copied_block_size = current_block_size;
      if (sorted_block_end_idx > ending_point + 1) {
        sorted_block_end_idx = ending_point + 1;
        copied_block_size = sorted_block_end_idx - sorted_block_start_idx;
      }

      if (sorted_block_end_idx >= sorted_block_start_idx) {
        std::memcpy(E + thread_offset, S + sorted_block_start_idx,
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
  return;
}

template <typename T>
void block_insertion_sort_seq(T *S, T *E, size_t l, size_t r, int k) {
  size_t n = r - l + 1;
  if (n <= (size_t)k) {
    return insertion_sort(S, l, r + 1);
  }
  for (size_t i = l; i < r; i += k) {
    size_t end = std::min(i + k, r + 1);
    insertion_sort(S, i, end);
  }
  size_t current_block_size = k;
  size_t dt = std::log(n) / std::log(k);
  dt = std::pow(k, dt);
  uint8_t potBIS = log2(k);

  // si la entrada es de tamanio potencia de k, el espacio adicional necesario
  // = n/k
  if (dt == n) {
    dt >>= potBIS;
  }
  T key;
  // current_block_size (k'): k, k^2, k^3... k^t-1 ?
  while (current_block_size < n) {
    for (size_t starting_point = l; starting_point < r + 1;
         starting_point += k * current_block_size) {
      size_t ending_point = starting_point + (k * current_block_size) - 1;
      if (ending_point > r) {
        ending_point = r;
      }

      size_t sorted_block_start_idx = starting_point + current_block_size;
      size_t sorted_block_end_idx = sorted_block_start_idx + current_block_size;
      size_t copied_block_size = current_block_size;
      if (sorted_block_end_idx > ending_point + 1) {
        sorted_block_end_idx = ending_point + 1;
        copied_block_size = sorted_block_end_idx - sorted_block_start_idx;
      }

      if (sorted_block_end_idx >= sorted_block_start_idx) {
        std::memcpy(E, S + sorted_block_start_idx,
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
        key = E[i];
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
  return;
}

template <typename T>
void ptask_bis_debug(T *S, size_t l, size_t r, int &k, int p,
                     experimental::MergeStrat merge_strat, size_t dt) {
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
  size_t e_start = 0;
  T *E = new T[dt];

  omp_set_num_threads(p);
  auto pstart = omp_get_wtime();
  // omp_set_nested(1);
#pragma omp parallel
#pragma omp single
#pragma omp taskloop num_tasks(k) untied
  for (size_t i = 0; i < (size_t)k; i++) {
    size_t il = l + i * current_block_size;
    size_t e_idx = e_start + i * next_block_size;
    size_t ir = std::min(il + current_block_size, r);
    // block_insertion_sort(S, E + e_idx, il, ir - 1, 1, k);
    block_insertion_sort_seq(S, E + e_idx, il, ir - 1, k);
  }
  auto pfinish = omp_get_wtime();
  DEBUG((pfinish - pstart));

  auto ls_ini = std::chrono::high_resolution_clock::now();
  // auto start = std::chrono::high_resolution_clock::now();
  if (merge_strat == experimental::MergeStrat::bis) {
    size_t second_block_start = l + current_block_size;
    size_t second_block_end = second_block_start + current_block_size;
    size_t copied_block_size = current_block_size;
    if (second_block_end > r) {
      second_block_end = r;
      copied_block_size = second_block_end - second_block_start;
    }

    if (second_block_end >= second_block_start) {
      memcpy(E + e_start, S + second_block_start,
             copied_block_size * sizeof(T));
      // std::copy(std::execution::par, S + second_block_start,
      //           S + second_block_end, E + e_start);
      // std::copy(S + second_block_start, S + second_block_end, E + e_start);
    }
    auto start = std::chrono::high_resolution_clock::now();
    bis::internal::insert_all_but_second_block(l, r, current_block_size, S);
    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    DEBUG(elapsed.count());
    start = std::chrono::high_resolution_clock::now();
    bis::internal::insert_second_block(l, r, copied_block_size, S, E + e_start);
    finish = std::chrono::high_resolution_clock::now();
    elapsed = finish - start;
    DEBUG(elapsed.count());
  } else {
    bis::helpers::mergeAdjacentSequences(S, l, r, current_block_size,
                                         merge_strat, p);
  }
  auto ls_end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> ls_elapsed = ls_end - ls_ini;
  DEBUG(ls_elapsed.count());
  delete[] E;
}

template <typename T>
void ptask_rec_bis_test3(T *S, size_t l, size_t r, int p, int &k,
                         experimental::MergeStrat merge_strat, size_t dt) {
  if (l >= r) {
    return;
  }
  auto tstart = omp_get_wtime();
  size_t n = r - l;
  auto strat = experimental::AllNames[merge_strat];
  DEBUG(n);
  // DEBUG(strat);
  if (n <= (size_t)k) {
    return insertion_sort(S, l, r);
  }
  size_t kpow = ceil(log(n) / log(k));
  size_t current_block_size = pow(k, kpow - 1);
  size_t next_block_size = pow(k, kpow - 2);
  size_t e_start = 0;
  T *E = new T[dt];

  // bool insertion_base = merge_strat ==

  // experimental::MergeStrat::inplace_mixed; report_num_threads(2);
  auto pstart = omp_get_wtime();

  // omp_set_nested(1);
#pragma omp parallel
#pragma omp single
#pragma omp taskloop num_tasks(p) untied
  for (size_t i = 0; i < (size_t)k; i++) {
    size_t il = l + i * current_block_size;
    size_t e_idx = e_start + i * next_block_size;
    size_t ir = std::min(il + current_block_size, r);
    // block_insertion_sort(S, E + e_idx, il, ir - 1, 1, k);
    block_insertion_sort_seq(S, E + e_idx, il, ir - 1, k);
  }

  auto pfinish = omp_get_wtime();
  DEBUG((pfinish - pstart));

  if (merge_strat == experimental::MergeStrat::bis) {
    auto start = std::chrono::high_resolution_clock::now();
    size_t second_block_start = l + current_block_size;
    size_t second_block_end = second_block_start + current_block_size;
    size_t copied_block_size = current_block_size;
    if (second_block_end > r) {
      second_block_end = r;
      copied_block_size = second_block_end - second_block_start;
    }

    if (second_block_end >= second_block_start) {
      // auto start = std::chrono::high_resolution_clock::now();
      memcpy(E + e_start, S + second_block_start,
             copied_block_size * sizeof(T));
      // std::copy(std::execution::par, S + second_block_start,
      //           S + second_block_end, E + e_start);
      // std::copy(S + second_block_start, S + second_block_end, E + e_start);
      // auto finish = std::chrono::high_resolution_clock::now();
      // std::chrono::duration<double> elapsed = finish - start;
      // DEBUG(elapsed.count());
    }
    bis::internal::insert_all_but_second_block(l, r, current_block_size, S);
    // start = std::chrono::high_resolution_clock::now();
    bis::internal::insert_second_block(l, r, copied_block_size, S, E + e_start);
    // finish = std::chrono::high_resolution_clock::now();
    auto finish = std::chrono::high_resolution_clock::now();
    // elapsed = finish - start;
    std::chrono::duration<double> elapsed = finish - start;
    DEBUG(elapsed.count());
  } else {
    auto start = std::chrono::high_resolution_clock::now();
    bis::helpers::mergeAdjacentSequences(S, l, r, current_block_size,
                                         merge_strat, p);
    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    DEBUG(elapsed.count());
  }
  delete[] E;
  auto tfinish = omp_get_wtime();
  DEBUG((tfinish - tstart));
  // auto finish = std::chrono::high_resolution_clock::now();
  // std::chrono::duration<double> elapsed = finish - start;
  // DEBUG(elapsed.count));
  // DEBUG(elapsed.count());
  // utils::print_array(S, n);
  // DEBUG2("merge done!");
  // DEBUG2("kek");
}

template <typename T>
void ptask_rec_bis_test2(T *S, size_t l, size_t r, int &k, int p,
                         experimental::MergeStrat merge_strat, size_t dt) {
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
  size_t e_start = 0;
  T *E = new T[dt];

  if (p > 1) {
    omp_set_num_threads(p);
#pragma omp parallel
    {
#pragma omp single
#pragma omp taskloop num_tasks(p) untied
      for (size_t i = 0; i < (size_t)k; i++) {
        size_t il = l + i * current_block_size;
        size_t e_idx = e_start + i * next_block_size;
        size_t ir = std::min(il + current_block_size, r);
        block_insertion_sort_seq(S, E + e_idx, il, ir - 1, k);
      }
    }
  } else {
    for (size_t i = 0; i < (size_t)k; i++) {
      size_t il = l + i * current_block_size;
      size_t e_idx = e_start + i * next_block_size;
      size_t ir = std::min(il + current_block_size, r);
      block_insertion_sort_seq(S, E + e_idx, il, ir - 1, k);
    }
  }

  // auto start = std::chrono::high_resolution_clock::now();
  if (merge_strat == experimental::MergeStrat::bis) {
    size_t second_block_start = l + current_block_size;
    size_t second_block_end = second_block_start + current_block_size;
    size_t copied_block_size = current_block_size;
    if (second_block_end > r) {
      second_block_end = r;
      copied_block_size = second_block_end - second_block_start;
    }

    if (second_block_end >= second_block_start) {
      memcpy(E + e_start, S + second_block_start,
             copied_block_size * sizeof(T));
    }
    bis::internal::insert_all_but_second_block(l, r, current_block_size, S);
    bis::internal::insert_second_block(l, r, copied_block_size, S, E + e_start);
  } else {
    bis::helpers::mergeAdjacentSequences(S, l, r, current_block_size,
                                         merge_strat, p);
  }
  delete[] E;
}

} // namespace
namespace experimental {

bool isPowerOfTwo(int num) { return num > 0 && (num & (num - 1)) == 0; }

template <typename T>
void rec_bis_general(T *S, size_t l, size_t r, int p, int k,
                     MergeStrat merge_strat) {
  // find the greates power of K that is <= n
  omp_set_num_threads(p);
  omp_set_nested(1);
  size_t n = r - l;
  size_t dt = std::log(n) / std::log(k);
  dt = std::pow(k, dt);
  if (dt == n) {
    return ptask_bis_debug(S, l, r, p, k, merge_strat, dt);
  }
  size_t n_kproblems = n / dt;
  size_t orphan_size = n - (dt * n_kproblems);
  size_t l_orphan = l + (dt * n_kproblems);
  int ppertaks = 1;
  // DEBUG(orphan_size);
  // DEBUG(n_kproblems);
  if (n_kproblems == 1) {
#pragma omp parallel
#pragma omp single
#pragma omp taskgroup
    {
      // #pragma omp task untied priority(dt / (orphan_size + (orphan_size ==
      // 0)))
#pragma omp task
      ptask_bis_debug(S, l, l + dt, k, p, merge_strat, dt);
#pragma omp task
      {
        if (orphan_size > 0) {
          T *E = new T[dt];
          block_insertion_sort(S, E, l_orphan, r - 1, 1, k);
          delete[] E;
        }
      }
    }
  } else {
#pragma omp parallel
#pragma omp single
    {
#pragma omp taskgroup
      {
#pragma omp taskloop untied
        for (size_t i = 0; i < n_kproblems; i++) {
          ptask_bis_debug(S, l + (i * dt), l + (i * dt) + dt, k, ppertaks,
                          merge_strat, dt);
        }
#pragma omp task untied
        if (orphan_size > 0) {
          T *E = new T[dt];
          block_insertion_sort(S, E, l_orphan, r - 1, ppertaks, k);
          delete[] E;
        }
      }
    }
  }

  experimental::MergeStrat default_strat = merge_strat;
  // n_kproblems >= 4 ? experimental::MergeStrat::inplace_pfor
  // : experimental::MergeStrat::inplace_exec_policy;
  if (orphan_size == 0 || (n_kproblems > 1 && isPowerOfTwo(n_kproblems + 1))) {
    // DEBUG2("merge nk directly");
    bis::helpers::mergeAdjacentSequences(S, l, r, dt, default_strat, p);
  } else if (n_kproblems == 1) {
    // DEBUG2("merge dt and orphan");
    std::inplace_merge(std::execution::par, S + l, S + l_orphan, S + r);
  } else {
    // DEBUG2("2 stages");
    bis::helpers::mergeAdjacentSequences(S, l, l_orphan, dt, default_strat, p);
    // DEBUG2("merge adjacent done");
    std::inplace_merge(std::execution::par, S + l, S + l_orphan, S + r);
    // DEBUG2("merge nk and orphan done");
  }
}

template <typename T>
void rec_bis_general2(T *S, size_t l, size_t r, int p, int k,
                      MergeStrat merge_strat) {
  // find the greates power of K that is <= n
  omp_set_num_threads(p);
  omp_set_nested(1);
  size_t n = r - l;
  size_t dt = std::log(n) / std::log(k);
  dt = std::pow(k, dt);
  if (dt == n) {
    return ptask_bis_debug(S, l, r, p, k, merge_strat, dt);
  }
  int n_kproblems = n / dt;
  size_t orphan_size = n - (dt * n_kproblems);
  size_t l_orphan = l + (dt * n_kproblems);

  bool nested = true;
  int n_threads = p;
  int n_tasks = n_kproblems;
  if (n_kproblems >= p) {
    nested = false;
    n_threads = 1;
  } else if (n_kproblems > 1) {
    int pspow = std::log(n_kproblems) / std::log(2);
    // DEBUG(pspow);
    n_tasks = std::pow(2, pspow);
    n_threads = std::max((p / n_tasks), 1);
  }

  int orphan_threads = n_threads * (orphan_size > 0);
  if (orphan_size > 0) {
    float orphan_ratio = (float)orphan_size / dt;
    if (orphan_ratio < 0.8) {
      size_t bdt = dt;
      while (orphan_threads > 1 && orphan_size < bdt) {
        orphan_threads /= 2;
        bdt /= 2;
      }
    }
  }

  // DEBUG(p);
  // DEBUG(n_kproblems);
  // DEBUG(n_threads);
  // DEBUG(n_tasks);
  // DEBUG(dt);
  // DEBUG(orphan_size);
  // DEBUG(orphan_threads);
  omp_set_nested(nested);
  if (n_kproblems == 1) {

#pragma omp parallel
#pragma omp single
#pragma omp taskgroup
    {
#pragma omp task untied
      ptask_rec_bis_test2(S, l, l + dt, k, p, merge_strat, dt);
#pragma omp task untied
      {
        if (orphan_size > 0) {
          T *E = new T[dt];
          block_insertion_sort(S, E, l_orphan, r - 1, orphan_threads, k);
          delete[] E;
        }
      }
    }
  } else {
#pragma omp parallel
#pragma omp single
    {
#pragma omp taskgroup
      {
#pragma omp taskloop untied num_tasks(n_tasks)
        for (size_t i = 0; i < (size_t)n_kproblems; i++) {
          ptask_rec_bis_test2(S, l + (i * dt), l + (i * dt) + dt, k, n_threads,
                              merge_strat, dt);
        }
#pragma omp task untied
        if (orphan_size > 0) {
          T *E = new T[dt];
          block_insertion_sort(S, E, l_orphan, r - 1, orphan_threads, k);
          delete[] E;
        }
      }
    }
  }
  experimental::MergeStrat default_strat = merge_strat;
  if (orphan_size == 0 || (n_kproblems > 1 && isPowerOfTwo(n_kproblems + 1))) {
    bis::helpers::mergeAdjacentSequences(S, l, r, dt, default_strat, p);
  } else if (n_kproblems == 1) {
    std::inplace_merge(std::execution::par, S + l, S + l_orphan, S + r);
  } else {
    bis::helpers::mergeAdjacentSequences(S, l, l_orphan, dt, default_strat, p);
    std::inplace_merge(std::execution::par, S + l, S + l_orphan, S + r);
  }
}

template <typename T>
void __run_bis_iterative(T *S, size_t n, int p, int k, int min_segments) {
  T *E = new T[n];
  if (p == 1) {
    block_insertion_sort_seq(S, E, 0, n - 1, k);
    delete[] E;
    return;
  }
  omp_set_num_threads(p);

#pragma omp parallel for
  for (size_t l = 0; l < n; l += k) {
    size_t r = std::min(l + k, n);
    bis::insertion_sort(S, l, r);
  }
  // current_block_size (k'): k, k^2, k^3... k^t-1 ?
  for (size_t current_block_size = k; current_block_size < n;
       current_block_size *= k) {
    size_t segment_size = current_block_size * k;
    int segment_count = (n + segment_size - 1) / segment_size;
    if (segment_count < min_segments) {
      // printf("n: %lu, segment_count: %d, k: %d, min_segments: %d\n", n,
      //        segment_count, k, min_segments);
      helpers::multiway_merge_output(S, E, 0, n, current_block_size, p);
      delete[] E;
      return;
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
      bis::internal::insert_all_but_second_block(starting_point, ending_point,
                                                 current_block_size, S);
      bis::internal::insert_second_block(starting_point, ending_point,
                                         copied_block_size, S, my_E);
    }
  }
  delete[] E;
  return;
}

template <typename T>
void pbis_task_split(T *S, size_t n, int p, int k, int method) {
  if (p < 2 || n <= bis::SEQUENTIAL_THRESHOLD) {
    return bis::block_insertion_sort(S, n, k);
  }
  int split_factor = std::max(std::min(method, 2), 1);
  int pp = p * split_factor;

  size_t np = std::ceil(n / pp);
  T *E = new T[n];
  omp_set_num_threads(p);
  omp_set_max_active_levels(2);
#pragma omp parallel for
  for (size_t i = 0; i < (size_t)pp; i++) {
    size_t il = i * np;
    size_t ir = std::min(np, n - il);
    block_insertion_sort_seq(S + il, E + il, 0, ir - 1, k);
  }
  helpers::multiway_merge_output(S, E, 0, n, np, pp);
  delete[] E;
}

template <typename T>
void rec_bis_general3(T *S, size_t l, size_t r, int p, int k) {
  omp_set_max_active_levels(3);
  int method = 1;
  size_t n = r - l;
  size_t dt = std::log(n) / std::log(k);
  dt = std::pow(k, dt);
  if (dt == n) {
    return __run_bis_iterative(S + l, n, p, k, method);
  }
  int n_kproblems = n / dt;
  size_t orphan_size = n - (dt * n_kproblems);
  size_t l_orphan = l + (dt * n_kproblems);

  int n_threads = p;
  int n_tasks = n_kproblems;
  if (n_kproblems >= p) {
    n_threads = 1;
  } else if (n_kproblems > 1) {
    int pspow = std::log(n_kproblems) / std::log(2);
    n_tasks = std::pow(2, pspow);
    n_threads = std::max((p / n_tasks), 1);
  }

  int orphan_threads = n_threads * (orphan_size > 0);
  if (orphan_size > 0) {
    float orphan_ratio = (float)orphan_size / dt;
    if (orphan_ratio < 0.8) {
      size_t bdt = dt;
      while (orphan_threads > 1 && orphan_size < bdt) {
        orphan_threads /= 2;
        bdt /= 2;
      }
    }
  }

  DEBUG(p);
  DEBUG(n_kproblems);
  DEBUG(n_threads);
  DEBUG(n_tasks);
  DEBUG(dt);
  DEBUG(l_orphan);
  DEBUG(orphan_size);
  DEBUG(orphan_threads);
  // omp_set_nested(nested);
  if (n_kproblems == 1) {

#pragma omp parallel
#pragma omp single
#pragma omp taskgroup
    {
#pragma omp task untied
      // ptask_rec_bis_test2(S, l, l + dt, k, p, merge_strat, dt);
      __run_bis_iterative(S + l, dt, p, k, method);
      std::cout << "[ " << l << ", " << (l + dt) << "]  issorted? "
                << std::is_sorted(S + l, S + l + dt) << std::endl;
#pragma omp task untied
      {
        if (orphan_size > 0) {
          // T *E = new T[dt];
          __run_bis_iterative(S + l_orphan, orphan_size, orphan_threads, k,
                              method);
          std::cout << "[ " << l_orphan << ", " << (l_orphan + orphan_size)
                    << "]  issorted? "
                    << std::is_sorted(S + l_orphan, S + l_orphan + orphan_size)
                    << std::endl;
          // block_insertion_sort(S, E, l_orphan, r - 1, orphan_threads, k);
          // delete[] E;
        }
      }
    }
  } else {

    DEBUG2("more than 1 kproble");

#pragma omp parallel
#pragma omp single
    {
#pragma omp taskgroup
      {
#pragma omp taskloop untied num_tasks(n_tasks)
        for (size_t i = 0; i < (size_t)n_kproblems; i++) {
          int lt = l + i * dt;
          __run_bis_iterative(S + lt, dt, n_threads, k, method);
          std::cout << "[ " << lt << ", " << (lt + dt) << "]  issorted? "
                    << std::is_sorted(S + lt, S + lt + dt) << std::endl;
          // ptask_rec_bis_test2(S, l + (i * dt), l + (i * dt) + dt, k,
          // n_threads,
          //                     merge_strat, dt);
        }
#pragma omp task untied
        if (orphan_size > 0) {
          // T *E = new T[dt];
          __run_bis_iterative(S + l_orphan, orphan_size, orphan_threads, k,
                              method);
          std::cout << "lt: " << l_orphan << "  issorted? "
                    << std::is_sorted(S + l_orphan, S + r) << std::endl;
          // block_insertion_sort(S, E, l_orphan, r - 1, orphan_threads, k);
          // delete[] E;
        }
      }
    }
  }
  bis::helpers::multiway_merge(S + l, 0, n, dt, p);
}

template <typename T>
void rec_block_insertion_sort(T *S, size_t l, size_t r, int p, int k,
                              MergeStrat merge_strat) {
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
  // if (p > 1) {
  // omp_set_num_threads(4);
  ptask_bis_debug(S, l, r, k, p, merge_strat, dt);
  // } else {
  // internal::rec_bis_clean(S, E, l, r, 0, k);
  // }
}

int get_num_segments(size_t n, int k, int mwm_threshold) {
  size_t dt = std::log(n) / std::log(k);
  size_t block_size = std::pow(k, dt);
  int block_count = (int)std::ceil(n / block_size * k);
  while (block_count >= mwm_threshold) {
    block_count = std::ceil(block_count / k);
  }
  return block_count;
}

std::pair<int, int> get_last_level_num_segments(size_t n, int p) {
  if (get_num_segments(n, 16, p / 2) < get_num_segments(n, 8, p)) {
    return std::make_pair(16, p / 2);
  }
  return std::make_pair(8, p);
}

template <typename T>
void bis_iterative(T *S, size_t l, size_t r, int p, int k, int method) {
  omp_set_max_active_levels(3);
  size_t n = r - l;
  if (n <= bis::SEQUENTIAL_THRESHOLD) {
    return bis::block_insertion_sort(S, n, k);
  }
  // size_t dt = std::log(n) / std::log(k);
  // dt = std::pow(k, dt);
  int _k = k, min_segments = 0;
  switch (method) {
  case 1: {
    min_segments = p;
    break;
  };
  case 2: {
    min_segments = p / 2;
    break;
  };
  case 3: {
    auto [__k, _min_segments] = get_last_level_num_segments(n, p);
    _k = __k;
    min_segments = _min_segments;
    break;
  };
  default: {
    break;
  };
  }
  min_segments = std::max(1, min_segments);
  // DEBUG(_k);
  // DEBUG(min_segments);
  // size_t potBIS = log2(k);
  // por mientras usaremos n mem adicional
  // if (dt == n) {
  //   dt >>= potBIS;
  // }
  // T *E = new T[dt];
  // //
  __run_bis_iterative(S + l, n, p, _k, min_segments);
}

template void rec_block_insertion_sort<int>(int *S, size_t l, size_t r, int p,
                                            int k, MergeStrat merge_strat);
template void rec_block_insertion_sort<long>(long *S, size_t l, size_t r, int p,
                                             int k, MergeStrat merge_strat);
template void rec_block_insertion_sort<float>(float *S, size_t l, size_t r,
                                              int p, int k,
                                              MergeStrat merge_strat);
template void rec_block_insertion_sort<double>(double *S, size_t l, size_t r,
                                               int p, int k,
                                               MergeStrat merge_strat);

template void rec_bis_general<int>(int *S, size_t l, size_t r, int p, int k,
                                   MergeStrat merge_strat);
template void rec_bis_general<long>(long *S, size_t l, size_t r, int p, int k,
                                    MergeStrat merge_strat);
template void rec_bis_general<double>(double *S, size_t l, size_t r, int p,
                                      int k, MergeStrat merge_strat);
template void rec_bis_general<float>(float *S, size_t l, size_t r, int p, int k,
                                     MergeStrat merge_strat);

template void rec_bis_general2<int>(int *S, size_t l, size_t r, int p, int k,
                                    MergeStrat merge_strat);

template void rec_bis_general2<float>(float *S, size_t l, size_t r, int p,
                                      int k, MergeStrat merge_strat);
template void rec_bis_general2<long>(long *S, size_t l, size_t r, int p, int k,
                                     MergeStrat merge_strat);
template void rec_bis_general2<double>(double *S, size_t l, size_t r, int p,
                                       int k, MergeStrat merge_strat);
template void rec_bis_general2<size_t>(size_t *S, size_t l, size_t r, int p,
                                       int k, MergeStrat merge_strat);

// template void rec_bis_general3<int>(int *S, size_t l, size_t r, int p, int k,
//                                     experimental::MergeStrat merge_strat);
template void rec_bis_general3<int>(int *S, size_t l, size_t r, int p, int k);
template void rec_bis_general3<float>(float *S, size_t l, size_t r, int p,
                                      int k);
template void rec_bis_general3<double>(double *S, size_t l, size_t r, int p,
                                       int k);
template void rec_bis_general3<long>(long *S, size_t l, size_t r, int p, int k);
template void rec_bis_general3<size_t>(size_t *S, size_t l, size_t r, int p,
                                       int k);

template void bis_iterative<int>(int *S, size_t l, size_t r, int p, int k,
                                 int method);
template void bis_iterative<float>(float *S, size_t l, size_t r, int p, int k,
                                   int method);
template void bis_iterative<double>(double *S, size_t l, size_t r, int p, int k,
                                    int method);
template void bis_iterative<long>(long *S, size_t l, size_t r, int p, int k,
                                  int method);
template void bis_iterative<size_t>(size_t *S, size_t l, size_t r, int p, int k,
                                    int method);

template void pbis_task_split(int *S, size_t n, int p, int k, int method);
template void pbis_task_split(float *S, size_t n, int p, int k, int method);
template void pbis_task_split(long *S, size_t n, int p, int k, int method);
template void pbis_task_split(double *S, size_t n, int p, int k, int method);
template void pbis_task_split(size_t *S, size_t n, int p, int k, int method);

} // namespace experimental
} // namespace bis
