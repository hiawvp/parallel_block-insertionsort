#ifndef PBIS_UTILS
#define PBIS_UTILS

#include <algorithm>
#include <iomanip>
#include <ios>
#include <iostream>
#include <stdio.h>
#include <string>
#include <string_view>

#define RAND_LIMIT 12377221

#ifndef TIMING_FLAG
#define TIMING_FLAG 1 // set timing mode
#endif

#if TIMING_FLAG
#define timeit(f, m)                                                           \
  do {                                                                         \
    auto start = std::chrono::high_resolution_clock::now();                    \
    f;                                                                         \
    auto finish = std::chrono::high_resolution_clock::now();                   \
    std::chrono::duration<double> elapsed = finish - start;                    \
    std::cout << "TIME  >>> " << std::setw(25) << m                            \
              << std::setiosflags(std::ios::fixed) << std::setprecision(8)     \
              << std::setw(15) << elapsed.count() << std::setw(4) << " [s]"    \
              << std::endl                                                     \
              << std::flush;                                                   \
  } while (0)
#else
#define timeit(f, m)
#endif

template <typename T> void gen_array(T *A, size_t n) {
  for (size_t i = 0; i < n; i++) {
    A[i] = rand() % RAND_LIMIT;
  }
}

template <typename T> void print_array(T *A, size_t n) {
  for (size_t i = 0; i < n - 1; i++) {
    std::cout << A[i] << ", ";
  }
  std::cout << std::endl;
}

template <typename T> int check_results(T *A, T *correct_output, size_t n) {
  if (std::equal(A, A + n, correct_output)) {
    std::cout << "Sort result Ok" << std::endl;
    return EXIT_SUCCESS;
  } else {
    std::cout << "Bad sort unlucky" << std::endl;
    return EXIT_FAILURE;
  }
}

#endif // !PBIS_UTILS
