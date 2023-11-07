#include "data_gen.hpp"
#include "utils.hpp"
#include <algorithm>
#include <cstddef>
#include <iostream>
#include <random>
#include <utility>

namespace dataGeneration {

namespace {
template <typename T>
void randomly_swap_elements(T *A, size_t n, int seed, float swap_proba) {
  std::mt19937 gen(seed + 1);
  std::uniform_real_distribution<float> distribution(0, 1);
  size_t half = n / 2;
  for (size_t i = 0; i < half; i++) {
    if (distribution(gen) <= swap_proba) {
      std::swap(A[i], A[i + half]);
    }
  }
}
} // namespace

template <typename T> std::pair<T, T> get_default_limits();
template <> std::pair<int, int> get_default_limits<int>() {
  return std::make_pair(DEFAULT_MIN_INT_VAL, DEFAULT_MAX_INT_VAL);
}
template <> std::pair<float, float> get_default_limits<float>() {
  return std::make_pair(DEFAULT_MIN_FLOAT_VAL, DEFAULT_MAX_FLOAT_VAL);
}
template <> std::pair<long, long> get_default_limits<long>() {
  return std::make_pair(DEFAULT_MIN_LONG_VAL, DEFAULT_MAX_LONG_VAL);
}
template <> std::pair<double, double> get_default_limits<double>() {
  return std::make_pair(DEFAULT_MIN_DOUBLE_VAL, DEFAULT_MAX_DOUBLE_VAL);
}
template <> std::pair<size_t, size_t> get_default_limits<size_t>() {
  return std::make_pair(DEFAULT_MIN_LONG_VAL, DEFAULT_MAX_LONG_VAL);
}

template <typename T>
void fill_uniform_array(std::mt19937 &generator, T *A, size_t n, T min_val,
                        T max_val);
template <>
void fill_uniform_array<int>(std::mt19937 &generator, int *A, size_t n,
                             int min_val, int max_val) {
  std::uniform_int_distribution<int> distribution(min_val, max_val);
  for (size_t i = 0; i < n; i++) {
    A[i] = distribution(generator);
  }
}

template <>
void fill_uniform_array<long>(std::mt19937 &generator, long *A, size_t n,
                              long min_val, long max_val) {
  std::uniform_int_distribution<long> distribution(min_val, max_val);
  for (size_t i = 0; i < n; i++) {
    A[i] = distribution(generator);
  }
}

template <>
void fill_uniform_array<size_t>(std::mt19937 &generator, size_t *A, size_t n,
                                size_t min_val, size_t max_val) {
  std::uniform_int_distribution<size_t> distribution(min_val, max_val);
  for (size_t i = 0; i < n; i++) {
    A[i] = distribution(generator);
  }
}

template <>
void fill_uniform_array<float>(std::mt19937 &generator, float *A, size_t n,
                               float min_val, float max_val) {
  std::uniform_real_distribution<float> distribution(min_val, max_val);
  for (size_t i = 0; i < n; i++) {
    A[i] = distribution(generator);
  }
}

template <>
void fill_uniform_array<double>(std::mt19937 &generator, double *A, size_t n,
                                double min_val, double max_val) {
  std::uniform_real_distribution<double> distribution(min_val, max_val);
  for (size_t i = 0; i < n; i++) {
    A[i] = distribution(generator);
  }
}
// namespace

template <typename T>
void generate_uniform_array(T *A, size_t n, int seed, size_t discard) {
  auto [min_val, max_val] = dataGeneration::get_default_limits<T>();
  generate_uniform_array(A, n, seed, discard, min_val, max_val);
}

// @returns void
// generate uniform distributed float array
template <typename T>
void generate_uniform_array(T *A, size_t n, int seed, size_t discard, T min_val,
                            T max_val) {
  std::mt19937 generator(seed);
  if (verbose) {
    std::printf("Initializing mt19937 generator with seed %d.\n", seed);
  }
  if (discard > 0) {
    generator.discard(discard);
    if (verbose) {
      std::printf("Skipping the first %zu values.\n", discard);
    }
  }
  if (min_val > max_val) {
    std::tie(min_val, max_val) = get_default_limits<T>();
  }
  if (verbose) {
    std::cout << "Generating " << n << " values in the range [" << min_val
              << "," << max_val << "]. " << std::endl;
  }
  fill_uniform_array<T>(generator, A, n, min_val, max_val);
}

template <typename T>
void generate_normal_array(T *A, size_t n, int seed, size_t discard) {
  std::mt19937 generator(seed);
  if (verbose) {
    std::printf("Initializing mt19937 generator with seed %d.\n", seed);
  }
  if (discard > 0) {
    generator.discard(discard);
    if (verbose) {
      std::printf("Skipping the first %zu values.\n", discard);
    }
  }
  // valores usados por el profe
  float mean = n / 2.0;
  float variance = n * 0.1;
  if (verbose) {
    std::printf(
        "Generating %zu INT normally distributed values with mean=%f and "
        "variance=%f.\n",
        n, mean, variance);
  }
  std::normal_distribution<> distribution{mean, variance};
  // std::function<T()> gen = [&distribution, &generator]() -> T {
  //   return distribution(generator);
  // };
  for (size_t i = 0; i < n; i++) {
    // A[i] = gen();
    A[i] = distribution(generator);
  }
}

template <typename T>
void generate_uniform_sorted_array(T *A, size_t n, int seed, size_t discard,
                                   bool reverse) {
  generate_uniform_array(A, n, seed, discard);
  if (reverse) {
    std::sort(A, A + n, std::greater<T>());
  } else {
    std::sort(A, A + n);
  }
  return;
}

// An element of the first half of the array has swap_proba probability
// of being swapped with an element of the second half
template <typename T>
void generate_uniform_almost_sorted_array(T *A, size_t n, int seed,
                                          size_t discard, bool reverse,
                                          float swap_proba) {
  generate_uniform_sorted_array(A, n, seed, discard, reverse);
  randomly_swap_elements(A, n, seed + 1, swap_proba);
}

template <typename T>
void generate_normal_sorted_array(T *A, size_t n, int seed, size_t discard,
                                  bool reverse) {
  generate_normal_array(A, n, seed, discard);
  if (reverse) {
    std::sort(A, A + n, std::greater<T>());
  } else {
    std::sort(A, A + n);
  }
  return;
}

template <typename T>
void generate_normal_almost_sorted_array(T *A, size_t n, int seed,
                                         size_t discard, bool reverse,
                                         float swap_proba) {
  generate_normal_sorted_array(A, n, seed, discard, reverse);
  randomly_swap_elements(A, n, seed + 1, swap_proba);
}

// #####################################
template void generate_uniform_array<int>(int *A, size_t n, int seed,
                                          size_t discard, int min_val,
                                          int max_val);
template void generate_uniform_array<float>(float *A, size_t n, int seed,
                                            size_t discard, float min_val,
                                            float max_val);
template void generate_uniform_array<size_t>(size_t *A, size_t n, int seed,
                                             size_t discard, size_t min_val,
                                             size_t max_val);
template void generate_uniform_array<double>(double *A, size_t n, int seed,
                                             size_t discard, double min_val,
                                             double max_val);

// #####################################
template void generate_uniform_array<int>(int *A, size_t n, int seed,
                                          size_t discard);

template void generate_uniform_array<float>(float *A, size_t n, int seed,
                                            size_t discard);

template void generate_uniform_array<size_t>(size_t *A, size_t n, int seed,
                                             size_t discard);
template void generate_uniform_array<double>(double *A, size_t n, int seed,
                                             size_t discard);

// #####################################
// uniform sorted
template void generate_uniform_sorted_array<int>(int *A, size_t n, int seed,
                                                 size_t discard, bool reverse);
template void generate_uniform_sorted_array<float>(float *A, size_t n, int seed,
                                                   size_t discard,
                                                   bool reverse);
template void generate_uniform_sorted_array<size_t>(size_t *A, size_t n,
                                                    int seed, size_t discard,
                                                    bool reverse);
template void generate_uniform_sorted_array<double>(double *A, size_t n,
                                                    int seed, size_t discard,
                                                    bool reverse);

// #####################################
template void generate_uniform_almost_sorted_array<int>(
    int *A, size_t n, int seed, size_t discard, bool reverse, float swap_proba);

template void generate_uniform_almost_sorted_array<size_t>(size_t *A, size_t n,
                                                           int seed,
                                                           size_t discard,
                                                           bool reverse,
                                                           float swap_proba);

template void generate_uniform_almost_sorted_array<float>(float *A, size_t n,
                                                          int seed,
                                                          size_t discard,
                                                          bool reverse,
                                                          float swap_proba);
template void generate_uniform_almost_sorted_array<double>(double *A, size_t n,
                                                           int seed,
                                                           size_t discard,
                                                           bool reverse,
                                                           float swap_proba);

// #####################################
template void generate_normal_array<int>(int *A, size_t n, int seed,
                                         size_t discard);
template void generate_normal_array<float>(float *A, size_t n, int seed,
                                           size_t discard);
template void generate_normal_array<long>(long *A, size_t n, int seed,
                                          size_t discard);
template void generate_normal_array<double>(double *A, size_t n, int seed,
                                            size_t discard);
template void generate_normal_array<size_t>(size_t *A, size_t n, int seed,
                                            size_t discard);

// normal sorted
template void generate_normal_sorted_array<int>(int *A, size_t n, int seed,
                                                size_t discard, bool reverse);
template void generate_normal_sorted_array<float>(float *A, size_t n, int seed,
                                                  size_t discard, bool reverse);
template void generate_normal_sorted_array<size_t>(size_t *A, size_t n,
                                                   int seed, size_t discard,
                                                   bool reverse);
template void generate_normal_sorted_array<double>(double *A, size_t n,
                                                   int seed, size_t discard,
                                                   bool reverse);

// #####################################
// normal almost sorted
template void generate_normal_almost_sorted_array<int>(int *A, size_t n,
                                                       int seed, size_t discard,
                                                       bool reverse,
                                                       float swap_proba);
template void generate_normal_almost_sorted_array<float>(float *A, size_t n,
                                                         int seed,
                                                         size_t discard,
                                                         bool reverse,
                                                         float swap_proba);
template void generate_normal_almost_sorted_array<double>(double *A, size_t n,
                                                          int seed,
                                                          size_t discard,
                                                          bool reverse,
                                                          float swap_proba);
template void generate_normal_almost_sorted_array<size_t>(size_t *A, size_t n,
                                                          int seed,
                                                          size_t discard,
                                                          bool reverse,
                                                          float swap_proba);
} // namespace dataGeneration
