#include "experiments.hpp"
#include "Block-InsertionSort/experimental_pbis.hpp"
#include "Block-InsertionSort/helpers.hpp"
#include "common/sequential_impls.hpp"
#include "data_gen.hpp"
#include "execution"
#include "impl_sorting_times.hpp"
#include "omp.h"
#include "utils.hpp"
#include <cstddef>
#include <cstdio>
#include <iostream>
#include <utility>

namespace exps {

namespace {
const bool verbose = true;

} // namespace

template <typename T>
void run_implementations(Experiment<T> &experiment,
                         GeneratorFunc<T> num_generator,
                         std::vector<int> sizes) {
  size_t max_n = sizes.back();
  T *base_array = new T[max_n];
  T *imp_sorted = new T[max_n];
  double total_time;
  int seed = experiment.get_seed();
  size_t reps = experiment.get_reps();
  size_t skip = 0;
  DEBUG(reps);
  for (auto n : sizes) {
    DEBUG(n);
    num_generator(base_array, n, seed, skip);
    std::sort(base_array, base_array + n);
    for (SortingImplementationDetails<T> &impl :
         experiment.get_implementations()) {
      DEBUG(impl.get_name());
      std::vector<double> times;
      for (size_t i = 0; i < reps; i++) {
        num_generator(imp_sorted, n, seed, skip);
        auto start = std::chrono::high_resolution_clock::now();
        impl.run(imp_sorted, n);
        auto finish = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = finish - start;
        // std::chrono::duration<double, std::milli> elapsed = finish - start;
        total_time = elapsed.count();
        if (std::equal(std::execution::par, base_array, base_array + n,
                       imp_sorted)) {
          times.push_back(total_time);
        } else {
          std::cout << "ERROR! >>> " << impl.get_name() << std::endl;
        }
      }
      impl.add_times(n, times);
    }
    skip += n;
  }
  delete[] base_array;
  delete[] imp_sorted;
  return;
}

template <typename T>
std::pair<InputDataInfo, GeneratorFunc<T>>
get_default_distribution(DataDist dist) {
  std::string input_type = utils::get_typename(utils::get_datatype<T>());
  std::string input_distribution;
  GeneratorFunc<T> num_generator;
  std::vector<std::pair<std::string, std::string>> extra_args;

  switch (dist) {
  case dataGeneration::normal: {
    input_distribution = "normal";
    extra_args = {std::make_pair("mean", "n/2"),
                  std::make_pair("std", "n*0.1")};
    num_generator = static_cast<void (*)(T *, size_t, int, size_t)>(
        &dataGeneration::generate_normal_array);
    break;
  }
  case dataGeneration::uniform: {
    input_distribution = "uniform";
    auto [minval, maxval] = dataGeneration::get_default_limits<T>();
    extra_args = {std::make_pair("min_val", std::to_string(minval)),
                  std::make_pair("max_val", std::to_string(maxval))};
    num_generator = static_cast<void (*)(T *, size_t, int, size_t)>(
        &dataGeneration::generate_uniform_array);
    break;
  }
  case dataGeneration::uniform_sorted: {
    input_distribution = "uniform sorted";
    auto [minval, maxval] = dataGeneration::get_default_limits<T>();
    extra_args = {std::make_pair("min_val", std::to_string(minval)),
                  std::make_pair("max_val", std::to_string(maxval))};
    std::function<void(T *, size_t, int, size_t)> f =
        [](T *A, size_t n, int seed, size_t discard) {
          return dataGeneration::generate_uniform_sorted_array(A, n, seed,
                                                               discard, false);
        };
    num_generator = f;
    break;
  }
  case dataGeneration::normal_sorted: {
    input_distribution = "normal sorted";
    extra_args = {std::make_pair("mean", "n/2"),
                  std::make_pair("std", "n*0.1")};
    std::function<void(T *, size_t, int, size_t)> f =
        [](T *A, size_t n, int seed, size_t discard) {
          return dataGeneration::generate_normal_sorted_array(A, n, seed,
                                                              discard, false);
        };
    num_generator = f;
    break;
  }
  case dataGeneration::uniform_reverse_sorted: {
    input_distribution = "uniform reverse sorted";
    auto [minval, maxval] = dataGeneration::get_default_limits<T>();
    extra_args = {std::make_pair("min_val", std::to_string(minval)),
                  std::make_pair("max_val", std::to_string(maxval))};
    std::function<void(T *, size_t, int, size_t)> f =
        [](T *A, size_t n, int seed, size_t discard) {
          return dataGeneration::generate_uniform_sorted_array(A, n, seed,
                                                               discard, true);
        };
    num_generator = f;
    break;
  }
  case dataGeneration::normal_reverse_sorted: {
    input_distribution = "normal reverse sorted";
    extra_args = {std::make_pair("mean", "n/2"),
                  std::make_pair("std", "n*0.1")};
    std::function<void(T *, size_t, int, size_t)> f =
        [](T *A, size_t n, int seed, size_t discard) {
          return dataGeneration::generate_normal_sorted_array(A, n, seed,
                                                              discard, true);
        };
    num_generator = f;
    break;
  }
  case dataGeneration::uniform_almost_sorted: {
    input_distribution = "uniform almost sorted";
    auto [minval, maxval] = dataGeneration::get_default_limits<T>();
    auto swap_proba = dataGeneration::DEFAULT_SWAP_PROBA;
    extra_args = {std::make_pair("min_val", std::to_string(minval)),
                  std::make_pair("max_val", std::to_string(maxval)),
                  std::make_pair("swap_proba", std::to_string(swap_proba))};
    std::function<void(T *, size_t, int, size_t)> f =
        [swap_proba](T *A, size_t n, int seed, size_t discard) {
          return dataGeneration::generate_uniform_almost_sorted_array(
              A, n, seed, discard, false, swap_proba);
        };
    num_generator = f;
    break;
  }
  case dataGeneration::normal_almost_sorted: {
    input_distribution = "normal almost sorted";
    auto swap_proba = dataGeneration::DEFAULT_SWAP_PROBA;
    extra_args = {std::make_pair("mean", "n/2"), std::make_pair("std", "n*0.1"),
                  std::make_pair("swap_proba", std::to_string(swap_proba))};
    std::function<void(T *, size_t, int, size_t)> f =
        [swap_proba](T *A, size_t n, int seed, size_t discard) {
          return dataGeneration::generate_normal_almost_sorted_array(
              A, n, seed, discard, false, swap_proba);
        };
    num_generator = f;
    break;
  }
  default: {
    std::cerr << dist << " Experiment not implemented!" << std::endl;
    exit(EXIT_FAILURE);
  }
  }

  InputDataInfo info;
  info.input_type = input_type;
  info.input_distribution = input_distribution;
  info.extra_args = extra_args;
  return std::make_pair(info, num_generator);
}

template <typename T>
Experiment<T> run_simple_experiment(int seed, DataDist dist,
                                    std::vector<int> sizes, int reps, size_t p,
                                    size_t k) {

  DEBUG(p);
  auto [info, num_generator] = get_default_distribution<T>(dist);
  Impls<T> impls = impls::sequential::get_all<T>(k);
  // Impls<T> pimpls = impls::parallel::get_all<T>(p, k);
  // impls.insert(impls.end(), pimpls.begin(), pimpls.end());
  // Impls<T> impls = impls::parallel::get_all<T>(p, k);
  Experiment<T> experiment =
      Experiment<T>(impls, info.input_type, info.input_distribution,
                    info.extra_args, seed, reps);
  run_implementations(experiment, num_generator, sizes);
  return experiment;
}

template <typename T>
Experiment<T> compare_implementations(Impls<T> &impls, int seed, DataDist dist,
                                      std::vector<int> sizes, int reps) {
  auto [info, num_generator] = get_default_distribution<T>(dist);
  Experiment<T> experiment =
      Experiment<T>(impls, info.input_type, info.input_distribution,
                    info.extra_args, seed, reps);
  run_implementations(experiment, num_generator, sizes);
  return experiment;
}

void wow(int reps) {
  using ttype = int;
  int seed = 0;
  int p = 32;
  size_t n = 16'000'000;
  int max_val = 1 << 19;
  ttype *std_sorted = new ttype[n];
  ttype *imp_sorted = new ttype[n];
  ttype *imp2_sorted = new ttype[n];
  int n_blocks;
  std::vector<size_t> steps;
  steps.push_back(8'000'000);
  steps.push_back(4'000'000);
  steps.push_back(2'000'000);
  steps.push_back(1'000'000);
  for (size_t step = 500'000; step > 32000; step -= 32000) {
    steps.push_back(step);
  }
  steps.push_back(32000);
  std::cout << "n, mwm, mas_ipx" << std::endl;
  for (auto step : steps) {
    for (int i = 0; i < reps; i++) {
      n_blocks = (n + step - 1) / step;
      std::cout << n_blocks << ",";
      dataGeneration::generate_uniform_array(std_sorted, n, seed, 0, 0,
                                             max_val);
      for (size_t l = 0; l < n; l += step) {
        size_t r = std::min(l + step, n);
        std::sort(std_sorted + l, std_sorted + r);
      }
      std::copy(std_sorted, std_sorted + n, imp_sorted);
      std::copy(std_sorted, std_sorted + n, imp2_sorted);
      timing2(bis::helpers::multiway_merge(imp_sorted, 0, n, step, p));
      timing2(bis::helpers::mergeAdjacentSequences(
          imp2_sorted, 0, n, step, bis::experimental::MergeStrat::inplace_mixed,
          p));
      // DEBUG(std::is_sorted(imp_sorted, imp_sorted + n));
      // DEBUG(std::is_sorted(imp2_sorted, imp2_sorted + n));
      std::cout << std::endl;
    }
  }
}

// TODO: Enable exps for other types
template Experiment<int> run_simple_experiment<int>(int seed, DataDist dist,
                                                    std::vector<int> sizes,
                                                    int reps, size_t p,
                                                    size_t k);

template Experiment<float> run_simple_experiment<float>(int seed, DataDist dist,
                                                        std::vector<int> sizes,
                                                        int reps, size_t p,
                                                        size_t k);

template Experiment<double>
run_simple_experiment<double>(int seed, DataDist dist, std::vector<int> sizes,
                              int reps, size_t p, size_t k);

template Experiment<size_t>
run_simple_experiment<size_t>(int seed, DataDist dist, std::vector<int> sizes,
                              int reps, size_t p, size_t k);

template Experiment<int> compare_implementations(Impls<int> &impls, int seed,
                                                 DataDist dist,
                                                 std::vector<int> sizes,
                                                 int reps);

template Experiment<size_t> compare_implementations(Impls<size_t> &impls,
                                                    int seed, DataDist dist,
                                                    std::vector<int> sizes,
                                                    int reps);

template Experiment<float> compare_implementations(Impls<float> &impls,
                                                   int seed, DataDist dist,
                                                   std::vector<int> sizes,
                                                   int reps);

template Experiment<double> compare_implementations(Impls<double> &impls,
                                                    int seed, DataDist dist,
                                                    std::vector<int> sizes,
                                                    int reps);
// template std::pair<InputDataInfo, GeneratorFunc<int>>
// get_default_distribution(DataDist dist);
} // namespace exps
