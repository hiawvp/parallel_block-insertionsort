
#ifndef EXPERIMENTS
#define EXPERIMENTS

#include "data_gen.hpp"
#include "impl_sorting_times.hpp"
#include "vector"
#include <cstddef>

namespace exps {

template <typename T>
using GeneratorFunc = std::function<void(T *, size_t, int, size_t)>;

template <typename T>
using Impls = std::vector<SortingImplementationDetails<T>>;

using DataDist = dataGeneration::DataDistribution;

struct InputDataInfo {
  std::string input_type;
  std::string input_distribution;
  std::vector<std::pair<std::string, std::string>> extra_args;
};

// template <typename T>
// void run_implementations(int seed, size_t reps, Experiment<T> experiment,
//                          Impls<T> &impls, GeneratorFunc<T> num_generator);
//
template <typename T>
Experiment<T> run_simple_experiment(int seed, DataDist dist,
                                    std::vector<int> sizes, int reps, size_t p,
                                    size_t k);

template <typename T>
Experiment<T> compare_implementations(Impls<T> &impls, int seed, DataDist dist,
                                      std::vector<int> sizes, int reps);

void wow(int reps = 10);
} // namespace exps

#endif /* BASIC_DRF_H_ */
