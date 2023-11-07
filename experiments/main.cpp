#include "Block-InsertionSort/bis.hpp"
#include "Block-InsertionSort/experimental_pbis.hpp"
#include "algorithm"
#include "common/data_gen.hpp"
#include "common/experiments.hpp"
#include "common/impl_sorting_times.hpp"
#include "common/impls.hpp"
#include "common/parallel_impls.hpp"
#include "common/sequential_impls.hpp"
#include "common/utils.hpp"
#include <cstddef>
#include <iostream>
#include <string>
#include <utility>

struct ExpConfig {
  int reps = 10;
  int seed = 0;
  dataGeneration::DataDistribution dist = dataGeneration::uniform;
  int p = 32;
  int k = bis::DEFAULT_BLOCK_SIZE;
  int kp = 16;
  int kh = 8;
  size_t min_size = 1'000'000;
  size_t max_size = 5'000'000;
  size_t step = 1'000'000;
};

void running_exps(int k, int reps, int p);
void running_new_exps(int k, int reps, int p);

template <typename T>
void compare_hpbis_strats(ExpConfig &exp_config, bool generate_pdf,
                          std::string filename_prefix);

int main(int argc, char *argv[]) {
  int k;
  int num_threads;
  int reps;
  if (argc < 4) {
    std::cout << "Invalid number of arguments!" << std::endl;
    std::cout << "Usage ./prog <num_threads> <k> <reps>" << std::endl;
    exit(EXIT_FAILURE);
  }

  num_threads = atoi(argv[1]);
  if (num_threads < 1 || num_threads > 64) {
    num_threads = 16;
    std::cout << "Ignoring num_threads, using num_threads=" << num_threads
              << std::endl;
  }
  k = atoi(argv[2]);
  if (k < 4) {
    k = bis::DEFAULT_BLOCK_SIZE;
    std::cout << "Ignoring k, using k=" << k << std::endl;
  }

  reps = atoi(argv[3]);
  if (reps < 1 || reps > 100) {
    reps = 5;
    std::cout << "Ignoring reps, using num_reps=" << reps << std::endl;
  }

  // running_new_exps(k, reps, num_threads);
  running_exps(k, reps, num_threads);

  return 0;
}

std::vector<int> get_default_exp_sizes(size_t max_size, size_t step,
                                       size_t min_size) {
  std::vector<int> sizes;
  const std::vector<int> small_sizes = {1'000, 10'000, 100'000, 500'000};
  for (auto var : small_sizes) {
    if ((size_t)var > min_size) {
      sizes.push_back(var);
    }
  }
  for (size_t i = min_size; i <= max_size; i += step) {
    sizes.push_back(i);
  }
  return sizes;
}

template <typename T>
void compare_p_and_k(ExpConfig &exp_config, bool generate_pdf = false,
                     std::string filename_suffix = "") {
  std::string output_file = "final_hbis_p_and_k" + filename_suffix;
  impls::Impls<T> impls = {};
  std::vector<int> kvalues;
  std::vector<int> pvalues;
  kvalues.push_back(8);
  kvalues.push_back(16);
  // kvalues.push_back(32);

  // cores
  pvalues.push_back(2);
  pvalues.push_back(4);
  // pvalues.push_back(8);
  // pvalues.push_back(16);
  // pvalues.push_back(32);

  auto sizes = get_default_exp_sizes(exp_config.max_size, exp_config.step,
                                     exp_config.min_size);
  impls.push_back(impls::sequential::get_implementation<T>(impls::BiS));
  for (auto p : pvalues) {
    for (auto k : kvalues) {
      impls.push_back(
          impls::parallel::get_implementation<T>(impls::PBiS, p, k));
      impls.push_back(
          impls::parallel::get_implementation<T>(impls::hybrid_PBiS, p, k));
    }
  }
  Experiment<T> result = exps::compare_implementations(
      impls, exp_config.seed, exp_config.dist, sizes, exp_config.reps);
  std::string filename =
      output_file + "_" + std::to_string(std::time(nullptr)) + ".json";
  result.dump_to_file(filename);
  if (generate_pdf) {
    result.generate_pdf_plot(output_file, true);
  }
}

template <typename T>
void compare_parallel_impls(ExpConfig &exp_config, bool generate_pdf,
                            std::string filename_prefix) {
  std::string pstring = std::to_string(exp_config.p);
  std::string distring =
      dataGeneration::getDataDistributionName(exp_config.dist);
  std::string output_file = filename_prefix + distring + "_" + pstring + "p";
  auto sizes = get_default_exp_sizes(exp_config.max_size, exp_config.step,
                                     exp_config.min_size);

  impls::Impls<T> impls =
      impls::parallel::get_all<T>(exp_config.p, exp_config.k);

  auto result = exps::compare_implementations(
      impls, exp_config.seed, exp_config.dist, sizes, exp_config.reps);
  std::string filename =
      output_file + "_" + std::to_string(std::time(nullptr)) + ".json";
  result.dump_to_file(filename);
  if (generate_pdf) {
    result.generate_pdf_plot(output_file, true);
  }
  return;
}

template <typename T>
void compare_hpbis_strats(ExpConfig &exp_config, bool generate_pdf,
                          std::string filename_prefix) {
  std::string pstring = std::to_string(exp_config.p);
  std::string distring =
      dataGeneration::getDataDistributionName(exp_config.dist);
  std::string output_file = filename_prefix + distring + "_" + pstring;
  auto sizes = get_default_exp_sizes(exp_config.max_size, exp_config.step,
                                     exp_config.min_size);

  impls::Impls<T> impls;

  impls.push_back(impls::parallel::get_implementation<T>(
      impls::AvailableImplementation::hybrid_PBiS, exp_config.p, 8, 1));
  impls.push_back(impls::parallel::get_implementation<T>(
      impls::AvailableImplementation::hybrid_PBiS, exp_config.p, 8, 2));
  impls.push_back(impls::parallel::get_implementation<T>(
      impls::AvailableImplementation::hybrid_PBiS, exp_config.p, 16, 1));
  impls.push_back(impls::parallel::get_implementation<T>(
      impls::AvailableImplementation::hybrid_PBiS, exp_config.p, 16, 2));
  impls.push_back(impls::parallel::get_implementation<T>(
      impls::AvailableImplementation::hybrid_PBiS_split, exp_config.p, 8, 1));
  impls.push_back(impls::parallel::get_implementation<T>(
      impls::AvailableImplementation::hybrid_PBiS_split, exp_config.p, 8, 2));

  auto result = exps::compare_implementations(
      impls, exp_config.seed, exp_config.dist, sizes, exp_config.reps);
  std::string filename =
      output_file + "_" + std::to_string(std::time(nullptr)) + ".json";
  result.dump_to_file(filename);
  if (generate_pdf) {
    result.generate_pdf_plot(output_file, true);
  }
}

template <typename T>
void find_sequential_threshold(ExpConfig &exp_config, bool generate_pdf,
                               std::string filename_prefix) {
  std::string pstring = std::to_string(exp_config.p);
  std::string distring =
      dataGeneration::getDataDistributionName(exp_config.dist);
  std::string output_file = filename_prefix + distring + "_" + pstring;
  std::vector<int> sizes;

  for (int i = 8; i < 18; i++) {
    sizes.push_back((2 << i));
  }

  impls::Impls<T> impls;

  impls.push_back(impls::parallel::get_implementation<T>(
      impls::AvailableImplementation::hybrid_PBiS, 2, exp_config.kh));

  impls.push_back(impls::parallel::get_implementation<T>(
      impls::AvailableImplementation::hybrid_PBiS, 4, exp_config.kh));

  impls.push_back(impls::parallel::get_implementation<T>(
      impls::AvailableImplementation::hybrid_PBiS, 8, exp_config.kh));

  impls.push_back(impls::sequential::get_implementation<T>(
      impls::AvailableImplementation::BiS, exp_config.k));

  auto result = exps::compare_implementations(
      impls, exp_config.seed, exp_config.dist, sizes, exp_config.reps);
  std::string filename =
      output_file + "_" + std::to_string(std::time(nullptr)) + ".json";
  result.dump_to_file(filename);
  if (generate_pdf) {
    result.generate_pdf_plot(output_file, true);
  }
}

void running_new_exps(int k, int reps, int p) {
  std::cout << "Running experiments" << std::endl;
  bool generate_pdf = p < 8; // false
  ExpConfig exp_config;
  exp_config.reps = reps;
  exp_config.k = k;
  exp_config.p = p;
  exp_config.max_size = 40'000'000;
  std::string filename_prefix = "final_strats_";
  compare_hpbis_strats<int>(exp_config, generate_pdf, filename_prefix);
  compare_p_and_k<int>(exp_config);
  // filename_prefix = "sambombaso_";
  // find_sequential_threshold<int>(exp_config, generate_pdf, filename_prefix);
}

void running_exps(int k, int reps, int p) {
  std::cout << "Running experiments" << std::endl;
  bool generate_pdf = false; // false
  ExpConfig exp_config;
  exp_config.reps = reps;
  exp_config.k = k;
  exp_config.p = p;

  std::string date_prefix = "0611";
  // compare_p_and_k<int>(exp_config, generate_pdf);

  std::string prefix;

  dataGeneration::DataDistribution main_dists[] = {
      dataGeneration::DataDistribution::normal,
      dataGeneration::DataDistribution::uniform,
  };
  dataGeneration::DataDistribution int_dists[] = {
      dataGeneration::DataDistribution::normal,
      dataGeneration::DataDistribution::normal_almost_sorted,
      dataGeneration::DataDistribution::normal_sorted,
      dataGeneration::DataDistribution::uniform,
      dataGeneration::DataDistribution::uniform_reverse_sorted,
  };

  prefix = date_prefix + "_float_";
  for (auto dist : main_dists) {
    exp_config.dist = dist;
    compare_parallel_impls<float>(exp_config, generate_pdf, prefix);
  }

  prefix = date_prefix + "_int_";
  for (auto dist : int_dists) {
    exp_config.dist = dist;
    compare_parallel_impls<int>(exp_config, generate_pdf, prefix);
  }

  prefix = date_prefix + "_ulong_";
  for (auto dist : int_dists) {
    exp_config.dist = dist;
    compare_parallel_impls<size_t>(exp_config, generate_pdf, prefix);
  }

  prefix = date_prefix + "_double_";
  for (auto dist : main_dists) {
    exp_config.dist = dist;
    compare_parallel_impls<double>(exp_config, generate_pdf, prefix);
  }
}
