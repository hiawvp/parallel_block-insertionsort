#ifndef IMPLS
#define IMPLS

#include "impl_sorting_times.hpp"
#include <cstddef>
#include <string>
#include <vector>

namespace impls {
template <typename T>
using Impls = std::vector<SortingImplementationDetails<T>>;
using ImplExtraArgs = std::vector<std::pair<std::string, int>>;

template <typename T> using Impl = SortingImplementationDetails<T>;
const int kdepana = 16;
const int peppy = 16;

enum AvailableImplementation {
  BiS,
  PBiS,
  hybrid_PBiS,
  hybrid_PBiS_split,
  omp_task_PBiS,
  std_sort_par,
  gnu_par_QS,
  gnu_par_QSB,
  gnu_par_MWMS,
  tbb_parallel_sort,
  ips4o_parallel_sort,
  ips2ra_parallel_sort,
  omp_task_mergesort,
  // recursive_BiS,
  // recursive_PBiS,
  // recursive_PBiS_test,
  MyBiS,
  std_sort,
  std_qsort,
  ips4o_sort,
  // parallel_mergesort,
};

} // namespace impls
#endif /* BASIC_DRF_H_ */
