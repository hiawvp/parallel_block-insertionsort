#include "parallel_impls.hpp"
#include "Block-InsertionSort/experimental_pbis.hpp"
#include "Block-InsertionSort/pbis.hpp"
// #include "algorithm"
#include "common/impls.hpp"
#include "common/utils.hpp"
#include "external/mergesort/mergesort.hpp"
#include "ips2ra.hpp"
#include "ips4o.hpp"
#include <algorithm>
#include <cstddef>
#include <execution>
#include <omp.h>
#include <oneapi/tbb/parallel_sort.h>
#include <parallel/algorithm>
#include <parallel/settings.h>
#include <parallel/tags.h>
#include <tbb/global_control.h>

namespace impls {

namespace parallel {

Impl<size_t> get_ips2ra(int p) {
  std::string impl_name = "ips2ra::parallel::sort";
  ImplExtraArgs extra_args;
  DEBUG(impl_name);
  std::function<void(size_t *, size_t)> f = [p](size_t *A, size_t n) {
    omp_set_num_threads(p);
    tbb::global_control c(tbb::global_control::max_allowed_parallelism, p);
    ips2ra::parallel::sort(A, A + n);
    return;
  };
  extra_args = {std::make_pair("p", p)};
  return SortingImplementationDetails(impl_name, extra_args, f);
}

template <typename T> Impls<T> __get_all(__PImplPool pool, int p, int k) {
  Impls<T> impls = {};
  switch (pool) {
  case NonRadixImpls: {
    for (const auto e : NonRadix) {
      impls.push_back(get_implementation<T>(e, p, k));
    }
    break;
  };
  case NonBiSNonRadixImpls: {
    for (const auto e : NonBiSNonRadix) {
      impls.push_back(get_implementation<T>(e, p, k));
    }
    break;
  };
  }
  return impls;
}

template <> Impls<size_t> __get_all(__PImplPool pool, int p, int k) {
  Impls<size_t> impls = {};
  switch (pool) {
  case NonRadixImpls: {
    for (const auto e : NonRadix) {
      impls.push_back(get_implementation<size_t>(e, p, k));
    }
    break;
  };
  case NonBiSNonRadixImpls: {
    for (const auto e : NonBiSNonRadix) {
      impls.push_back(get_implementation<size_t>(e, p, k));
    }
    break;
  };
  }
  impls.push_back(get_ips2ra(p));
  return impls;
}

template <typename T> Impls<T> get_all(int p, int k) {
  return __get_all<T>(NonRadixImpls, p, k);
}

template <typename T> Impls<T> get_non_bis_impls(int p) {
  return __get_all<T>(NonBiSNonRadixImpls, p, 1);
}

template <typename T>
Impl<T> get_implementation(AvailableImplementation val, int p, int k,
                           int method) {
  std::string impl_name;
  ImplExtraArgs extra_args;
  switch (val) {
  case PBiS: {
    impl_name = "Parallel BiS";
    DEBUG(impl_name);
    std::function<void(T *, size_t)> f = [k, p](T *A, size_t n) {
      // return pbis::pblock_insertion_sort(A, n, p, k);
      return bis::parallel::block_insertion_sort(A, n, p, k);
    };
    extra_args = {std::make_pair("k", k), std::make_pair("p", p)};
    return SortingImplementationDetails(impl_name, extra_args, f);
  }
  case hybrid_PBiS: {
    impl_name = "Hybrid PBiS";
    DEBUG(impl_name);
    std::function<void(T *, size_t)> f = [k, p, method](T *A, size_t n) {
      return bis::experimental::bis_iterative(A, 0, n, p, k, method);
    };
    extra_args = {std::make_pair("k", k), std::make_pair("p", p)};
    if (method != 0) {
      extra_args.push_back(std::make_pair("st", method));
    }
    return SortingImplementationDetails(impl_name, extra_args, f);
  }
  case hybrid_PBiS_split: {
    impl_name = "Hybrid PBiS PSplit";
    DEBUG(impl_name);
    std::function<void(T *, size_t)> f = [k, p, method](T *A, size_t n) {
      return bis::experimental::pbis_task_split(A, n, p, k, method);
    };
    extra_args = {std::make_pair("k", k), std::make_pair("p", p)};
    if (method != 0) {
      extra_args.push_back(std::make_pair("st", method));
    }
    return SortingImplementationDetails(impl_name, extra_args, f);
  }
  case omp_task_PBiS: {
    impl_name = "OMP Task PBiS";
    DEBUG(impl_name);
    std::function<void(T *, size_t)> f = [k, p](T *A, size_t n) {
      return bis::experimental::rec_bis_general2(A, 0, n, p, k);
      // return bis::parallel::rec_block_insertion_sort(A, 0, n, p, k);
    };
    extra_args = {std::make_pair("k", k), std::make_pair("p", p)};
    return SortingImplementationDetails(impl_name, extra_args, f);
  }
  case std_sort_par: {
    impl_name = "std::sort execution::par";
    DEBUG(impl_name);
    std::function<void(T *, size_t)> f = [p](T *A, size_t n) {
      tbb::global_control c(tbb::global_control::max_allowed_parallelism, p);
      std::sort(std::execution::par, A, A + n);
      return;
    };
    extra_args = {std::make_pair("p", p)};
    return SortingImplementationDetails(impl_name, extra_args, f);
  }

  case gnu_par_QS: {
    impl_name = "gnu_parallel::QS";
    DEBUG(impl_name);
    std::function<void(T *, size_t)> f = [p](T *A, size_t n) {
      // tbb::global_control c(tbb::global_control::max_allowed_parallelism, p);
      std::__parallel::sort(A, A + n, __gnu_parallel::quicksort_tag(p));
      return;
    };
    extra_args = {std::make_pair("p", p)};
    return SortingImplementationDetails(impl_name, extra_args, f);
  }

  case gnu_par_QSB: {
    impl_name = "gnu_parallel::QSB";
    DEBUG(impl_name);
    std::function<void(T *, size_t)> f = [p](T *A, size_t n) {
      // tbb::global_control c(tbb::global_control::max_allowed_parallelism, p);
      std::__parallel::sort(A, A + n,
                            __gnu_parallel::balanced_quicksort_tag(p));
      return;
    };
    extra_args = {std::make_pair("p", p)};
    return SortingImplementationDetails(impl_name, extra_args, f);
  }

  case gnu_par_MWMS: {
    impl_name = "gnu_parallel::MWMS";
    DEBUG(impl_name);
    std::function<void(T *, size_t)> f = [p](T *A, size_t n) {
      // tbb::global_control c(tbb::global_control::max_allowed_parallelism, p);
      std::__parallel::sort(A, A + n,
                            __gnu_parallel::multiway_mergesort_tag(p));
      return;
    };
    extra_args = {std::make_pair("p", p)};
    return SortingImplementationDetails(impl_name, extra_args, f);
  }
  case tbb_parallel_sort: {
    impl_name = "tbb::parallel_sort";
    DEBUG(impl_name);
    std::function<void(T *, size_t)> f = [p](T *A, size_t n) {
      tbb::global_control c(tbb::global_control::max_allowed_parallelism, p);
      tbb::parallel_sort(A, A + n);
      return;
    };
    extra_args = {std::make_pair("p", p)};
    return SortingImplementationDetails(impl_name, extra_args, f);
  }

  case ips4o_parallel_sort: {
    impl_name = "ips4o::parallel::sort";
    DEBUG(impl_name);
    std::function<void(T *, size_t)> f = [p](T *A, size_t n) {
      omp_set_num_threads(p);
      tbb::global_control c(tbb::global_control::max_allowed_parallelism, p);
      ips4o::parallel::sort(A, A + n);
      return;
    };
    extra_args = {std::make_pair("p", p)};
    return SortingImplementationDetails(impl_name, extra_args, f);
  }

  case omp_task_mergesort: {
    impl_name = "OMP Task MergeSort";
    DEBUG(impl_name);
    std::function<void(T *, size_t)> f = [p](T *A, size_t n) {
      return collected::mergeSortImproved(A, n, p);
    };
    extra_args = {std::make_pair("p", p)};
    return SortingImplementationDetails(impl_name, extra_args, f);
  }

  default: {
    std::cerr << "ERROR! Implementacion no disponible!" << std::endl;
    exit(EXIT_FAILURE);
  }
  }
}

template Impl<int> get_implementation<int>(AvailableImplementation val, int p,
                                           int k, int method);

template Impl<float> get_implementation<float>(AvailableImplementation val,
                                               int p, int k, int method);

template Impl<long> get_implementation<long>(AvailableImplementation val, int p,
                                             int k, int method);
template Impl<double> get_implementation<double>(AvailableImplementation val,
                                                 int p, int k, int method);
template Impl<size_t> get_implementation<size_t>(AvailableImplementation val,
                                                 int p, int k, int method);

template Impls<int> get_all<int>(int p, int k);
template Impls<float> get_all<float>(int p, int k);
template Impls<double> get_all<double>(int p, int k);
template Impls<long> get_all<long>(int p, int k);
template Impls<size_t> get_all<size_t>(int p, int k);

template Impls<int> get_non_bis_impls<int>(int p);
template Impls<float> get_non_bis_impls<float>(int p);
template Impls<double> get_non_bis_impls<double>(int p);
template Impls<long> get_non_bis_impls<long>(int p);
template Impls<size_t> get_non_bis_impls<size_t>(int p);

} // namespace parallel
} // namespace impls
