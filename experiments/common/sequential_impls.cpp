#include "sequential_impls.hpp"
#include "Block-InsertionSort/bis.hpp"
#include "Block-InsertionSort/original_bis/Sorting.hpp"
#include "algorithm"
#include "common/impls.hpp"
#include "cstdlib"
#include "ips4o.hpp"
#include <cstddef>
#include <cstdlib>

namespace impls {
namespace {
int compare(const void *a, const void *b) { return (*(int *)a - *(int *)b); }
} // namespace

namespace sequential {

template <typename T> Impls<T> get_all(int k) {
  DEBUG(k);
  Impls<T> impls = {};
  for (const auto e : All) {
    impls.push_back(get_implementation<T>(e, k));
  }
  return impls;
}

template <typename T>
Impl<T> get_implementation(AvailableImplementation val, int k) {
  ImplExtraArgs extra_args;
  std::string impl_name;
  switch (val) {
  case BiS: {
    impl_name = "Original BiS";
    DEBUG2(impl_name);
    std::function<void(T *, size_t)> f = [](T *A, size_t n) {
      return srt::blockInsertionSort<T>(A, 0, n - 1);
    };
    extra_args = {std::make_pair("k", srt::blockIS)};
    return SortingImplementationDetails(impl_name, extra_args, f);
  }

  case std_qsort: {
    impl_name = "std::qsort";
    DEBUG2(impl_name);
    std::function<void(T *, size_t)> f = [](T *A, size_t n) {
      return std::qsort(A, n, sizeof(T), compare);
    };
    extra_args = {};
    return SortingImplementationDetails(impl_name, extra_args, f);
  }

  case std_sort: {
    impl_name = "std::sort";
    DEBUG2(impl_name);
    std::function<void(T *, size_t)> f = [](T *A, size_t n) {
      return std::sort(A, A + n);
    };
    extra_args = {};
    return SortingImplementationDetails(impl_name, extra_args, f);
  }

  case ips4o_sort: {
    impl_name = "ips4o::sort";
    DEBUG2(impl_name);
    std::function<void(T *, size_t)> f = [](T *A, size_t n) {
      return ips4o::sort(A, A + n);
    };
    extra_args = {};
    return SortingImplementationDetails(impl_name, extra_args, f);
  }
  default: {
    std::cerr << "ERROR! Implementacion no disponible!" << std::endl;
    exit(EXIT_FAILURE);
  }
  }
}
template Impl<int> get_implementation<int>(AvailableImplementation val, int k);

template Impl<float> get_implementation<float>(AvailableImplementation val,
                                               int k);

template Impl<long> get_implementation<long>(AvailableImplementation val,
                                             int k);
template Impl<double> get_implementation<double>(AvailableImplementation val,
                                                 int k);
template Impl<size_t> get_implementation<size_t>(AvailableImplementation val,
                                                 int k);
template Impls<int> get_all(int k);
template Impls<long> get_all(int k);
template Impls<double> get_all(int k);
template Impls<float> get_all(int k);
template Impls<size_t> get_all(int k);
} // namespace sequential
} // namespace impls
