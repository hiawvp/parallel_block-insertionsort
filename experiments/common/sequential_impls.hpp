#ifndef SIMPLS
#define SIMPLS

#include "impl_sorting_times.hpp"
#include "impls.hpp"
#include "utils.hpp"
#include <cstddef>
#include <string>
#include <vector>

namespace impls {
namespace sequential {
static const AvailableImplementation All[] = {
    BiS, std_qsort, ips4o_sort, ips4o_sort, std_sort,
};

template <typename T>
Impl<T> get_implementation(AvailableImplementation val, int k = kdepana);

template <typename T> Impls<T> get_all(int k);

} // namespace sequential
} // namespace impls
#endif /* BASIC_DRF_H_ */
