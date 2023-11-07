#ifndef PIMPLS
#define PIMPLS

#include "common/impl_sorting_times.hpp"
#include "common/impls.hpp"

namespace impls {
namespace parallel {
static const AvailableImplementation All[] = {
    PBiS,
    hybrid_PBiS,
    hybrid_PBiS_split,
    gnu_par_QS,
    gnu_par_MWMS,
    // omp_task_PBiS,
    // std_sort_par,
    // gnu_par_QSB,
    tbb_parallel_sort,
    ips4o_parallel_sort,
    omp_task_mergesort,
    ips2ra_parallel_sort,
};

static const AvailableImplementation NonRadix[] = {
    PBiS,
    hybrid_PBiS,
    hybrid_PBiS_split,
    gnu_par_QS,
    gnu_par_MWMS,
    // gnu_par_QSB,
    // omp_task_PBiS,
    // std_sort_par,
    tbb_parallel_sort,
    ips4o_parallel_sort,
    omp_task_mergesort,
};

static const AvailableImplementation NonBiSNonRadix[] = {
    // gnu_par_QS, gnu_par_MWMS,
    // ips4o_parallel_sort,
    // tbb_parallel_sort,
    // ips4o_parallel_sort,
    // omp_task_mergesort,
};

enum __PImplPool {
  NonRadixImpls,
  NonBiSNonRadixImpls,
};

template <typename T>
Impl<T> get_implementation(AvailableImplementation val, int p = peppy,
                           int k = kdepana, int method = 0);
template <typename T> Impls<T> get_all(int p = peppy, int k = kdepana);
template <typename T> Impls<T> get_non_bis_impls(int p = peppy);

Impl<unsigned long> get_ips2ra(int p);

} // namespace parallel
} // namespace impls
#endif /* BASIC_DRF_H_ */
