
#ifndef OMP_MERGESORT
#define OMP_MERGESORT
#include <cstddef>
namespace collected {

void mergesort_omp(int *S, size_t n, int p);

template <typename T> void mergeSortImproved(T *A, size_t n, int p);
template <typename T> void mergeSort(T *A, size_t n, int p);
} // namespace collected
#endif /* BASIC_DRF_H_ */
