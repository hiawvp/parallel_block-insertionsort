#ifndef PARALLEL_BIS
#define PARALLEL_BIS
#include "bis.hpp"
#include <iostream>

// TODO: templatizar bis
// implementaciones personales de bis
namespace bis {
namespace parallel {

template <typename T>
void block_insertion_sort(T *S, size_t n, int p = 1,
                          int k = DEFAULT_BLOCK_SIZE);
template <typename T>
void rec_block_insertion_sort(T *S, size_t l, size_t r, int p = 1,
                              int k = DEFAULT_BLOCK_SIZE);
} // namespace parallel
} // namespace bis
#endif /* BASIC_DRF_H_ */
