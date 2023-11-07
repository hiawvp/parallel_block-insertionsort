#ifndef EXPERIMENTAL_PBIS
#define EXPERIMENTAL_PBIS
#include "bis.hpp"
#include <cstddef>
#include <iostream>

// TODO: templatizar bis
// implementaciones personales de bis
namespace bis {
namespace experimental {

enum MergeStrat {
  // classic,
  bis,
  inplace,
  inplace_pfor,
  inplace_exec_policy,
  inplace_mixed,
  inplace_mixed_adaptive,
  out_of_place,
  gnu_multiway,
  // ptask,
};

static const MergeStrat All[] = {
    bis,           inplace,
    inplace_pfor,  inplace_exec_policy,
    inplace_mixed, inplace_mixed_adaptive,
    out_of_place,  gnu_multiway,
    // ptask
};

static const std::string AllNames[] = {
    "bis",           "inplace",
    "inplace_pfor",  "inplace_exec_policy",
    "inplace_mixed", "inplace_mixed_adaptive",
    "out_of_place",  "gnu_multiway",
    // ptask
};
template <typename T>
void rec_block_insertion_sort(
    T *S, size_t l, size_t r, int p = 1, int k = DEFAULT_BLOCK_SIZE,
    MergeStrat merge_strat = MergeStrat::inplace_mixed);

template <typename T>
void rec_bis_general(T *S, size_t l, size_t r, int p, int k,
                     MergeStrat merge_strat = MergeStrat::inplace_mixed);

template <typename T>
void rec_bis_general2(T *S, size_t l, size_t r, int p, int k,
                      MergeStrat merge_strat = MergeStrat::inplace_mixed);

template <typename T>
void rec_bis_general3(T *S, size_t l, size_t r, int p, int k);

template <typename T>
void bis_iterative(T *S, size_t l, size_t r, int p, int k, int method = 0);

template <typename T>
void pbis_task_split(T *S, size_t n, int p, int k, int method = 0);

template <typename T>
void rec_bis_strats(T *S, size_t l, size_t r, int p, int k);

} // namespace experimental
} // namespace bis
#endif /* BASIC_DRF_H_ */
