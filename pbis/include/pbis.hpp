#ifndef EXPERIMENTAL_PBIS
#define EXPERIMENTAL_PBIS
#include "pbis/bis_par.hpp"
#include "pbis/bis_seq.hpp"
#include <cstddef>
#include <functional>
#include <iostream>
#include <iterator>

namespace bis {

namespace {

template <typename RandomIT, typename Compare>
void invoke_parallel_sort(RandomIT first, RandomIT last, Compare comp) {
  size_t n = std::distance(first, last);
  if (n <= bis::SEQUENTIAL_THRESHOLD) {
    // return bis::sequential_sorter::block_insertion_sort(first, n, comp);
    return bis::sequential_sorter::block_insertion_sort(&(*first), n, comp);
  } else {
    // return bis::parallel_sorter::block_insertion_sort(&(*first), n, comp);
    return bis::parallel_sorter::hybrid_pbis_psplip(&(*first), n, comp);
  }
}

template <typename RandomIT, typename Compare>
void invoke_sequential_sort(RandomIT first, RandomIT last, Compare comp) {
  size_t n = std::distance(first, last);
  // return bis::sequential_sorter::block_insertion_sort(first, n, comp);
  return bis::sequential_sorter::block_insertion_sort(&(*first), n, comp);
}

} // namespace

namespace parallel {
template <typename RandomIT> void sort(RandomIT first, RandomIT last) {
  using ValueType = typename std::iterator_traits<RandomIT>::value_type;
  auto comp = std::less<ValueType>();
  bis::invoke_parallel_sort(first, last, comp);
}

template <typename RandomIT, typename Compare>
void sort(RandomIT first, RandomIT last, Compare comp) {
  bis::invoke_parallel_sort(first, last, comp);
}
} // namespace parallel

template <typename RandomIT> void sort(RandomIT first, RandomIT last) {
  using ValueType = typename std::iterator_traits<RandomIT>::value_type;
  auto comp = std::less<ValueType>();
  bis::invoke_sequential_sort(first, last, comp);
}
template <typename RandomIT, typename Compare>
void sort(RandomIT first, RandomIT last, Compare comp) {
  bis::invoke_sequential_sort(first, last, comp);
}
} // namespace bis

#endif /* BASIC_DRF_H_ */
