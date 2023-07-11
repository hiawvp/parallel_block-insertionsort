#include "pbis.hpp"
#include "pbis_example_utils.hpp"
#include <algorithm>
#include <cstddef>
#include <cstdlib>
#include <execution>
#include <functional>
#include <iostream>
#include <parallel/multiway_merge.h>
#include <parallel/settings.h>
#include <parallel/types.h>

#define DEFAULT_SIZE 7'777'777
#define DEFAULT_SEED 0

int main(int argc, char *argv[]) {
  srand(DEFAULT_SEED);
  size_t n = DEFAULT_SIZE;
  using datatype = float;
  datatype *A = new datatype[n];
  datatype *sorted_copy = new datatype[n];
  datatype exit_code = 0;
  auto less_comp = std::less<datatype>();
  auto greater_comp = std::greater<datatype>();

  gen_array(A, n);
  std::copy_n(A, n, sorted_copy);
  timeit(std::sort(sorted_copy, sorted_copy + n), "std::sort");
  timeit(bis::sort(A, A + n), "bis::sort");
  exit_code = exit_code || check_results(A, sorted_copy, n);

  gen_array(A, n);
  std::copy_n(A, n, sorted_copy);
  timeit(std::sort(sorted_copy, sorted_copy + n), "std::sort");
  timeit(bis::parallel::sort(A, A + n), "bis::parallel::sort");
  exit_code = exit_code || check_results(A, sorted_copy, n);

  gen_array(A, n);
  std::copy_n(A, n, sorted_copy);
  timeit(std::sort(sorted_copy, sorted_copy + n, less_comp), "std::sort");
  timeit(bis::parallel::sort(A, A + n, less_comp), "bis::parallel::sort");
  exit_code = exit_code || check_results(A, sorted_copy, n);

  std::vector<float> V(n);
  std::vector<float> SC(n);
  gen_array(V.data(), n);
  std::copy_n(V.begin(), n, SC.begin());
  timeit(std::sort(SC.begin(), SC.end()), "std::sort");
  timeit(bis::parallel::sort(V.begin(), V.end()), "bis::sort");
  exit_code = exit_code || check_results(V.data(), SC.data(), n);

#if defined(_USETBB)

  std::cout << "tbb enabled" << std::endl;
  __gnu_parallel::_Settings s;
  s.sort_algorithm = __gnu_parallel::_SortAlgorithm::QS;
  __gnu_parallel::_Settings::set(s);

  gen_array(A, n);
  std::copy_n(A, n, sorted_copy);
  timeit(std::sort(std::execution::par, sorted_copy, sorted_copy + n,
                   greater_comp),
         "std::sort(gnu_par::QS)");
  timeit(bis::parallel::sort(A, A + n, greater_comp), "bis::parallel::sort");
  exit_code = exit_code || check_results(A, sorted_copy, n);

  s.sort_algorithm = __gnu_parallel::_SortAlgorithm::QS_BALANCED;
  __gnu_parallel::_Settings::set(s);

  gen_array(A, n);
  std::copy_n(A, n, sorted_copy);
  timeit(std::sort(std::execution::par, sorted_copy, sorted_copy + n,
                   greater_comp),
         "std::sort(gnu_par::QSB)");
  timeit(bis::parallel::sort(A, A + n, greater_comp), "bis::parallel::sort");
  exit_code = exit_code || check_results(A, sorted_copy, n);

  s.sort_algorithm = __gnu_parallel::_SortAlgorithm::MWMS;
  __gnu_parallel::_Settings::set(s);

  gen_array(A, n);
  std::copy_n(A, n, sorted_copy);
  timeit(std::sort(std::execution::par, sorted_copy, sorted_copy + n,
                   greater_comp),
         "std::sort(gnu_par::MWMS)");
  timeit(bis::parallel::sort(A, A + n, greater_comp), "bis::parallel::sort");
  exit_code = exit_code || check_results(A, sorted_copy, n);

#endif // _USETBB

  // print_array(A, n);
  // print_array(sorted_copy, n);

  // std::cout << "Press Enter to close" << std::endl;
  // std::cin.ignore();

  return exit_code;
}
