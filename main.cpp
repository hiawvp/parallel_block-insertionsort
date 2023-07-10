#include "pbis.hpp"
#include "pbis_example_utils.hpp"
#include <algorithm>
#include <cstddef>
#include <cstdlib>
#include <execution>
#include <functional>
#include <iostream>

#define DEFAULT_SIZE 7'777'777
#define DEFAULT_SEED 0

int main(int argc, char *argv[]) {
  srand(DEFAULT_SEED);
  size_t n = DEFAULT_SIZE;
  int *A = new int[n];
  int *sorted_copy = new int[n];
  int exit_code = 0;
  auto less_comp = std::less<int>();
  auto greater_comp = std::greater<int>();

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

#if defined(_USETBB)
  gen_array(A, n);
  std::copy_n(A, n, sorted_copy);
  timeit(std::sort(std::execution::par, sorted_copy, sorted_copy + n,
                   greater_comp),
         "std::sort(exec::par)");
  timeit(bis::parallel::sort(A, A + n, greater_comp), "bis::parallel::sort");
  exit_code = exit_code || check_results(A, sorted_copy, n);

#endif // _USETBB
  gen_array(A, n);
  std::copy_n(A, n, sorted_copy);
  timeit(std::sort(sorted_copy, sorted_copy + n, less_comp), "std::sort");
  timeit(bis::parallel::sort(A, A + n, less_comp), "bis::parallel::sort");
  exit_code = exit_code || check_results(A, sorted_copy, n);

  // forced bad result:
  gen_array(A, n);
  std::copy_n(A, n, sorted_copy);
  timeit(std::sort(sorted_copy, sorted_copy + n, greater_comp), "std::sort");
  timeit(bis::sort(A, A + n, less_comp), "bis::sort");
  exit_code = exit_code || check_results(A, sorted_copy, n);

  // print_array(A, n);
  // print_array(sorted_copy, n);

  std::cout << "Press Enter to close" << std::endl;
  std::cin.ignore();

  return exit_code;
}
