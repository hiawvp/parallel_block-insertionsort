#include "utils.hpp"
#include <iostream>
#include <random>
#include <vector>

namespace utils {

void print_sep() {
  if (!DEBUG_FLAG) {
    return;
  }
  std::cout << '\n';
  for (size_t i = 0; i < 55; i++) {
    std::cout << '-';
  }
  std::cout << "\n\n";
}

std::string get_typename(DataTypes val) {
  switch (val) {
  case int_data:
    return "int";
  case float_data:
    return "float";
  case long_data:
    return "long";
  case double_data:
    return "double";
  default:
    return "unknown";
  }
}

template <typename T> void print_array(T *A, size_t n, std::string hint) {
  if (!DEBUG_FLAG) {
    return;
  }
  if (n > MAX_PRINT || n < 1) {
    // std::cout << "Array too big qq" << endl;
    return;
  }
  std::cout << hint << std::endl;
  if (n == 1) {
    std::cout << "[" << A[0] << "]" << std::endl;
    return;
  }
  std::cout << "[ ";
  for (size_t i = 0; i < n - 1; i++) {
    std::cout << A[i] << ", ";
  }
  std::cout << A[n - 1] << " ]\n";
}

template <typename T>
bool check_sort_result(T *std_sort, T *my_sort, size_t n) {
  for (size_t i = 0; i < n; i++) {
    if (std_sort[i] != my_sort[i]) {
      std::cout << "ERROR EN EL INDICE " << i << std::endl;
      return false;
    }
  }
  return true;
}

template void print_array<int>(int *A, size_t n, std::string hint);
template void print_array<float>(float *A, size_t n, std::string hint);
template void print_array<long>(long *A, size_t n, std::string hint);
template void print_array<double>(double *A, size_t n, std::string hint);
template void print_array<size_t>(size_t *A, size_t n, std::string hint);

template bool check_sort_result<int>(int *std_sort, int *my_sort, size_t n);
template bool check_sort_result<float>(float *std_sort, float *my_sort,
                                       size_t n);
template bool check_sort_result<long>(long *std_sort, long *my_sort, size_t n);
template bool check_sort_result<double>(double *std_sort, double *my_sort,
                                        size_t n);
template bool check_sort_result<size_t>(size_t *std_sort, size_t *my_sort,
                                        size_t n);

} // namespace utils
