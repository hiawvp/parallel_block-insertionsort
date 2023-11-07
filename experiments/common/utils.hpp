
#ifndef UTILS_XD
#define UTILS_XD

#include <iomanip>
#include <ios>
#include <iostream>
#include <stdio.h>
#include <string>
#include <string_view>

#ifndef DEBUG_FLAG
#define DEBUG_FLAG 1 // set debug mode
#endif

#ifndef TIMING_FLAG
#define TIMING_FLAG 1 // set timing mode
#endif

#if DEBUG_FLAG
#define DEBUG(x)                                                               \
  do {                                                                         \
    std::cerr << "DEBUG >>> [" << __FUNCTION__ << "()][Line " << __LINE__      \
              << "] " << std::setw(20) << #x << ": " << x << std::endl;        \
  } while (0)
#else
#define DEBUG(x)
#endif

#if DEBUG_FLAG
#define DEBUG2(x)                                                              \
  do {                                                                         \
    std::cerr << "DEBUG >>> [" << __FUNCTION__ << "()][Line " << __LINE__      \
              << "] " << std::setw(20) << x << std::endl;                      \
  } while (0)
#else
#define DEBUG2(x)
#endif

#if TIMING_FLAG
#define timing(f)                                                              \
  do {                                                                         \
    auto start = std::chrono::high_resolution_clock::now();                    \
    f;                                                                         \
    auto finish = std::chrono::high_resolution_clock::now();                   \
    std::chrono::duration<double> elapsed = finish - start;                    \
    std::cout << "TIME  >>> " << std::setw(65) << std::string_view(#f, 50)     \
              << "...) " << std::setiosflags(std::ios::fixed)                  \
              << std::setprecision(8) << std::setw(20) << elapsed.count()      \
              << std::setw(4) << " [s]" << std::endl;                          \
  } while (0)
#else
#define timing(a)
#endif

#if TIMING_FLAG
#define timing2(f)                                                             \
  do {                                                                         \
    auto start = std::chrono::high_resolution_clock::now();                    \
    f;                                                                         \
    auto finish = std::chrono::high_resolution_clock::now();                   \
    std::chrono::duration<double> elapsed = finish - start;                    \
    std::cout << std::setprecision(8) << std::setw(20) << elapsed.count()      \
              << ",";                                                          \
  } while (0)
#else
#define timing(a)
#endif
namespace utils {

enum DataTypes { int_data, float_data, long_data, double_data, unknown };
const int MAX_PRINT = 256;
const ulong PRINT_LIMIT = 220;

void print_sep();

// prints elements in a single line
// void print_array(int *A, size_t n);
template <typename T>
void print_array(T *A, size_t n, std::string hint = "array:");

template <typename T> DataTypes get_datatype() {
  DataTypes data_type = unknown;
  if (std::is_same<T, int>::value) {
    data_type = int_data;
  } else if (std::is_same<T, long>::value) {
    data_type = long_data;
  } else if (std::is_same<T, float>::value) {
    data_type = float_data;
  } else if (std::is_same<T, double>::value) {
    data_type = double_data;
  }
  return data_type;
}

std::string get_typename(DataTypes val);

/*Returns false at first mismatch, true if all elements are equal. */
template <typename T> bool check_sort_result(T *std_sort, T *my_sort, size_t n);
size_t generate_random_number(int min, int max, int seed);
} // namespace utils

#endif /* BASIC_DRF_H_ */
