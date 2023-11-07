#ifndef DATA_GEN
#define DATA_GEN

#include "string"
#include <functional>

namespace dataGeneration {

const bool verbose = false;
enum DataDistribution {
  normal,
  normal_sorted,
  normal_reverse_sorted,
  normal_almost_sorted,
  uniform,
  uniform_sorted,
  uniform_reverse_sorted,
  uniform_almost_sorted,
};

static const DataDistribution All[] = {
    normal,  normal_sorted,  normal_reverse_sorted,  normal_almost_sorted,
    uniform, uniform_sorted, uniform_reverse_sorted, uniform_almost_sorted,
};

inline std::string getDataDistributionName(DataDistribution distribution) {
  std::string name;

  switch (distribution) {
  case DataDistribution::normal:
    name = "normal";
    break;
  case DataDistribution::normal_sorted:
    name = "normal-sorted";
    break;
  case DataDistribution::normal_reverse_sorted:
    name = "normal-reverse-sorted";
    break;
  case DataDistribution::normal_almost_sorted:
    name = "normal-almost-sorted";
    break;
  case DataDistribution::uniform:
    name = "uniform";
    break;
  case DataDistribution::uniform_sorted:
    name = "uniform-sorted";
    break;
  case DataDistribution::uniform_reverse_sorted:
    name = "uniform-reverse-sorted";
    break;
  case DataDistribution::uniform_almost_sorted:
    name = "uniform-almost-sorted";
    break;
  default:
    name = "unknown";
    break;
  }

  return name;
}

template <typename T>
using GeneratorFunc = std::function<void(T *, size_t, int, size_t)>;

const int DEFAULT_MAX_INT_VAL = 1 << 27;
// const int DEFAULT_MAX_INT_VAL = 1 << 18;
const int DEFAULT_MIN_INT_VAL = 0;

const float DEFAULT_MAX_FLOAT_VAL = 1 << 27;
const float DEFAULT_MIN_FLOAT_VAL = 0.0f;

const float DEFAULT_SWAP_PROBA = 0.05f;

const double DEFAULT_MAX_DOUBLE_VAL = 1 << 27;
const double DEFAULT_MIN_DOUBLE_VAL = 0.0f;

const long DEFAULT_MIN_LONG_VAL = 0;
const long DEFAULT_MAX_LONG_VAL = 1L << 54;

template <typename T> std::pair<T, T> get_default_limits();

template <typename T>
void generate_uniform_array(T *A, size_t n, int seed, size_t discard);

template <typename T>
void generate_uniform_sorted_array(T *A, size_t n, int seed, size_t discard,
                                   bool reverse);

template <typename T>
void generate_uniform_almost_sorted_array(T *A, size_t n, int seed,
                                          size_t discard, bool reverse,
                                          float swap_proba);

template <typename T>
void generate_uniform_array(T *A, size_t n, int seed, size_t discard, T min_val,
                            T max_val);

template <typename T>
void generate_normal_array(T *A, size_t n, int seed, size_t discard);

template <typename T>
void generate_normal_sorted_array(T *A, size_t n, int seed, size_t discard,
                                  bool reverse);

template <typename T>
void generate_normal_almost_sorted_array(T *A, size_t n, int seed,
                                         size_t discard, bool reverse,
                                         float swap_proba);
} // namespace dataGeneration

#endif /* BASIC_DRF_H_ */
