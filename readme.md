# Parallel Block-InsertionSort

## Original Work

https://doi.org/10.1016/j.jocs.2022.101866

Héctor Ferrada,
A sorting algorithm based on ordered block insertions,
Journal of Computational Science,

### PBIS - Titulo pendiente lol
sorting on multicore and manycore, an OpenMP-based parallel
implementation of Block-InsertionSort, etc
jvasquez INFO298 

## Requirements

All tests were done using g++ (GCC) 13.1.1 

In order to run in parallel the program must be compiled with the `-fopenmp`
flag, othewise, a sequential Block-InsertionSort algorithm will run. To compare
it with a parallel sort from the standard library you can either enable
parallel mode with `-D_GLIBCXX_PARALLEL` or include these flags when compiling:
`-ltbb -D_USETBB` to be able to use `std::execution::par` as execution policy.


Please refer to the official sources for installation instructions of oneAPI Threading Building Blocks (oneTBB)

- [OneTBB install](https://github.com/oneapi-src/oneTBB/blob/master/INSTALL.md) 

- [Installing on linux (Intel)](https://www.intel.com/content/www/us/en/docs/onetbb/get-started-guide/2021-9/overview.html) 


## Usage

### Sample code

```cpp
#include "pbis.hpp"
#include <vector>

int main(){
  int n = 1 << 18;
  int *A = new int[n];
  // fill your array with data

  //sequential
  bis::sort(A, A + n);

  //parallel
  bis::parallel::sort(A, A + n);

  // custom Comparison function
  bis::parallel::sort(A, A + n, std::less<int>());

  // vector
  std::vector<float> V(n);
  bis::parallel::sort(V.begin(), V.end())
}
```

### Compiling and running example `main.cpp`
```shell
foo@bar:~$ g++ -Wall -std=c++17 -O3 -fopenmp -I./pbis/include main.cpp -o prog
foo@bar:~$ ./prog
TIME  >>>                 std::sort     0.65004909 [s]
TIME  >>>                 bis::sort     0.68824921 [s]
Sort result Ok
TIME  >>>                 std::sort     0.64751839 [s]
TIME  >>>       bis::parallel::sort     0.28056250 [s]
Sort result Ok
```
