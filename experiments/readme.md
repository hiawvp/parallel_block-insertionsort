# Experiments

this folder contains all the implementations used during development. The output names for some implementations may differ with the ones used in the paper

check the `ExpConfig` struct default values in `main.cpp` for seed, input sizes and more.

**Note: requires linking to [tbb](https://github.com/oneapi-src/oneTBB)**

## Random number generators used:

`./common/data_gen.cpp`

## Actual imports and setup of parallel implementations:

`./common/parallel_impls.cpp`

## Debug logging and other:

`./common/utils.hpp`

## Sequential Block-InsertionSort from the original work (Ferrada, H. 2022)

`./Block-InsertionSort/original_bis/`

# Usage

```shell
make
./prog <num_threads> <k> <reps>
```

this will generate, for each experiment, a json file that has the following schema:

```json
{
  "experiment_info": {
    "input_type": {
      "type": "string",
    },
    "input_distribution": "string",
    "extra_args": [
      {
        "name": "string",
        "value": "string"
      }
    ]
  },
  "results": [
    {
      "name": "string",
      "extra_args": [
        {
          "name": "string",
          "value": "number"
        },
      ],
      "exec_times": [
        {
          "n": "number",
          "times": [
            "number"
          ]
        }
      }
    }
  ]
}
```

