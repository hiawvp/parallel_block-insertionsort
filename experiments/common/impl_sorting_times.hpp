#ifndef IMPL_SORTING_TIMES_HPP
#define IMPL_SORTING_TIMES_HPP

// TODO: exportar pdfs con fecha en el nombre
// TODO: exportar pdfs y jsons en la carpeta results

#include "utils.hpp"
#include <cstddef>
#include <cstdio>
#include <filesystem>
#include <functional>
#include <string>
#include <utility>
#include <variant>
#include <vector>

template <typename T> class SortingImplementationDetails {
public:
  using arg_type = std::pair<std::string, int>;
  using arg_list = std::vector<arg_type>;
  using time_list = std::vector<std::pair<size_t, std::vector<double>>>;
  using fn_type = std::function<void(T *, size_t)>;

  void add_times(size_t n, const std::vector<double> &times) {
    times_.push_back({n, times});
  }

  void run(T *items, size_t n) { fn_(items, n); }

  std::string get_name() { return name_; }

  SortingImplementationDetails(const std::string &name,
                               const arg_list &extra_args, const fn_type &fn)
      : name_(name), extra_args_(extra_args), fn_(fn) {}

  std::string to_json() const {
    std::string json_str;

    json_str += "{\n";
    // Serialize experiment info
    // json_str += "  \"experiment_info\": " + info_.to_json() + ",\n";

    json_str += "  \"name\": \"" + name_ + "\",\n";
    // Serialize extra arguments
    //
    if (extra_args_.empty()) {
      json_str += "  \"extra_args\": [],\n";
    } else {
      json_str += "  \"extra_args\": [\n";
      for (size_t i = 0; i < extra_args_.size(); i++) {
        json_str += "    { \"name\": \"" + extra_args_[i].first +
                    "\", \"value\": " + std::to_string(extra_args_[i].second) +
                    " }";
        if (i < extra_args_.size() - 1) {
          json_str += ",\n";
        } else {
          json_str += "\n";
        }
      }
      json_str += "  ],\n";
    }
    // Serialize execution times
    json_str += "  \"exec_times\": [\n";
    for (size_t i = 0; i < times_.size(); i++) {
      json_str +=
          "    { \"n\": " + std::to_string(times_[i].first) + ", \"times\": [";
      for (size_t j = 0; j < times_[i].second.size(); j++) {
        json_str += std::to_string(times_[i].second[j]);
        if (j < times_[i].second.size() - 1) {
          json_str += ",";
        }
      }
      json_str += "] }";
      if (i < times_.size() - 1) {
        json_str += ",\n";
      } else {
        json_str += "\n";
      }
    }
    json_str += "  ]\n";

    json_str += "}";

    return json_str;
  }

private:
  std::string name_;
  arg_list extra_args_;
  fn_type fn_;
  time_list times_;
};

// TODO: add generator and seed
template <typename T> class Experiment {
public:
  using arg_type = std::pair<std::string, std::string>;
  using arg_list = std::vector<arg_type>;
  using Impls = std::vector<SortingImplementationDetails<T>>;
  // ke kek ekk e
  Experiment(const Impls &implementations, const std::string &input_type,
             const std::string &input_distribution, const arg_list &extra_args,
             const int &seed, const int &reps)
      : implementations_(implementations), input_type_(input_type),
        input_distribution_(input_distribution), extra_args_(extra_args),
        seed_(seed), reps_(reps) {}

  std::string to_json() const {
    std::string json_str;

    json_str += "{\n";
    json_str += "  \"experiment_info\": {\n";
    // json_str += "  \"n\": " + std::to_string(n_) + ",\n";
    json_str += "  \"input_type\": \"" + input_type_ + "\",\n";
    json_str += "  \"input_distribution\": \"" + input_distribution_ + "\",\n";
    if (extra_args_.empty()) {
      json_str += "  \"extra_args\": []\n";
    } else {
      json_str += "  \"extra_args\": [\n";
      for (size_t i = 0; i < extra_args_.size(); i++) {
        json_str += "    { \"name\": \"" + extra_args_[i].first +
                    "\", \"value\": \"" + extra_args_[i].second + "\" }";
        if (i < extra_args_.size() - 1) {
          json_str += ",\n";
        } else {
          json_str += "\n";
        }
      }
      json_str += "  ]\n";
    }
    json_str += "},\n";
    json_str += "  \"results\": ";
    json_str += "[\n ";
    DEBUG(implementations_.size());
    if (!implementations_.empty()) {
      std::string sep;
      for (size_t i = 0; i < implementations_.size(); ++i) {
        sep = (i == implementations_.size() - 1 ? "" : ",");
        json_str += implementations_[i].to_json() + sep + "\n";
      }
    }
    json_str += "  ]\n";
    json_str += "}";
    return json_str;
  }

  void dump_to_file(std::string filename = "") {
    DEBUG(filename);
    if (filename.size() == 0) {
      filename = "output_temp.json";
    }
    std::cout << "output: " << filename << std::endl;
    std::string json_str = to_json();
    std::string command =
        "echo '" + json_str + "' | python -m json.tool > " + filename;
    std::system(command.c_str());

    auto output_path =
        std::filesystem::current_path().string() + "/" + filename;
    dumped_results_file_path_ = output_path;
    DEBUG(output_path);
  }

  void generate_pdf_plot(std::string outputFilename, bool open = false) {
    DEBUG(outputFilename);
    bool temp_json_required = false;

    if (dumped_results_file_path_.size() == 0) {
      temp_json_required = true;
      std::string filename =
          "temp_output_" + std::to_string(std::time(nullptr)) + ".json";
      dump_to_file(filename);
    }
    DEBUG(dumped_results_file_path_);

    std::string title = "Execution time " + input_distribution_ + " " +
                        input_type_ + " data, " + std::to_string(reps_) +
                        " reps.";
    if (!extra_args_.empty()) {
      title += "\n (";
      for (size_t i = 0; i < extra_args_.size(); i++) {
        title += " " + extra_args_[i].first + "=" + extra_args_[i].second;
        title += i < extra_args_.size() - 1 ? ", " : "";
      }
      title += ")";
    }
    DEBUG(title);
    std::string command = "node ./plot_chartjs/index.js -i " +
                          dumped_results_file_path_ + " -o " + outputFilename +
                          " -t \"" + title + "\"";
    std::cout << std::system(command.c_str()) << std::endl;
    if (temp_json_required) {
      DEBUG("deleting dumped results file");
      std::cout << std::system(("rm " + dumped_results_file_path_).c_str())
                << std::endl;
    }
    if (open) {
      std::system(("evince " + outputFilename + ".pdf &").c_str());
    }
  }

  int get_seed() { return seed_; };
  int get_reps() { return reps_; };
  Impls &get_implementations() { return implementations_; };

private:
  // size_t n_;
  Impls implementations_;
  std::string input_type_;
  std::string input_distribution_;
  arg_list extra_args_;
  int seed_;
  int reps_;
  std::string generator_ = "std::mt19937";
  std::string dumped_results_file_path_ = "";
};

#endif // IMPL_SORTING_TIMES_HPP
