#pragma once
#include <spdlog/spdlog.h>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>
#include <unordered_map>
#include <stdexcept>
#include <algorithm>

#include "types.h"
#include "output/include/Logger.hpp"

extern int test_current_step;

class ThermoStats {
private:
  bool columns_determined_and_header_printed = false;
  std::vector<std::string> active_keys;

  const std::vector<std::string> canonical_key_order = {
      "temperature", "pressure", "vdwl", "coul", "kspace",
    "bond", "angle", "dihedral", "improper" , "total-potential-energy"};

  const std::unordered_map<std::string, std::string> key_to_display_name = {
      {"temperature", "temperature"},
      {"pressure", "pressure"},
      {"vdwl", "E_vdwl"},
      {"coul", "E_coul"},
     {"kspace", "E_kspace"},
      {"bond", "E_bond"},
      {"angle", "E_angle"},
      {"dihedral", "E_dihedral"},
      {"improper", "E_improper"},
      {"total-potential-energy", "E_pe"}
  };

  const std::unordered_map<std::string, int> key_to_width = {
      {"temperature", 13},
      {"pressure",    12},
      {"vdwl",        15},
      {"coul",        14},
     {"kspace",        14},
      {"bond",        12},
      {"angle",       12},
      {"dihedral",     12},
      {"improper",    12},
      {"total-potential-energy",  12}
  };

  const std::unordered_map<std::string, int> key_to_precision = {
      {"temperature", 4},
      {"pressure", 4},
      {"vdwl", 4},
      {"coul", 4},
     {"kspace", 4},
      {"bond", 4},
      {"angle", 4},
      {"dihedral", 4},
      {"improper", 4},
      {"total-potential-energy", 4}
  };

  std::unordered_map<std::string, rbmd::Real> thermo_data;
  int current_step = -1;

  ThermoStats() = default;

  void PrintHeader() {
    if (active_keys.empty()) {
      return;
    }

    std::stringstream header_ss;
    header_ss << std::setw(8) << "Step";

    for (const auto& key : active_keys) {
      std::string display_name = key_to_display_name.count(key)
                                   ? key_to_display_name.at(key)
                                   : key;  //
      // MODIFIED: Reduced default width for other keys
      int width = key_to_width.count(key) ? key_to_width.at(key) : 15;
      header_ss << std::setw(width) << display_name;
    }

    Logger::Instance().info("{}", header_ss.str());
  }

  void DetermineAndPrintHeaderIfNeeded() {
    if (columns_determined_and_header_printed || thermo_data.empty()) {
      return;
    }

    // 1. Add existing keys according to the canonical order
    for (const auto& canonical_key : canonical_key_order) {
      if (thermo_data.count(canonical_key)) {
        if (std::find(active_keys.begin(), active_keys.end(), canonical_key) ==
            active_keys.end()) {
          active_keys.push_back(canonical_key);
        }
      }
    }

    // 2. Add other keys present in thermo_data but not in the canonical order list
    std::vector<std::string> additional_keys;
    for (const auto& pair : thermo_data) {
      const std::string& key = pair.first;
      if (std::find(active_keys.begin(), active_keys.end(), key) == active_keys.
          end()) {
        additional_keys.push_back(key);
      }
    }
    active_keys.insert(active_keys.end(), additional_keys.begin(),
                       additional_keys.end());

    if (!active_keys.empty()) {
      PrintHeader();
      columns_determined_and_header_printed = true;
    } else {
      Logger::Instance().warn(
          "ThermoStats: No keys found in the first data set (step {}). Header cannot be printed.",
          current_step);
    }
  }

public:
  static ThermoStats& Instance() {
    static ThermoStats instance;
    return instance;
  }

  ThermoStats(const ThermoStats&) = delete;

  ThermoStats& operator=(const ThermoStats&) = delete;

  void SetStep(int step) {
    current_step = step;
  }

  void AddThermoData(const std::string& key, rbmd::Real value) {
    thermo_data[key] = value;
  }

  bool ShouldOutput()
  {
    auto interval = DataManager::getInstance().getConfigData()->Get<rbmd::Id>(
"interval", "outputs", "thermo_out");
    return test_current_step % interval == 0;
  }

  void OutputRow() {
    if (ShouldOutput()) {
      if (!columns_determined_and_header_printed && !thermo_data.empty()) {
        DetermineAndPrintHeaderIfNeeded();
      }

      if (!columns_determined_and_header_printed || active_keys.empty()) {
        if (current_step != -1 && !thermo_data.empty()) {
          Logger::Instance().warn(
              "ThermoStats: Header not initialized, skipping data row for step {}.",
              current_step);
        } else if (current_step != -1 && thermo_data.empty() && !
                   columns_determined_and_header_printed) {
          Logger::Instance().debug(
              "ThermoStats: No data to determine header or print for step {}.",
              current_step);
                   }
        thermo_data.clear();
        return;
      }

      std::stringstream ss;
      ss << std::fixed;

      ss << std::setw(8) << current_step;

      for (const auto& key : active_keys) {
        // MODIFIED: Reduced default width for other keys
        int width = key_to_width.count(key) ? key_to_width.at(key) : 15;
        int precision = key_to_precision.count(key)
                          ? key_to_precision.at(key)
                          : 6; // Default precision for other keys

        if (thermo_data.count(key)) {
          ss << std::setw(width) << std::setprecision(precision) << thermo_data.
              at(key);
        } else {
          std::string na_str = "N/A";
          std::stringstream temperature_ss_na; // This stringstream is not strictly necessary here
          // could directly write to ss.
          temperature_ss_na << std::setw(width) << na_str;
          ss << temperature_ss_na.str();
        }
      }

      Logger::Instance().info("{}", ss.str());
      thermo_data.clear();
    }
  }

};