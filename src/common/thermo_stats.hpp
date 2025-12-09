#pragma once
#include <spdlog/spdlog.h>
#include <iomanip>
#include <memory>
#include <sstream>
#include <string>
#include <vector>
#include <unordered_map>
#include <algorithm>
#include <cmath>

#include "types.h"
#include "output/include/Logger.hpp"
#include "thermo_view.h"

extern int test_current_step;

class ThermoStats {
private:
  bool columns_determined_and_header_printed = false;
  std::vector<std::string> active_keys;
  int row_count = 0;

  // ============ 可配置参数 ============
  static constexpr int HEADER_REPEAT_INTERVAL = 101;
  static constexpr int COLUMN_WIDTH = 14;
  static constexpr int STEP_COLUMN_WIDTH = 12;      // Step列宽度
  static constexpr int SCIENTIFIC_PRECISION = 4;  // 科学计数法精度

  const std::vector<std::string> canonical_key_order = {
      "temperature", "pressure", "vdwl", "coul", "kspace",
      "bond", "angle", "dihedral", "improper", "total-potential-energy"};

  const std::unordered_map<std::string, std::string> key_to_display_name = {
      {"temperature", "Temp"},
      {"pressure", "Press"},
      {"vdwl", "E_vdwl"},
      {"coul", "E_coul"},
      {"kspace", "E_kspace"},
      {"bond", "E_bond"},
      {"angle", "E_angle"},
      {"dihedral", "E_dihed"},
      {"improper", "E_improp"},
      {"total-potential-energy", "PotEng"}
  };

  std::unordered_map<std::string, rbmd::Real> thermo_data;
  int current_step = -1;
  std::shared_ptr<IThermoView> thermo_view;

  ThermoStats() = default;

  // ============ 格式化工具函数 ============
   /**
   * @brief 格式化步数
   * - 小于1千万：普通整数
   * - 大于等于1千万：科学计数法
   */
  std::string FormatStep(int64_t step) const {
    std::stringstream ss;
    if (step < 10000000) {  // 小于1千万，直接显示整数
      ss << step;
    } else {  // 大步数用科学计数法
      ss << std::scientific << std::setprecision(2) << static_cast<double>(step);
    }
    return ss.str();
  }
  
  /**
   * @brief 全部使用科学计数法格式化
   */
  std::string FormatValue(rbmd::Real value) const {
    std::stringstream ss;
    ss << std::scientific << std::setprecision(SCIENTIFIC_PRECISION) << value;
    return ss.str();
  }

  /**
   * @brief 生成分隔线
   */
  std::string GenerateSeparator(char fill = '-') const {
    std::stringstream ss;
    ss << "+" << std::string(STEP_COLUMN_WIDTH, fill);  // Step列
    for (size_t i = 0; i < active_keys.size(); ++i) {
      ss << "+" << std::string(COLUMN_WIDTH, fill);     // 数据列
    }
    ss << "+";
    return ss.str();
  }

  /**
   * @brief 打印表头
   */
  void PrintHeader() {
    if (active_keys.empty()) return;

    std::string thick_sep = GenerateSeparator('=');
    Logger::Instance().info("{}", thick_sep);

    std::stringstream header_ss;
    header_ss << "|" << std::setw(STEP_COLUMN_WIDTH) << std::right << "Step";  // Step列
    
    for (const auto& key : active_keys) {
      std::string display_name = key_to_display_name.count(key)
                                     ? key_to_display_name.at(key)
                                     : key;
      if (display_name.length() > static_cast<size_t>(COLUMN_WIDTH - 1)) {
        display_name = display_name.substr(0, COLUMN_WIDTH - 2) + "~";
      }
      header_ss << "|" << std::setw(COLUMN_WIDTH) << std::right << display_name;
    }
    header_ss << "|";

    Logger::Instance().info("{}", header_ss.str());
    Logger::Instance().info("{}", thick_sep);
  }

  void DetermineAndPrintHeaderIfNeeded() {
    if (columns_determined_and_header_printed || thermo_data.empty()) {
      return;
    }

    for (const auto& canonical_key : canonical_key_order) {
      if (thermo_data.count(canonical_key)) {
        if (std::find(active_keys.begin(), active_keys.end(), canonical_key) ==
            active_keys.end()) {
          active_keys.push_back(canonical_key);
        }
      }
    }

    std::vector<std::string> additional_keys;
    for (const auto& pair : thermo_data) {
      const std::string& key = pair.first;
      if (std::find(active_keys.begin(), active_keys.end(), key) ==
          active_keys.end()) {
        additional_keys.push_back(key);
      }
    }
    std::sort(additional_keys.begin(), additional_keys.end());
    active_keys.insert(active_keys.end(), additional_keys.begin(),
                       additional_keys.end());

    if (!active_keys.empty()) {
      if (!thermo_view) {
        Logger::Instance().info("");
        PrintHeader();
      }
      columns_determined_and_header_printed = true;
    }
  }

  ThermoFrame BuildFrame() const {
    ThermoFrame frame;
    frame.step = current_step;
    frame.ordered_keys = active_keys;
    frame.values = thermo_data;
    for (const auto& key : active_keys) {
      auto it = key_to_display_name.find(key);
      frame.labels[key] = it != key_to_display_name.end() ? it->second : key;
    }
    return frame;
  }

  void PrintRow(const ThermoFrame& frame) {
    std::stringstream ss;
    ss << "|" << std::setw(STEP_COLUMN_WIDTH) << std::right << FormatStep(frame.step);

    for (const auto& key : frame.ordered_keys) {
      ss << "|";
      auto value_it = frame.values.find(key);
      if (value_it != frame.values.end()) {
        ss << std::setw(COLUMN_WIDTH) << std::right << FormatValue(value_it->second);
      } else {
        ss << std::setw(COLUMN_WIDTH) << std::right << "---";
      }
    }
    ss << "|";

    Logger::Instance().info("{}", ss.str());
  }

public:
  static ThermoStats& Instance() {
    static ThermoStats instance;
    return instance;
  }

  ThermoStats(const ThermoStats&) = delete;
  ThermoStats& operator=(const ThermoStats&) = delete;

  void SetView(std::shared_ptr<IThermoView> view) {
    thermo_view = std::move(view);
    if (thermo_view) {
      thermo_view->OnReset();
    }
  }

  bool HasView() const { return static_cast<bool>(thermo_view); }

  void SetStep(int step) { current_step = step; }

  void AddThermoData(const std::string& key, rbmd::Real value) {
    thermo_data[key] = value;
  }

  bool ShouldOutput() {
    auto interval = DataManager::getInstance().getConfigData()->Get<rbmd::Id>(
        "interval", "outputs", "thermo_out");
    return test_current_step % interval == 0;
  }

  void OutputRow() {
    if (!ShouldOutput()) return;

    if (!columns_determined_and_header_printed && !thermo_data.empty()) {
      DetermineAndPrintHeaderIfNeeded();
    }

    if (!columns_determined_and_header_printed || active_keys.empty()) {
      thermo_data.clear();
      return;
    }

    if (row_count > 0 && row_count % HEADER_REPEAT_INTERVAL == 0) {
      if (!thermo_view) {
        Logger::Instance().info("{}", GenerateSeparator('-'));
        PrintHeader();
      }
    }

    ThermoFrame frame = BuildFrame();
    if (thermo_view) {
      thermo_view->OnFrame(frame);
    } else {
      PrintRow(frame);
    }

    row_count++;
    thermo_data.clear();
  }

  void Finalize() {
    if (thermo_view) {
      thermo_view->OnShutdown();
      return;
    }
    if (columns_determined_and_header_printed && !active_keys.empty()) {
      Logger::Instance().info("{}", GenerateSeparator('='));
    }
  }

  void Reset() {
    columns_determined_and_header_printed = false;
    active_keys.clear();
    thermo_data.clear();
    row_count = 0;
    current_step = -1;
    if (thermo_view) {
      thermo_view->OnReset();
    }
  }
};
