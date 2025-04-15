#pragma once
#include <cmath>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include "../common/object.h"
#include "../common/types.h"
#include "common/json.hpp"

class ConfigData : public Object {
 public:
  /**
   * @brief constructor
   * @param file json config file
   */
  ConfigData(const std::string& file) {
    if (!IsJsonFile(file)) {
      //_console->error("{} is not json file!", file);
      return;
    }

    ParseJsonFile(file);
  }

  ~ConfigData() = default;

 public:
  /**
   * @brief get value
   * @tparam T
   * @tparam ...Args
   * @param key
   * @param ...args
   * @return value
   */
  template <typename T, typename... Args>
  T Get(std::string key, Args&&... args) {
    nlohmann::ordered_json json_node = _json_node;

    auto getNode = [this, &json_node](const auto& arg) {
      if (json_node.contains(arg) && json_node[arg].is_object()) {
        json_node = json_node[arg];
      } else {
        //_console->error("{} is not a object!", arg);
        return;
      }
    };

    (getNode(std::forward<Args>(args)), ...);

    try {
      if (json_node.contains(key)) {
        return json_node[key].get<T>();
      } else {
        throw std::runtime_error("no key named: " + key);
      }
    } catch (const std::exception&) {
      // log
      //_console->error("no key named: {}", key);
      throw;
    }
  }

  template <typename T, typename... Args>
  std::vector<T> GetArray(std::string key, Args&&... args) {
    nlohmann::ordered_json json_node = _json_node;

    auto getNode = [this, &json_node](const auto& arg) {
      if (json_node.contains(arg) && json_node[arg].is_object()) {
        json_node = json_node[arg];
      } else {
        //_console->error("{} is not a object!", arg);
        return;
      }
    };

    (getNode(std::forward<Args>(args)), ...);
    std::cout << "Checking key: " << key << " in node: " << json_node.dump(4)
              << std::endl;

    if (json_node.contains(key)) {
      nlohmann::ordered_json value = json_node[key];
      if (value.is_array()) {
        std::vector<T> result;
        for (const auto& item : value) {
          result.push_back(item.get<T>());  // 将数组元素转换为 T 类型
        }
        return result;
      } else {
        throw std::runtime_error(key + " is not an array");
      }
    } else {
      throw std::runtime_error("no key named: " + key);
    }
  }

  /**
   * @brief get json node
   * @param key
   * @return json node
   */
  auto& GetJsonNode(const std::string& key) {
    if (!_json_node.contains(key) || !_json_node[key].is_object()) {
      //_console->warn("Can not find key: {}", key);
    }

    return _json_node[key];
  }

  /**
   * @brief Check whether there is key Node
   * @param key node name
   * @return true or false
   */
  bool HasNode(const std::string& key) { return _json_node.contains(key); }

 private:
  /**
   * @brief Check whether the file is json file
   * @param file file path
   * @return true or false
   */
  bool IsJsonFile(const std::string& file) {
    auto length = file.length();
    return (length >= 5 && file.substr(length - 5) == ".json");
  }

  /**
   * @brief parse json file
   * @param file
   */
  void ParseJsonFile(const std::string& file) {
    try {
      std::ifstream filestream(file);
      if (!filestream.is_open()) {
        //_console->error("failed to open file: {}", file);
        return;
      }

      // 读取整个文件内容
      std::string json_str((std::istreambuf_iterator<char>(filestream)),
                           std::istreambuf_iterator<char>());
      filestream.close();

      // 解析JSON，允许注释
      // 最后一个参数true表示允许解析包含注释的JSON
      _json_node = nlohmann::json::parse(json_str, nullptr, false, true);
    } catch (const nlohmann::json::parse_error& e) {
      //_console->error("Error parsing JSON: {}", e.what());
    } catch (const std::exception& e) {
      //_console->error("Error parsing JSON: {}", e.what());
      }
    }

private:
    nlohmann::ordered_json _json_node;
};
