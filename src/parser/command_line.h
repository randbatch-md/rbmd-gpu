#pragma once

#include "common/object.h"
#include "cxxopts.hpp"
#include "version.h"

class CommandLine : public Object {
 public:
  CommandLine(int argc, char* argv[]);

  //
  CommandLine() = delete;
  CommandLine(const CommandLine&) = delete;
  CommandLine& operator=(const CommandLine&) = delete;

  virtual ~CommandLine() = default;

 public:
  bool RunApplication();
  static void Initialize();

  //
  std::string GetConfigPath() const  {
    if (!_co.count("j"))
      throw std::runtime_error("Usage: " " -j input.json");
    return _co["j"].as<std::string>();
  }

  //Help
  static std::string GetHelpText() {
    return R"(
    Usage:  -j input.json [options]

    Options:
      -j, --json FILE    Specify configuration file (required)
      -h, --help         Show this help message
      -v, --version      Show version information
    )";
  }
  //Version
  static std::string GetVersionText() {
    return VERSION;
  }

 private:
  void ParseCommand();

  static CommandLine& getInstance() {
    static CommandLine instance(0, nullptr);
    return instance;
  }

 private:
    cxxopts::ParseResult _co;
    cxxopts::Options _opts;
};