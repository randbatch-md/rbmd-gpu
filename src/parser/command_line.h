#pragma once

#include "common/object.h"
#include "cxxopts.hpp"

class CommandLine : public Object {
 public:
  CommandLine(int argc, char* argv[]);
  virtual ~CommandLine() = default;

 public:
  bool RunApplication();
  static void Initialize();

 private:
  void ParseCommand();

 private:
  cxxopts::ParseResult _co;
  static cxxopts::Options _opts;
};