#pragma once
#include <iostream>
#include <memory>

#include "spdlog/sinks/stdout_color_sinks.h"
#include "spdlog/spdlog.h"

class Object {
 public:
  Object() {
    //_console = spdlog::get("console");
    // if (!_console) {
    //  _console = spdlog::stdout_color_mt("console");
    //}
  }

  virtual ~Object() = default;

 protected:
  // std::shared_ptr<spdlog::logger> _console;
};