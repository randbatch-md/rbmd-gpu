#pragma once
#include "spdlog/spdlog.h"

#include <memory>
#include <string>
#include <utility>

class RankLogger {
 private:
  int _current_rank = -1;
  std::shared_ptr<spdlog::logger> _handle_logger = nullptr;
  bool _is_configured = false;

  RankLogger() = default;
  ~RankLogger() = default;

 public:
  RankLogger(const RankLogger &) = delete;
  RankLogger &operator=(const RankLogger &) = delete;

  static RankLogger &getInstance() {
    static RankLogger instance;
    return instance;
  }

  void Configure(int current_rank, std::shared_ptr<spdlog::logger> spd_logger) {
    if (_is_configured) {
      spdlog::warn(
          "RankLogger is already configured for rank {}. Re-configuration "
          "attempt ignored.",
          _current_rank);
      return;
    }

    if (!spd_logger) {
      spdlog::error(
          "RankLogger configuration failed: provided spd_logger is null.");
      return;
    }

    _current_rank = current_rank;
    _handle_logger = std::move(spd_logger);

    _handle_logger->set_pattern("%^[%l rank: " + std::to_string(_current_rank) +
                                "]%$ %v");

    _is_configured = true;
  }

  bool isConfigured() const { return _is_configured; }

  template <typename... Args>
  void log(spdlog::level::level_enum level, const char *fmt,
           const Args &...args) {
    if (_is_configured && _handle_logger) {
      _handle_logger->log(level, fmt, args...);
    } else if (!_is_configured) {
      if (spdlog::default_logger_raw()) {
        spdlog::default_logger_raw()->warn(
            "RankLogger used before being configured. Original message: {}",
            fmt::format(fmt, args...));
      }
    }
  }

  template <typename... Args>
  void trace(const char *fmt, const Args &...args) {
    log(spdlog::level::trace, fmt, args...);
  }

  template <typename... Args>
  void debug(const char *fmt, const Args &...args) {
    log(spdlog::level::debug, fmt, args...);
  }

  template <typename... Args>
  void info(const char *fmt, const Args &...args) {
    log(spdlog::level::info, fmt, args...);
  }

  template <typename... Args>
  void warn(const char *fmt, const Args &...args) {
    log(spdlog::level::warn, fmt, args...);
  }

  template <typename... Args>
  void error(const char *fmt, const Args &...args) {
    log(spdlog::level::err, fmt, args...);
  }

  template <typename... Args>
  void critical(const char *fmt, const Args &...args) {
    log(spdlog::level::critical, fmt, args...);
  }
};