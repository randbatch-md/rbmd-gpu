#pragma once
#include "spdlog/spdlog.h"
#include <memory>
#include <utility>

class Logger {
 private:
  int _current_rank = -1;  // 默认值，表示未配置或无效
  std::shared_ptr<spdlog::logger> _handle_logger = nullptr;
  bool _is_configured = false;

  Logger() = default;
  ~Logger() = default;

 public:
  Logger(const Logger &) = delete;
  Logger &operator=(const Logger &) = delete;

  static Logger &Instance() {
    static Logger instance;
    return instance;
  }

  void Configure(int current_rank, std::shared_ptr<spdlog::logger> spd_logger) {
    if (_is_configured) {
      spdlog::warn(
          "Logger is already configured. Re-configuration attempt "
          "ignored.");
      return;
    }

    _current_rank = current_rank;
    _handle_logger = std::move(spd_logger);

    if (_handle_logger) {
      _is_configured = true;
    } else {
      spdlog::error(
          "Logger configuration failed: provided spd_logger is null.");
    }
  }

  bool IsConfigured() const { return _is_configured; }

  template <typename... Args>
  void log(spdlog::level::level_enum level, const char *fmt,
           const Args &...args) {
    if (_is_configured && _handle_logger && _current_rank == 0) {
      _handle_logger->log(level, fmt, args...);
    } else if (!_is_configured) {
      if (spdlog::default_logger_raw()) {
        spdlog::default_logger_raw()->warn(
            "Logger used before being configured. Original message: {}",
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
