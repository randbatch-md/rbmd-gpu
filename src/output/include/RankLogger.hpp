#pragma once
#include <spdlog/spdlog.h>

#include <memory>
#include <utility>

class RankLogger {
 protected:
  int _current_rank;
  std::shared_ptr<spdlog::logger> _handle_logger;

 public:
  explicit RankLogger(int current_rank,
                      std::shared_ptr<spdlog::logger> spd_logger);

  template <typename... Args>
  void log(spdlog::level::level_enum level, const char *fmt,
           const Args &...args);

  template <typename... Args>
  void trace(const char *fmt, const Args &...args);

  template <typename... Args>
  void debug(const char *fmt, const Args &...args);

  template <typename... Args>
  void info(const char *fmt, const Args &...args);

  template <typename... Args>
  void warn(const char *fmt, const Args &...args);

  template <typename... Args>
  void error(const char *fmt, const Args &...args);

  template <typename... Args>
  void critical(const char *fmt, const Args &...args);
};

inline RankLogger::RankLogger(const int current_rank,
                              std::shared_ptr<spdlog::logger> spd_logger) {
  _current_rank = current_rank;
  _handle_logger = std::move(spd_logger);
  _handle_logger->set_pattern("%^[%l rank: " + std::to_string(_current_rank) +
                              "]%$ %v");
}

template <typename... Args>
void RankLogger::log(spdlog::level::level_enum level, const char *fmt,
                     const Args &...args) {
  _handle_logger->log(level, fmt, args...);
}

template <typename... Args>
void RankLogger::trace(const char *fmt, const Args &...args) {
  log(spdlog::level::trace, fmt, args...);
}

template <typename... Args>
void RankLogger::debug(const char *fmt, const Args &...args) {
  log(spdlog::level::debug, fmt, args...);
}

template <typename... Args>
void RankLogger::info(const char *fmt, const Args &...args) {
  log(spdlog::level::info, fmt, args...);
}

template <typename... Args>
void RankLogger::warn(const char *fmt, const Args &...args) {
  log(spdlog::level::warn, fmt, args...);
}

template <typename... Args>
void RankLogger::error(const char *fmt, const Args &...args) {
  log(spdlog::level::err, fmt, args...);
}

template <typename... Args>
void RankLogger::critical(const char *fmt, const Args &...args) {
  log(spdlog::level::critical, fmt, args...);
}