#pragma once
#include <memory>
#include <utility>
#include <spdlog/spdlog.h>

class GlobalLogger {
protected:
    int _current_rank;
    std::shared_ptr<spdlog::logger> _handle_logger;

public:
    explicit GlobalLogger(int current_rank, std::shared_ptr<spdlog::logger> spd_logger);

    template<typename... Args>
    void log(spdlog::level::level_enum level, const char *fmt, const Args &... args);

    template<typename... Args>
    void trace(const char *fmt, const Args &... args);

    template<typename... Args>
    void debug(const char *fmt, const Args &... args);

    template<typename... Args>
    void info(const char *fmt, const Args &... args);

    template<typename... Args>
    void warn(const char *fmt, const Args &... args);

    template<typename... Args>
    void error(const char *fmt, const Args &... args);

    template<typename... Args>
    void critical(const char *fmt, const Args &... args);
};

inline GlobalLogger::GlobalLogger(const int current_rank, std::shared_ptr<spdlog::logger> spd_logger) {
    _current_rank = current_rank;
    _handle_logger = std::move(spd_logger);
}

template<typename... Args>
void GlobalLogger::log(spdlog::level::level_enum level, const char *fmt, const Args &... args) {
    if (_current_rank == 0) {
        _handle_logger->log(level, fmt, args...);
    }
}

template<typename... Args>
void GlobalLogger::trace(const char *fmt, const Args &... args) {
    log(spdlog::level::trace, fmt, args...);
}

template<typename... Args>
void GlobalLogger::debug(const char *fmt, const Args &... args) {
    log(spdlog::level::debug, fmt, args...);
}

template<typename... Args>
void GlobalLogger::info(const char *fmt, const Args &... args) {
    log(spdlog::level::info, fmt, args...);
}

template<typename... Args>
void GlobalLogger::warn(const char *fmt, const Args &... args) {
    log(spdlog::level::warn, fmt, args...);
}

template<typename... Args>
void GlobalLogger::error(const char *fmt, const Args &... args) {
    log(spdlog::level::err, fmt, args...);
}

template<typename... Args>
void GlobalLogger::critical(const char *fmt, const Args &... args) {
    log(spdlog::level::critical, fmt, args...);
}