#include <spdlog/spdlog.h>
#include <spdlog/sinks/base_sink.h>
#include <fcntl.h>
#include <unistd.h>
#include <sys/stat.h>
#include <errno.h>
#include <string.h>
#include <thread>
#include <chrono>

template<typename Mutex>
class linux_pipe_sink : public spdlog::sinks::base_sink<Mutex> {
public:
    linux_pipe_sink(const std::string& pipe_name)
        : pipe_name_(pipe_name), fd_(-1) {

        // Check if the pipe exists
        struct stat st;
        if (stat(pipe_name.c_str(), &st) == 0) {
            // If it exists, delete it
            if (unlink(pipe_name.c_str()) != 0) {
                throw spdlog::spdlog_ex("Failed to delete existing pipe: " + pipe_name + ", error: " + strerror(errno));
            }
        }

        // Create a new named pipe
        if (mkfifo(pipe_name.c_str(), 0666) != 0) {
            throw spdlog::spdlog_ex("Failed to create named pipe: " + pipe_name + ", error: " + strerror(errno));
        }

        //! Open the pipe in blocking + read-write mode (this is the key modification)
        //! Opening for both read and write simultaneously can avoid waiting for another process to open the other end
        fd_ = open(pipe_name.c_str(), O_RDWR);

        if (fd_ == -1) {
            std::string error_msg = "Failed to open named pipe: " + pipe_name + ", error: " + strerror(errno);
            // If open fails, try to remove the created pipe before throwing
            unlink(pipe_name_.c_str()); // Attempt to clean up
            throw spdlog::spdlog_ex(error_msg);
        }
    }

    ~linux_pipe_sink() {
        if (fd_ != -1) {
            close(fd_);
        }
        // Delete the pipe
        unlink(pipe_name_.c_str());
    }

protected:
    void sink_it_(const spdlog::details::log_msg& msg) override {
        if (fd_ == -1) {
            return; // Pipe is invalid, skip writing
        }

        spdlog::memory_buf_t formatted;
        this->formatter_->format(msg, formatted);

        // Write to the named pipe
        ssize_t written = write(fd_, formatted.data(), formatted.size());
        if (written < 0) {
            // Handle write error
            std::cerr << "Failed to write to pipe: " << strerror(errno) << std::endl;
        }
    }

    void flush_() override {
        // Named pipes do not require explicit flushing
    }

private:
    std::string pipe_name_;
    int fd_;
};

using linux_pipe_sink_mt = linux_pipe_sink<std::mutex>;
using linux_pipe_sink_st = linux_pipe_sink<spdlog::details::null_mutex>;