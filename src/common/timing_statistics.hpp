#include <chrono>
#include <iomanip>
#include <string>
#include <unordered_map>
#include <cmath>
#include <limits>
#include <vector>
#include <sstream>
#include <spdlog/spdlog.h>

class TimingStatistics {
private:
    struct Stats {
        double min = std::numeric_limits<double>::max();
        double max = std::numeric_limits<double>::lowest();
        double sum = 0.0;
        double sum_sq = 0.0;  // For variance calculation
        size_t count = 0;
        
        void update(double value) {
            min = std::min(min, value);
            max = std::max(max, value);
            sum += value;
            sum_sq += value * value;
            count++;
        }
        
        double avg() const {
            return count > 0 ? sum / count : 0.0;
        }
        
        double variance() const {
            if (count <= 1) return 0.0;
            double mean = avg();
            return (sum_sq / count) - (mean * mean);
        }
        
        double stddev() const {
            return std::sqrt(variance());
        }
    };
    
    std::unordered_map<std::string, Stats> categories;
    double total_time = 0.0;
    
    // Private constructor for singleton
    TimingStatistics() = default;
    
public:
    // Delete copy constructor and assignment operator
    TimingStatistics(const TimingStatistics&) = delete;
    TimingStatistics& operator=(const TimingStatistics&) = delete;
    
    // Get singleton instance
    static TimingStatistics& Instance() {
        static TimingStatistics instance;
        return instance;
    }
    
    // Record time for a category
    void record(const std::string& category, double time_seconds) {
        categories[category].update(time_seconds);
        total_time += time_seconds;
    }
    
  // Print statistics summary
  void print_summary() const {
      // Calculate average total time across all categories
      double avg_total = 0.0;
      for (const auto& [name, stats] : categories) {
        avg_total += stats.avg();
      }

      // Build output string
      std::stringstream ss;
      ss << "\n";
      ss << "┌──────────────────┬──────────────────┐\n";
      ss << "│     Category     │      Avg Time    │\n";
      ss << "├──────────────────┼──────────────────┤\n";

      // Set output format
      ss << std::fixed << std::setprecision(6);

      // Define output order
      const std::vector<std::string> order = {
        "Neighbor-List", "Short-Range", "Long-Range", "Bond","Angle","Dihedral","Improper"
    };

      for (const auto& section : order) {
        auto it = categories.find(section);
        if (it != categories.end()) {
          const Stats& stats = it->second;

          ss << "│ " << std::setw(16) << std::left << section << " │ "
               << std::setw(16) << std::right << stats.avg() << " │\n";
        }
      }

      ss << "└──────────────────┴──────────────────┘\n";
      // Log using spdlog
      Logger::Instance().info(ss.str().c_str());
    }
    
    // Reset statistics
    void reset() {
        categories.clear();
        total_time = 0.0;
    }
};

// Convenience macro for timing code blocks
#define TIME_SECTION(section_name, code_block) \
    do { \
        auto start = std::chrono::high_resolution_clock::now(); \
        code_block \
        auto end = std::chrono::high_resolution_clock::now(); \
        std::chrono::duration<double> duration = end - start; \
        TimingStatistics::getInstance().record(section_name, duration.count()); \
    } while(0)

