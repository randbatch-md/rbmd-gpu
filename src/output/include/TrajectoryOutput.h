#pragma once
#include <vector>
#include <string>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <atomic>
#include "output.h"
#include "rbmd_define.h"

// A wrapper for device pointers to clarify intent
struct DeviceDataPointers {
  const rbmd::Real* d_px{};
  const rbmd::Real* d_py{};
  const rbmd::Real* d_pz{};
  const rbmd::Real* d_vx{};
  const rbmd::Real* d_vy{};
  const rbmd::Real* d_vz{};
  const rbmd::Id* d_atoms_type{};
};

// Internal frame data, stored as a raw pointer to fixed memory
struct TrajectoryFrame {
  rbmd::Real* h_px{};
  rbmd::Real* h_py{};
  rbmd::Real* h_pz{};
  rbmd::Real* h_vx{};
  rbmd::Real* h_vy{};
  rbmd::Real* h_vz{};
  rbmd::Id* h_atoms_type{};
  rbmd::Real* h_charge{};
  rbmd::Id* h_atom_id_to_idx{};
  int timestep{};
  size_t num_atoms{};
  Box box_snapshot;
  EVENT_T copy_complete_event{};
};

class TrajectoryOutput:public Output {
public:
  // Constructor parameter: the desired ratio of available memory to use
  explicit TrajectoryOutput(double available_memory_usage_ratio = 0.25);
  ~TrajectoryOutput() override;
  // Disallow copying
  TrajectoryOutput(const TrajectoryOutput&) = delete;
  TrajectoryOutput& operator=(const TrajectoryOutput&) = delete;
  TrajectoryOutput(TrajectoryOutput&&) = delete;
  TrajectoryOutput& operator=(TrajectoryOutput&&) = delete;

  // Initializes the object: calculates buffer size, allocates resources, and starts the I/O thread
  void Init() override;

  // Main execution function, called by the main simulation loop.
  void Execute(int current_timestep) override;

private:
  void OutputWorker();
  void AllocateRingBuffer(size_t num_atoms);
  void DeallocateRingBuffer();
  [[nodiscard]]  size_t CalculateSingleFrameMemory(size_t num_atoms) const;
  [[nodiscard]]     bool _initialized = false;

  // A ring buffer.
  std::vector<TrajectoryFrame> _ring_buffer;
  size_t _ring_buffer_size{0};

  // Index for the blocking queue.
  std::atomic<size_t> _write_index{0};
  std::atomic<size_t> _read_index{0};

  // Synchronization primitive.
  std::thread _output_thread;
  std::mutex _mutex;
  std::condition_variable _cv_not_full;
  std::condition_variable _cv_not_empty;
  std::atomic<bool> _stop_flag{false};

  STREAM _stream;

  // Configuration.
  double _memory_usage_ratio;
  rbmd::Id _interval{};
  size_t _num_atoms{0};
};