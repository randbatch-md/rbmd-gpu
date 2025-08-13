#include "TrajectoryOutput.h"
#include "memory_utils.h"
#include <fstream>
#include <algorithm>
#include <stdexcept>
#include "common/rbmd_define.h"
#include "Logger.hpp"


TrajectoryOutput::TrajectoryOutput(double available_memory_usage_ratio)
  : _memory_usage_ratio(available_memory_usage_ratio), _output_thread(), _mutex(), _cv_not_full(), _cv_not_empty(), _stream(nullptr) {}

TrajectoryOutput::~TrajectoryOutput() {
  _stop_flag = true;
  _cv_not_empty.notify_all();
  _cv_not_full.notify_all();
  if (_output_thread.joinable()) {
    _output_thread.join();
  }
  DeallocateRingBuffer();
  if (_stream) {
    CHECK_RUNTIME(STREAM_DESTORY(_stream));
  }
}

void TrajectoryOutput::Init() {
  _num_atoms = *(_structure_info_data->_num_atoms);
  _interval =  DataManager::getInstance().getConfigData()->Get<rbmd::Id>
  ("interval", "outputs", "trajectory_out");

  if (_num_atoms == 0) {
    Logger::Instance().warn(
        "TrajectoryOutput::Init called with 0 atoms. Output will be disabled.");
    return;
  }

  size_t frame_mem = CalculateSingleFrameMemory(_num_atoms);
  size_t available_mem = rbmd::utils::get_available_memory();
  size_t memory_to_use = static_cast<size_t>(
    available_mem * _memory_usage_ratio);
  size_t calculated_size = (memory_to_use > frame_mem)
                             ? (memory_to_use / frame_mem)
                             : 1;

  _ring_buffer_size = std::max(static_cast<size_t>(2),
                               std::min(calculated_size,
                                        static_cast<size_t>(64)));

  Logger::Instance().info("Trajectory Ring Buffer Initializing...");
  Logger::Instance().info("  Available Memory: {:.2f} GB",
                          available_mem / (1024.0 * 1024.0 * 1024.0));
  Logger::Instance().info("  Single Frame Memory: {:.2f} MB",
                          frame_mem / (1024.0 * 1024.0));
  Logger::Instance().info("  Calculated Ring Buffer Size: {}",
                          _ring_buffer_size);

  CHECK_RUNTIME(STREAM_CREATE(&_stream));
  AllocateRingBuffer(_num_atoms);
  _output_thread = std::thread(&TrajectoryOutput::OutputWorker, this);
}

void TrajectoryOutput::AllocateRingBuffer(size_t num_atoms) {
  _ring_buffer.resize(_ring_buffer_size);
  for (auto& frame : _ring_buffer) {
    frame.num_atoms = num_atoms;
    CHECK_RUNTIME(
        MALLOCHOST(reinterpret_cast<void**>(&frame.h_px), num_atoms * sizeof(
          rbmd::Real)));
    CHECK_RUNTIME(
        MALLOCHOST(reinterpret_cast<void**>(&frame.h_py), num_atoms * sizeof(
          rbmd::Real)));
    CHECK_RUNTIME(
        MALLOCHOST(reinterpret_cast<void**>(&frame.h_pz), num_atoms * sizeof(
          rbmd::Real)));
    CHECK_RUNTIME(
        MALLOCHOST(reinterpret_cast<void**>(&frame.h_vx), num_atoms * sizeof(
          rbmd::Real)));
    CHECK_RUNTIME(
        MALLOCHOST(reinterpret_cast<void**>(&frame.h_vy), num_atoms * sizeof(
          rbmd::Real)));
    CHECK_RUNTIME(
        MALLOCHOST(reinterpret_cast<void**>(&frame.h_vz), num_atoms * sizeof(
          rbmd::Real)));
    CHECK_RUNTIME(
        MALLOCHOST(reinterpret_cast<void**>(&frame.h_atoms_type), num_atoms *
          sizeof(rbmd::Id)));
    CHECK_RUNTIME(NEW_FLAG_EVENT(&frame.copy_complete_event, EVENT_DISABLE));
  }
}

void TrajectoryOutput::DeallocateRingBuffer() {
  for (auto& frame : _ring_buffer) {
    CHECK_RUNTIME(FREE(frame.h_px));
    CHECK_RUNTIME(FREE(frame.h_py));
    CHECK_RUNTIME(FREE(frame.h_pz));
    CHECK_RUNTIME(FREE(frame.h_vx));
    CHECK_RUNTIME(FREE(frame.h_vy));
    CHECK_RUNTIME(FREE(frame.h_vz));
    CHECK_RUNTIME(FREE(frame.h_atoms_type));
    CHECK_RUNTIME(EVENT_DESTORY(frame.copy_complete_event));
  }
  _ring_buffer.clear();
}

size_t TrajectoryOutput::CalculateSingleFrameMemory(size_t num_atoms) const {
  size_t mem = 0;
  mem += num_atoms * sizeof(rbmd::Real) * 3; // px, py, pz
  mem += num_atoms * sizeof(rbmd::Real) * 3; // vx, vy, vz
  mem += num_atoms * sizeof(rbmd::Id); // atoms_type
  return mem;
}

void TrajectoryOutput::Execute(int current_timestep) {
  if (current_timestep % _interval != 0 || _ring_buffer_size == 0)
    return;

  std::unique_lock<std::mutex> lock(_mutex);
  _cv_not_full.wait(lock, [this] {
    size_t next_write = (_write_index.load() + 1) % _ring_buffer_size;
    return next_write != _read_index.load();
  });

  size_t write_idx = _write_index.load();
  auto& target_frame = _ring_buffer[write_idx];
  target_frame.timestep = current_timestep;
  target_frame.box_snapshot = *(DataManager::getInstance().getMDData()->_box);

  // Release the lock before the asynchronous copy
  lock.unlock();

  // Asynchronous copy: Device to Host (Pinned Memory)
  const size_t pos_bytes = _num_atoms * sizeof(rbmd::Real);
  const size_t type_bytes = _num_atoms * sizeof(rbmd::Id);
  CHECK_RUNTIME(
      MEMCPY_ASYNC(target_frame.h_px, raw_ptr(_device_data->_d_px), pos_bytes,
        D2H, _stream));
  CHECK_RUNTIME(
      MEMCPY_ASYNC(target_frame.h_py, raw_ptr(_device_data->_d_py), pos_bytes,
        D2H, _stream));
  CHECK_RUNTIME(
      MEMCPY_ASYNC(target_frame.h_pz, raw_ptr(_device_data->_d_pz), pos_bytes,
        D2H, _stream));
  CHECK_RUNTIME(
      MEMCPY_ASYNC(target_frame.h_vx, raw_ptr(_device_data->_d_vx), pos_bytes,
        D2H, _stream));
  CHECK_RUNTIME(
      MEMCPY_ASYNC(target_frame.h_vy, raw_ptr(_device_data->_d_vy), pos_bytes,
        D2H, _stream));
  CHECK_RUNTIME(
      MEMCPY_ASYNC(target_frame.h_vz, raw_ptr(_device_data->_d_vz), pos_bytes,
        D2H, _stream));
  CHECK_RUNTIME(
      MEMCPY_ASYNC(target_frame.h_atoms_type, raw_ptr(_device_data->
        _d_atoms_type), type_bytes, D2H, _stream));

  // Record an event in the stream
  CHECK_RUNTIME(EVENT_RECORD(target_frame.copy_complete_event, _stream));

  lock.lock(); // Reacquire the lock to update the index
  _write_index.store((write_idx + 1) % _ring_buffer_size);
  _cv_not_empty.notify_one();
}

void TrajectoryOutput::OutputWorker() {
  std::ofstream trj_file("rbmd.trj", std::ios::out | std::ios::trunc);
  if (!trj_file.is_open()) {
    Logger::Instance().error("Failed to open trajectory file: rbmd.trj");
    return;
  }

  while (true) {
    std::unique_lock<std::mutex> lock(_mutex);
    _cv_not_empty.wait(lock, [this] {
      return (_read_index.load() != _write_index.load()) || _stop_flag.load();
    });

    if (_stop_flag.load() && (_read_index.load() == _write_index.load())) {
      break;
    }

    size_t read_idx = _read_index.load();
    const auto& source_frame = _ring_buffer[read_idx];
    lock.unlock();

    // Wait for the asynchronous copy of this frame to complete
    CHECK_RUNTIME(EVENT_SYNC(source_frame.copy_complete_event));

    // Write the frame data to a file.
    trj_file << "ITEM: TIMESTEP\n" << source_frame.timestep << "\n";
    trj_file << "ITEM: NUMBER OF ATOMS\n" << source_frame.num_atoms << "\n";
    trj_file << "ITEM: BOX BOUNDS pp pp pp\n";
    trj_file << source_frame.box_snapshot._coord_min[0] << " " << source_frame.box_snapshot._coord_max[0] << "\n";
    trj_file << source_frame.box_snapshot._coord_min[1] << " " << source_frame.box_snapshot._coord_max[1] << "\n";
    trj_file << source_frame.box_snapshot._coord_min[2] << " " << source_frame.box_snapshot._coord_max[2] << "\n";
    trj_file << "ITEM: ATOMS id type x y z vx vy vz\n";

    for (size_t i = 0; i < source_frame.num_atoms; ++i) {
      trj_file << i + 1 << " "
          << source_frame.h_atoms_type[i] << " "
          << source_frame.h_px[i] << " " << source_frame.h_py[i] << " " <<
          source_frame.h_pz[i] << " "
          << source_frame.h_vx[i] << " " << source_frame.h_vy[i] << " " <<
          source_frame.h_vz[i] << "\n";
    }

    lock.lock();
    _read_index.store((read_idx + 1) % _ring_buffer_size);
    _cv_not_full.notify_one();
  }
}