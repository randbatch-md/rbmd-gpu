#include "TrajectoryOutput.h"

#include <algorithm>
#include <fstream>
#include <stdexcept>
#include <string>
#include "Logger.hpp"
#include "common/rbmd_define.h"
#include "linked_cell_locator.h"
#include "memory_utils.h"
#include "spdlog/fmt/fmt.h"

TrajectoryOutput::TrajectoryOutput(double available_memory_usage_ratio)
    : _memory_usage_ratio(available_memory_usage_ratio),
      _output_thread(),
      _mutex(),
      _cv_not_full(),
      _cv_not_empty(),
      _stream(nullptr),
      _ring_buffer_size(0),
      _num_atoms(0),
      _interval(1),
      _write_index(0),
      _read_index(0),
      _stop_flag(false),
      _initialized(false) {}

// TrajectoryOutput::~TrajectoryOutput() {
//   _stop_flag = true;
//   _cv_not_empty.notify_all();
//   _cv_not_full.notify_all();
//   if (_output_thread.joinable()) {
//     _output_thread.join();
//   }
//   DeallocateRingBuffer();
//   if (_stream) {
//     CHECK_RUNTIME(STREAM_DESTORY(_stream));
//   }
// }

TrajectoryOutput::~TrajectoryOutput() {
  using Clock = std::chrono::steady_clock;
  auto t_all0 = Clock::now();

  _stop_flag = true;

  auto t0 = Clock::now();
  _cv_not_empty.notify_all();
  _cv_not_full.notify_all();
  auto t1 = Clock::now();
  double notify_s = std::chrono::duration<double>(t1 - t0).count();

  double join_s = 0.0;
  if (_output_thread.joinable()) {
    auto tj0 = Clock::now();
    _output_thread.join();
    auto tj1 = Clock::now();
    join_s = std::chrono::duration<double>(tj1 - tj0).count();
  }

  double stream_sync_s = 0.0;
  if (_stream) {
    auto ts0 = Clock::now();
    CHECK_RUNTIME(STREAM_SYNC(_stream));
    auto ts1 = Clock::now();
    stream_sync_s = std::chrono::duration<double>(ts1 - ts0).count();
  }

  auto td0 = Clock::now();
  DeallocateRingBuffer();
  auto td1 = Clock::now();
  double dealloc_s = std::chrono::duration<double>(td1 - td0).count();

  double destroy_s = 0.0;
  if (_stream) {
    auto tx0 = Clock::now();
    CHECK_RUNTIME(STREAM_DESTORY(_stream));
    auto tx1 = Clock::now();
    destroy_s = std::chrono::duration<double>(tx1 - tx0).count();
  }

  auto t_all1 = Clock::now();
  double total_s = std::chrono::duration<double>(t_all1 - t_all0).count();

  Logger::Instance().info(
      "TrajectoryOutput dtor timing: notify={:.3f} s, join={:.3f} s, stream_sync={:.3f} s, deallocate={:.3f} s, destroy_stream={:.3f} s, total={:.3f} s",
      notify_s, join_s, stream_sync_s, dealloc_s, destroy_s, total_s);
}


void TrajectoryOutput::Init() {
  _num_atoms = *(_structure_info_data->_num_atoms);
  _interval = DataManager::getInstance().getConfigData()->Get<rbmd::Id>(
      "interval", "outputs", "trajectory_out");

  if (_num_atoms == 0) {
    Logger::Instance().warn(
        "TrajectoryOutput::Init called with 0 atoms. Output will be disabled.");
    return;
  }

  size_t frame_mem = CalculateSingleFrameMemory(_num_atoms);
  size_t available_mem = rbmd::utils::get_available_memory();
  size_t memory_to_use =
      static_cast<size_t>(available_mem * _memory_usage_ratio);
  size_t calculated_size =
      (memory_to_use > frame_mem) ? (memory_to_use / frame_mem) : 1;

  _ring_buffer_size =
      std::max(static_cast<size_t>(2),
               std::min(calculated_size, static_cast<size_t>(64)));

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
  _initialized = true;
}

void TrajectoryOutput::AllocateRingBuffer(size_t num_atoms) {
  _ring_buffer.resize(_ring_buffer_size);
  for (auto& frame : _ring_buffer) {
    frame.num_atoms = num_atoms;
    frame.h_px = nullptr;
    frame.h_py = nullptr;
    frame.h_pz = nullptr;
    frame.h_vx = nullptr;
    frame.h_vy = nullptr;
    frame.h_vz = nullptr;
    frame.h_atoms_type = nullptr;
    frame.copy_complete_event = nullptr;
    CHECK_RUNTIME(MALLOCHOST(reinterpret_cast<void**>(&frame.h_px),
                             num_atoms * sizeof(rbmd::Real)));
    CHECK_RUNTIME(MALLOCHOST(reinterpret_cast<void**>(&frame.h_py),
                             num_atoms * sizeof(rbmd::Real)));
    CHECK_RUNTIME(MALLOCHOST(reinterpret_cast<void**>(&frame.h_pz),
                             num_atoms * sizeof(rbmd::Real)));
    CHECK_RUNTIME(MALLOCHOST(reinterpret_cast<void**>(&frame.h_vx),
                             num_atoms * sizeof(rbmd::Real)));
    CHECK_RUNTIME(MALLOCHOST(reinterpret_cast<void**>(&frame.h_vy),
                             num_atoms * sizeof(rbmd::Real)));
    CHECK_RUNTIME(MALLOCHOST(reinterpret_cast<void**>(&frame.h_vz),
                             num_atoms * sizeof(rbmd::Real)));
    CHECK_RUNTIME(MALLOCHOST(reinterpret_cast<void**>(&frame.h_atoms_type),
                             num_atoms * sizeof(rbmd::Id)));
    CHECK_RUNTIME(MALLOCHOST(reinterpret_cast<void**>(&frame.h_charge),
                             num_atoms * sizeof(rbmd::Real)));
    CHECK_RUNTIME(MALLOCHOST(reinterpret_cast<void**>(&frame.h_atom_id_to_idx),
                             num_atoms * sizeof(rbmd::Id)));
    CHECK_RUNTIME(NEW_FLAG_EVENT(&frame.copy_complete_event, EVENT_DISABLE));
  }
}

void TrajectoryOutput::DeallocateRingBuffer() {
  for (auto& frame : _ring_buffer) {
    CHECK_RUNTIME(FREE_PINNED_HOST(frame.h_px));
    CHECK_RUNTIME(FREE_PINNED_HOST(frame.h_py));
    CHECK_RUNTIME(FREE_PINNED_HOST(frame.h_pz));
    CHECK_RUNTIME(FREE_PINNED_HOST(frame.h_vx));
    CHECK_RUNTIME(FREE_PINNED_HOST(frame.h_vy));
    CHECK_RUNTIME(FREE_PINNED_HOST(frame.h_vz));
    CHECK_RUNTIME(FREE_PINNED_HOST(frame.h_atoms_type));
    CHECK_RUNTIME(EVENT_DESTORY(frame.copy_complete_event));
    CHECK_RUNTIME(FREE_PINNED_HOST(frame.h_charge));
    CHECK_RUNTIME(FREE_PINNED_HOST(frame.h_atom_id_to_idx));
  }
  _ring_buffer.clear();
}

size_t TrajectoryOutput::CalculateSingleFrameMemory(size_t num_atoms) const {
  size_t mem = 0;
  mem += num_atoms * sizeof(rbmd::Real) * 3;  // px, py, pz
  mem += num_atoms * sizeof(rbmd::Real) * 3;  // vx, vy, vz
  mem += num_atoms * sizeof(rbmd::Real);      // charge
  mem += num_atoms * sizeof(rbmd::Id);        // atoms_type
  mem += num_atoms * sizeof(rbmd::Id);        // atom_id_to_idx
  return mem;
}

size_t EstimateOptimalBufferSize(size_t single_frame_memory) {
  size_t ideal_size = single_frame_memory * 2;

  const size_t min_buffer = 2 * 1024 * 1024;   // min 2MB
  const size_t max_buffer = 512 * 1024 * 1024;  // max 512MB

  if (ideal_size < min_buffer) {
    ideal_size = min_buffer;
  } else if (ideal_size > max_buffer) {
    ideal_size = max_buffer;
  }

  size_t aligned_size = 1;
  while (aligned_size < ideal_size) {
    aligned_size <<= 1;
  }

  // 确保对齐后的值仍在范围内
  if (aligned_size > max_buffer) {
    aligned_size = max_buffer;
  }

  return aligned_size;
}

void TrajectoryOutput::Execute(int current_timestep) {
  if (current_timestep % _interval != 0 || _ring_buffer_size == 0) return;

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
  CHECK_RUNTIME(MEMCPY_ASYNC(target_frame.h_px, raw_ptr(_device_data->_d_px),
                             pos_bytes, D2H, _stream));
  CHECK_RUNTIME(MEMCPY_ASYNC(target_frame.h_py, raw_ptr(_device_data->_d_py),
                             pos_bytes, D2H, _stream));
  CHECK_RUNTIME(MEMCPY_ASYNC(target_frame.h_pz, raw_ptr(_device_data->_d_pz),
                             pos_bytes, D2H, _stream));
  CHECK_RUNTIME(MEMCPY_ASYNC(target_frame.h_vx, raw_ptr(_device_data->_d_vx),
                             pos_bytes, D2H, _stream));
  CHECK_RUNTIME(MEMCPY_ASYNC(target_frame.h_vy, raw_ptr(_device_data->_d_vy),
                             pos_bytes, D2H, _stream));
  CHECK_RUNTIME(MEMCPY_ASYNC(target_frame.h_vz, raw_ptr(_device_data->_d_vz),
                             pos_bytes, D2H, _stream));
  CHECK_RUNTIME(MEMCPY_ASYNC(target_frame.h_atoms_type,
                             raw_ptr(_device_data->_d_atoms_type), type_bytes,
                             D2H, _stream));
  CHECK_RUNTIME(MEMCPY_ASYNC(target_frame.h_charge,
                             raw_ptr(_device_data->_d_charge), pos_bytes, D2H,
                             _stream));
  CHECK_RUNTIME(MEMCPY_ASYNC(
      target_frame.h_atom_id_to_idx,
      raw_ptr(
          LinkedCellLocator::GetInstance().GetLinkedCell()->_atom_id_to_idx),
      type_bytes, D2H, _stream));

  // Record an event in the stream
  CHECK_RUNTIME(EVENT_RECORD(target_frame.copy_complete_event, _stream));

  lock.lock();  // Reacquire the lock to update the index
  _write_index.store((write_idx + 1) % _ring_buffer_size);
  _cv_not_empty.notify_one();
}

void TrajectoryOutput::OutputWorker() {
  std::ofstream trj_file("rbmd.trj", std::ios::out | std::ios::trunc);
  if (!trj_file.is_open()) {
    Logger::Instance().error("Failed to open trajectory file: rbmd.trj");
    return;
  }
  const size_t buffer_size = EstimateOptimalBufferSize(CalculateSingleFrameMemory(_num_atoms));
  std::unique_ptr<char[]> write_buffer(new char[buffer_size]);
  size_t buffer_pos = 0;

  // Lambda function to flush buffer
  auto flush_buffer = [&]() {
    if (buffer_pos > 0) {
      trj_file.write(write_buffer.get(), buffer_pos);
      buffer_pos = 0;
    }
  };

  //  Lambda function to write string to buffer
  auto write_to_buffer = [&](const std::string& str) {
    const char* data = str.c_str();
    size_t len = str.length();

    // If single string exceeds buffer size, write directly to file
    if (len > buffer_size) {
      flush_buffer();
      trj_file.write(data, len);
      return;
    }

    // If current data would cause buffer overflow, flush buffer first
    if (buffer_pos + len > buffer_size) {
      flush_buffer();
    }

    // Write to buffer
    std::memcpy(write_buffer.get() + buffer_pos, data, len);
    buffer_pos += len;
  };

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
    // Use fmt::format for all output formatting
    write_to_buffer(fmt::format("ITEM: TIMESTEP\n{}\n", source_frame.timestep));

    write_to_buffer(fmt::format("ITEM: NUMBER OF ATOMS\n{}\n", source_frame.num_atoms));

    write_to_buffer("ITEM: BOX BOUNDS pp pp pp\n");

    // Format box bounds using fmt
    write_to_buffer(fmt::format("{} {}\n{} {}\n{} {}\n",
                                source_frame.box_snapshot._coord_min[0],
                                source_frame.box_snapshot._coord_max[0],
                                source_frame.box_snapshot._coord_min[1],
                                source_frame.box_snapshot._coord_max[1],
                                source_frame.box_snapshot._coord_min[2],
                                source_frame.box_snapshot._coord_max[2]));

    write_to_buffer("ITEM: ATOMS id type q x y z vx vy vz\n");

    // Batch process atom data using fmt
    const size_t atoms_per_batch = 1000;
    std::string batch_buffer;
    batch_buffer.reserve(atoms_per_batch * 200);

    for (size_t i = 0; i < source_frame.num_atoms; ++i) {
      const auto idx = source_frame.h_atom_id_to_idx[i];

      // Use fmt::format for atom data - much cleaner than snprintf
      batch_buffer += fmt::format("{} {} {} {} {} {} {} {} {}\n",
                                  i + 1,
                                  source_frame.h_atoms_type[idx] + 1,
                                  source_frame.h_charge[idx],
                                  source_frame.h_px[idx],
                                  source_frame.h_py[idx],
                                  source_frame.h_pz[idx],
                                  source_frame.h_vx[idx],
                                  source_frame.h_vy[idx],
                                  source_frame.h_vz[idx]);

      // Write batch when full or at the end
      if ((i + 1) % atoms_per_batch == 0 || i == source_frame.num_atoms - 1) {
        write_to_buffer(batch_buffer);
        batch_buffer.clear();
      }
    }


    // Disabling it yields higher performance, but risks losing frame data on
    // crashes
    // flush_buffer();

    lock.lock();
    _read_index.store((read_idx + 1) % _ring_buffer_size);
    _cv_not_full.notify_one();
  }

  // Finally flush the buffer
  flush_buffer();
  trj_file.close();
}
