#pragma once

namespace device {
    struct DEVICE_CPU;
    struct DEVICE_GPU;
} // namespace device

// enum class DeviceType
//{
//	UnKnown = 0,
//	CpuDevice = 1,
//	GpuDevice = 2,
// };
//
// template <typename T>
// struct DeviceTypeToEnum{
//	static constexpr DeviceType value = {};
// };
//
// template<>
// struct DeviceTypeToEnum<DEVICE_CPU> {
//	static constexpr DeviceType value = DeviceType::CpuDevice;
// };
//
// template<>
// struct DeviceTypeToEnum<DEVICE_GPU> {
//	static constexpr DeviceType value = DeviceType::GpuDevice;
// };
