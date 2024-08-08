#pragma once
#include <memory>

struct MDData;
struct PostProcessData;
class ConfigData;
struct DeviceData;
class MemoryScheduler;

class DataManager
{
public:
	/**
	 * @brief The copy constructor and assignment operator are not open to the public
	 * @param  
	*/
	DataManager(const DataManager&) = delete;
	DataManager& operator=(const DataManager&) = delete;

private:
	/**
	 * @brief Constructors and destructors are not open to the public
	*/
	DataManager() = default;
	~DataManager() = default;

public:
	static DataManager& getInstance()
	{
		static DataManager data_manager;
		return data_manager;
	}

	auto& getMDData() const { return _md_data; }
	auto& getPostProcessData() const { return _postProcess_data; }
	auto& getConfigData() const { return _config_data; }
	auto& getDeviceData() const { return _device_data; }
	auto& getMemoryScheduler() const { return _memory_scheduler; }

private:
	std::shared_ptr<MDData> _md_data;
	std::shared_ptr<PostProcessData> _postProcess_data;
	std::shared_ptr<ConfigData> _config_data;
	std::shared_ptr<DeviceData> _device_data;
	std::shared_ptr<MemoryScheduler> _memory_scheduler;
};
