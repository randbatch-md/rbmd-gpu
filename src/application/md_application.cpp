#include "md_application.h"

#include <spdlog/sinks/basic_file_sink.h>
#include <spdlog/sinks/stdout_color_sinks.h>
#include <memory>

#include "../common/json.hpp"
#include "../simulate_pipeline/include/npt_ensemble.h"
#include "../simulate_pipeline/include/nve_ensemble.h"
#include "../simulate_pipeline/include/nvt_ensemble.h"
#include "atomic_reader.h"
#include "command_line.h"
#include "cvff_memory_scheduler.h"
#include "lj_memory_scheduler.h"
#include "memory_scheduler.h"
#include "output/include/TrajectoryOutput.h"
#include "output/include/linux_pipe_sink.hpp"
#include "output/include/Logger.hpp"
MDApplication::MDApplication(int argc, char* argv[])
  : Application(argc, argv),
  _cmd(std::make_shared<CommandLine>(argc, argv))
{
  DataManager::Initialize(_cmd->GetConfigPath());

  const std::string divider = "================================================";
  const std::string LOGO = R"(
         _____    ____    __  __   _____
        |  __ \  |  _ \  |  \/  | |  __ \
        | |__) | | |_) | | \  / | | |  | |
        |  _  /  |  _ <  | |\/| | | |  | |
        | | \ \  | |_) | | |  | | | |__| |
        |_|  \_\ |____/  |_|  |_| |_____/

  )";
  std::string Logo_output ="\n" + divider + "\n" + LOGO + "\n" + divider + "\n";
  auto console_sink = std::make_shared<spdlog::sinks::stdout_color_sink_mt>();
  auto file_sink = std::make_shared<spdlog::sinks::basic_file_sink_mt>("./rbmd.log",true);
  std::vector<spdlog::sink_ptr> sinks;
  // Create a named pipe sink
#ifdef WITH_GUI
  auto pipe_sink = std::make_shared<linux_pipe_sink_mt>("/tmp/rbmd_log_pipe");
  sinks.push_back(pipe_sink);
#endif
  sinks.push_back(console_sink);
  sinks.push_back(file_sink);
  auto combined_logger = std::make_shared<spdlog::logger>("logger", sinks.begin(), sinks.end());
  combined_logger->set_pattern("[%Y-%m-%d %H:%M:%S.%02e] [%^%l%$] %v");
  combined_logger->set_level(spdlog::level::info);
  Logger::Instance().Configure(0,combined_logger);
  Logger::Instance().info(Logo_output.c_str());
  //
  ReadMDData();
}

int MDApplication::Execute() {


  //ReadMDData();
  std::string ensemble_type = DataManager::getInstance().getConfigData()->Get
    <std::string>("ensemble", "execution");
  if("NVE" == ensemble_type)
  {
    _simulate_pipeline = std::make_shared<NVEensemble>();
  }
  else if("NVT" == ensemble_type)
  {
    _simulate_pipeline = std::make_shared<NVTensemble>();
  }
  else if("NPT" == ensemble_type)
  {
    _simulate_pipeline = std::make_shared<NPTensemble>();
  }
  _output = std::make_shared<TrajectoryOutput>();
  _simulate = std::make_shared<Simulate>(_simulate_pipeline,_output);

  _simulate->Init();

  _simulate->Execute();

  const std::string divider = "================================================";
  std::string finish_output =divider + "        Finish        " + divider ;
  Logger::Instance().info(finish_output.c_str());

  DataManager::getInstance().unloadDeviceData();
  return 0;
}

void MDApplication::AddSimulate() {
  auto execution_node = _config_data->GetJsonNode("execution");
  std::vector<std::string> simulate_pipelines;
  if (!execution_node.is_object()) {
    return;
  }

  //for (const auto& [simulate_pipeline, simulate_child_node] : execution_node.items()) {
  for (const auto& item : execution_node.items())
  {
    const std::string& simulate_pipeline = item.key();
    const auto& simulate_child_node = item.value();
    auto type = simulate_child_node["type"].get<std::string>();
    std::shared_ptr<Ensemble> ensemble;

    if ("NVT" == type) {
      ensemble = std::make_shared<NVTensemble>();
    } else if ("NPT" == type) {
      ensemble = std::make_shared<NPTensemble>();
    } else if ("NVE" == type) {
      ensemble = std::make_shared<NVEensemble>();
    } else {
      Logger::Instance().error(" the type of execution of json file is wrong");
    }
    _simulate_pipelines.push_back(std::move(ensemble));
    _simulate_nodes.emplace_back(simulate_child_node);
  }
}

int MDApplication::ReadMDData() {
  ////auto& md_data = std::dynamic_pointer_cast<MDSystem>(_system)->GetMDData();
  // auto md_data = DataManager::getInstance().getMDData().get();
  // std::shared_ptr<BaseReader> reader;
  //_config_data = DataManager::getInstance().getConfigData();
  // auto atom_style = _config_data->Get<std::string>("atom_style",
  // "init_configuration", "read_data"); if ("atomic" == atom_style)
  //{
  //	reader = std::make_shared<AtomicReader>("rbmd.data", *md_data);
  // }
  // else if ("charge" == atom_style)
  //{
  //	//reader = std::make_shared<Charge_Reader>("rbmd.data", md_data);
  // }
  // else if ("full" == atom_style)
  //{
  //	//reader = std::make_shared<FullReader>("rbmd.data", md_data);
  // }
  // else
  //{
  //	//log
  //	//_console->error("ilLegal atom style!");
  //	return -1;
  // }
  std::shared_ptr<BaseReader> reader;
  std::shared_ptr<MDData> md_data = DataManager::getInstance().getMDData();

  auto config_data = DataManager::getInstance().getConfigData();
  std::string file_path = config_data->Get<std::string>("file", "init_configuration", "read_data");
  reader = std::make_shared<AtomicReader>(file_path, *md_data);
  //reader = std::make_shared<AtomicReader>("rbmd.data", *md_data);
  reader->Execute();

  std::shared_ptr<MemoryScheduler> memory_scheduler;
  auto force_type = DataManager::getInstance().getConfigData()->Get<std::string>("type", "hyper_parameters", "force_field");
  if ("CVFF" == force_type) {
      memory_scheduler = std::make_shared<CVFFMemoryScheduler>();
  }
  else
  {
      memory_scheduler = std::make_shared<LJMemoryScheduler>();
  }

  DataManager::getInstance().Fill2Device(memory_scheduler);
  return 0;
}
