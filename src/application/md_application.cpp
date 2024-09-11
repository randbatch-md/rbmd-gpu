#include "md_application.h"
#include "command_line.h"
//#include "md_system.h"
#include "atomic_reader.h"
//#include "JsonParser.h"
#include "executioner.h"
#include <memory>
#include "../simulate_pipeline/include/nve_ensemble.h"
#include "../simulate_pipeline/include/nvp_ensemble.h"
#include "../simulate_pipeline/include/nvt_ensemble.h"

MDApplication::MDApplication(int argc, char* argv[]) : 
	Application(argc,argv)
{
	//_system = std::make_shared<MDSystem>();
}

int MDApplication::Execute()
{
	try
	{
		if (-1 == ReadMDData())
		{
			//log
			return -1;
		}

		if (!_config_data->HasNode("execution"))
		{
			return -1;
		}
		//_executioner = std::make_shared<Executioner>(_parser->GetJsonNode("execution"), _system);
		//_executioner = std::make_shared<Executioner>(_config_data->GetJsonNode("execution"), _simulate_pipelines);

		for (size_t i = 0; i < _simulate_pipelines.size(); i++)
		{
			_executioner = std::make_shared<Executioner>(_simulate_nodes[i], _simulate_pipelines[i]);

			_executioner->Init();

			if (-1 == _executioner->Execute())
			{
				//log
				_console->error("execute failed!");
			}
		}

	}
	catch (const std::exception&)
	{
		//log
		return -1;
	}

	return 0;
}

void MDApplication::AddSimulate()
{
	auto execution_node = _config_data->GetJsonNode("execution");
	std::vector<std::string> simulate_pipelines;
	if (execution_node.isObject()) 
	{
		simulate_pipelines = execution_node.getMemberNames();
	}

	for (const auto& simulate_pipeline : simulate_pipelines)
	{
		auto& simulate_child_node = execution_node[simulate_pipeline.c_str()];
		auto type = simulate_child_node["type"].asString(); 
		std::shared_ptr<Ensemble> ensemble;

		if ("NVT"== type)
		{
			ensemble= std::make_shared<NVTensemble>();
		}
		else if ("NVP" == type)
		{
			ensemble = std::make_shared<NVPensemble>();
		}
		else if ("NVE" == type)
		{
			ensemble = std::make_shared<NVEensemble>();
		}
		else
		{
			std::cout << " the type of execution of json file is wrong" << std::endl;
		}
		_simulate_pipelines.push_back(ensemble);
		_simulate_nodes.push_back(simulate_child_node);
	}

}

int MDApplication::ReadMDData()
{
	//auto& md_data = std::dynamic_pointer_cast<MDSystem>(_system)->GetMDData();
	std::shared_ptr<MDData> md_data = DataManager::getInstance().getMDData();
	std::shared_ptr<BaseReader> reader;
	_config_data = DataManager::getInstance().getConfigData();
	auto atom_style = _config_data->Get<std::string>("atom_style", "init_configuration", "read_data");
	if ("atomic" == atom_style)
	{
		reader = std::make_shared<AtomicReader>("rbmd.data", *md_data);
	}
	else if ("charge" == atom_style)
	{
		//reader = std::make_shared<Charge_Reader>("rbmd.data", md_data);
	}
	else if ("full" == atom_style)
	{
		//reader = std::make_shared<FullReader>("rbmd.data", md_data);
	}
	else
	{
		//log
		_console->error("ilLegal atom style!");
		return -1;
	}
	reader->Execute();

	return 0;
}
