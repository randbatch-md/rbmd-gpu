#include "md_application.h"
#include "command_line.h"
#include "md_system.h"
#include "atomic_reader.h"
//#include "full_reader .h"
//#include "charge_reader.h"
#include "JsonParser.h"
#include "executioner.h"
#include <memory>

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
		_executioner = std::make_shared<Executioner>(_config_data->GetJsonNode("execution"), _simulate_pipeline);

		
		_executioner->Init();

		if (-1 == _executioner->Execute())
		{
			//log
			_console->error("execute failed!");
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

}

int MDApplication::ReadMDData()
{
	auto& md_data = std::dynamic_pointer_cast<MDSystem>(_system)->GetMDData();
	std::shared_ptr<BaseReader> reader;
	auto atom_style = _config_data->Get<std::string>("atom_style", "init_configuration", "read_data");
	if ("atomic" == atom_style)
	{
		reader = std::make_shared<AtomicReader>("rbmd.data", md_data);
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
