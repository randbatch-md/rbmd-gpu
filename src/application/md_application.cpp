#include "md_application.h"
#include "command_line.h"
#include "md_system.h"
#include "atomic_reader.h"
#include "JsonParser.h"
#include "executioner.h"
#include <memory>

MDApplication::MDApplication(int argc, char* argv[]) : 
	Application(argc,argv)
{
	_system = std::make_shared<MDSystem>();
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

		if (!_parser->HasNode("execution"))
		{
			return -1;
		}
		_executioner = std::make_shared<Executioner>(_parser->GetJsonNode("execution"), _system);

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

int MDApplication::ReadMDData()
{
	auto& md_data = std::dynamic_pointer_cast<MDSystem>(_system)->GetMDData();
	std::shared_ptr<BaseReader> reader;
	auto atom_style = _parser->Get<std::string>("atom_style", "init_configuration", "read_data");
	if ("atomic" == atom_style)
	{
		reader = std::make_shared<AtomicReader>("rbmd.data", md_data);
	}
	else if ("charge" == atom_style)
	{
		//todo
	}
	else if ("full" == atom_style)
	{
		//todo
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
