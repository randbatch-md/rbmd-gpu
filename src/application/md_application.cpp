#include "md_application.h"
#include "command_line.h"
#include "md_system.h"
#include "atomic_reader.h"

MDApplication::MDApplication(int argc, char* argv[])
{
	_command_line = std::make_shared<CommandLine>(argc, argv);

	//if (_command_line->IsRunApplication())
	//{
		//_parser = std::
		//_reader = std::
		_system = std::make_shared<MDSystem>();
		auto& md_data = std::dynamic_pointer_cast<MDSystem>(_system)->GetMDData();

		_reader = std::make_shared<AtomicReader>("./rbmd_atomic.data", md_data);
		_reader->Execute();
	//}


}

void MDApplication::Run()
{
	if (!_command_line->IsRunApplication())
	{
		//log
		return;
	}


}