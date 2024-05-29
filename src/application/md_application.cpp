#include "md_application.h"
#include "command_line.h"
#include "md_system.h"
#include "atomic_reader.h"

MDApplication::MDApplication(int argc, char* argv[]) : 
	Application(argc,argv)
{

}

int MDApplication::Execute()
{
	try
	{
		_command_line->GetFile();
		//_parser = std::
		//_reader = std::
		_system = std::make_shared<MDSystem>();
		auto& md_data = std::dynamic_pointer_cast<MDSystem>(_system)->GetMDData();

		_reader = std::make_shared<AtomicReader>("./rbmd_atomic.data", md_data);
		_reader->Execute();
	}
	catch (const std::exception&)
	{
		//log
		return -1;
	}

	return 0;
}