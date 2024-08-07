#include "md_application.h"

int main(int argc, char*argv[])
{
	std::shared_ptr<Application> app = std::make_shared<MDApplication>(argc, argv);
	app->Run();

	return 0;
}
// 注释：“工具”->选项->“doxygen” 修改为doxygen 、*；
// 使用：Ctrl+/ 详细注释；“///" 单个注释