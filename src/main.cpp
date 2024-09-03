#include "md_application.h"

int main(int argc, char* argv[]) {
  std::shared_ptr<Application> app =
      std::make_shared<MDApplication>(argc, argv);
  app->Run();

  return 0;
}

// ע�ͣ������ߡ�->ѡ��->��doxygen�� �޸�Ϊdoxygen
// ��*�� ʹ�ã�Ctrl+/ ��ϸע�ͣ���///" ����ע��