#include "base_reader.h"

class ConfReader : public BaseReader
{
public:
	ConfReader(const std::string& filePath); 
	~ConfReader() = default;

	void Execute() override;
};