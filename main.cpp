#include "./src/reader/atomic_reader.h"
#include "md_data.h"

int main()
{
	MDData data;

	BaseReader* reader = new AtomicReader("./rbmd_atomic.data", data);
	reader->Execute();

	//potential


	return 0;
}