#include "./src/reader/atomic_reader.h"
#include "md_data.h"

int main()
{
	MDData data;

	MmapReader* reader = new AtomicReader("./lj_6w.data", data);
	reader->Execute();

	//potential


	return 0;
}