#include "./src/reader/atomic_reader.h"

int main()
{
	MmapReader* reader = new AtomicReader("./lj_6w.data");
	reader->Execute();

	return 0;
}