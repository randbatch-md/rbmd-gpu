#include <hip/hip_runtime.h>
#include "lj_memory_scheduler.h"

const int lj_num_streams = 7

int LJMemoryScheduler::AsyncMemoryH2D()
{
	hipStream_t streams[lj_num_streams];




	for (int i = 0; i < lj_num_streams; ++i)
	{
		hipStreamSynchronize(streams[i]);
	}

}

int LJMemoryScheduler::AsyncMemoryD2H()
{

}
