#include "./include/lj_memory_scheduler.h"

//const int lj_num_streams = 7

bool LJMemoryScheduler::asyncMemoryH2D()
{
	//hipStream_t streams[lj_num_streams];




	//for (int i = 0; i < lj_num_streams; ++i)
	//{
	//	hipStreamSynchronize(streams[i]);
	//}
	return true;
}

bool LJMemoryScheduler::asyncMemoryD2H()
{
	return true;
}
