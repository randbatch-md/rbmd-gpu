#pragma once
#include <hip/hip_runtime.h>

#define hipErrorCheck(res) {\
	if(res != hipSuccess){\
		fprintf(stderr, "unexpected Device Error %s:%d: %s, %s\n", __FILE__, __LINE__, hipGetErrorName(res), hipGetErrorString(res));\
		exit(res);\
	}\
}