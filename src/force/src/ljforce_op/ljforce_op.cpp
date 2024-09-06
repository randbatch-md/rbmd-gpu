#include "ljforce_op.h"


namespace op
{
	//template <typename FPTYPE>
	//struct ComputeLJforceOp<FPTYPE, device::DEVICE_CPU>
	//{
	//	//to do openmp
	//	void operator()(Box& box,
	//		            const rbmd::Id& N,
	//		            const rbmd::Id* atoms_type,
	//		            const rbmd::Id* molecular_type,
	//		            const rbmd::Real* sigma,
	//		            const rbmd::Real* eps,
	//		            const rbmd::Id* start_id,
	//	                const rbmd::Id* end_id,
	//	                const rbmd::Id* id_verletlist,
	//		            const rbmd::Real* px,
	//		            const rbmd::Real* py,
	//		            const rbmd::Real* pz,
	//		            rbmd::Real* force_x,
	//		            rbmd::Real* force_y,
	//		            rbmd::Real* force_z,
	//		            rbmd::Real* evdwl);
	//	{
	//	}
	//};

	template <typename T>
    void DeviceReduceSum(const T* d_data, T* d_result, int num_elements)
    {
        int temp_storage_bytes = 0;
        void* d_temp_storage = nullptr;

        // ��һ�ε��ã�ȷ����Ҫ����ʱ�洢��С
        hipcub::DeviceReduce::Reduce(d_temp_storage, temp_storage_bytes,
                                     d_data, d_result, num_elements, hipcub::Sum(), 0.0);
        
        // ������ʱ�洢
        hipMalloc(&d_temp_storage, temp_storage_bytes);

        // �ڶ��ε��ã�ִ�й�Լ
        hipcub::DeviceReduce::Reduce(d_temp_storage, temp_storage_bytes,
                                     d_data, d_result, num_elements, hipcub::Sum(), 0.0);

        // �ͷ���ʱ�洢
        hipFree(d_temp_storage);
    }

    // ��ȷʵ����
    template void DeviceReduceSum<float>(const float*, float*, int);
    template void DeviceReduceSum<double>(const double*, double*, int);
}