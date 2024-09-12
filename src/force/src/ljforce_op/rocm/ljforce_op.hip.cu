#include "ljforce_op.h"
#include "model/box.h"
#include <hip/hip_runtime.h>
#include "../common/rbmd_define.h"

namespace op
{

//#define THREADS_PER_BLOCK 256


	//force kernel 
    //LJForce
    __device__ void lj126(
		    rbmd::Real cut_off,
		    rbmd::Real px12,
		    rbmd::Real py12,
		    rbmd::Real pz12,
		    rbmd::Real eps_ij,
		    rbmd::Real sigma_ij,
			rbmd::Real& f_ij,
			rbmd::Real& e_ij)
	{
		const rbmd::Real  small_value = 0.0001;
		const rbmd::Real dis_2 = px12 * px12 + py12 * py12 + pz12 * pz12;
		const rbmd::Real cut_off_2 = cut_off * cut_off;

		if (dis_2 < cut_off_2 && dis_2 > small_value)
		{
			rbmd::Real sigmaij_6 = sigma_ij * sigma_ij * sigma_ij * sigma_ij * sigma_ij * sigma_ij;
			rbmd::Real dis_6 = dis_2 * dis_2 * dis_2;
			rbmd::Real sigmaij_dis_6 = sigmaij_6 / dis_6;


			 f_ij = -24 * eps_ij * ((2 * sigmaij_dis_6 - 1) * sigmaij_dis_6) / dis_2 ;
			 e_ij = 0.5 * (4 * eps_ij * (sigmaij_6 / dis_6 - 1) * (sigmaij_6 / dis_6));
		}
	}

	__global__ void ComputeLJForce(
						Box* box,
			            const rbmd::Real cut_off,
			            const rbmd::Id num_atoms,
			            const rbmd::Id* atoms_type,
			            const rbmd::Id* molecular_type,
			            const rbmd::Real* sigma,
			            const rbmd::Real* eps,
			            const rbmd::Id* start_id,
		                const rbmd::Id* end_id,
		                const rbmd::Id* id_verletlist,
			            const rbmd::Real* px,
			            const rbmd::Real* py,
			            const rbmd::Real* pz,
			            rbmd::Real* force_x,
			            rbmd::Real* force_y,
			            rbmd::Real* force_z,
			            rbmd::Real* evdwl)
	{
		rbmd::Real sum_fx = 0;
		rbmd::Real sum_fy = 0;
	    rbmd::Real sum_fz = 0;
		rbmd::Real sum_eij = 0;

		unsigned int tid1 = blockIdx.x * blockDim.x + threadIdx.x;
		if (tid1 < num_atoms)
		{
			rbmd::Id typei = atoms_type[tid1];
			rbmd::Id molecular_id_i=  molecular_type[tid1];
			rbmd::Real eps_i = eps[typei];
			rbmd::Real sigma_i = sigma[typei];

			rbmd::Real x1 = px[tid1];
			rbmd::Real y1 = py[tid1];
			rbmd::Real z1 = pz[tid1];
			printf("--------test-6-------\n");

			for (int j = start_id[tid1]; j < end_id[tid1]; ++j)
			{
				printf("--------test-7-------\n");

				rbmd::Id tid2 = id_verletlist[j];
				printf("--------test-7-1----tid2:%d---\n", tid2);

				rbmd::Id typej = atoms_type[tid2];
				printf("--------test-7-2----typej:%d---\n", typej);

				rbmd::Id molecular_id_j = molecular_type[tid2];
				printf("--------test-7-3----molecular_id_j:%d---\n", molecular_id_j);

				rbmd::Real eps_j = eps[typej];
				printf("--------test-7-4----eps_j:%f---\n", eps_j);

				rbmd::Real sigma_j = sigma[typej];
				printf("--------test-8-----sigma_j:%f---\n", sigma_j);

				//mix
				rbmd::Real eps_ij = sqrt(eps_i * eps_j);
				rbmd::Real sigma_ij = (sigma_i + sigma_j) / 2;

				rbmd::Real x2 = px[tid2];
				rbmd::Real y2 = py[tid2];
				rbmd::Real z2 = pz[tid2];
				rbmd::Real px12 = x2 - x1;
				rbmd::Real py12 = y2 - y1;
				rbmd::Real pz12 = z2 - z1;
				printf("--------test-9-------\n");

				//if (molecular_id_i == molecular_id_j)
					//continue; 

				MinImageDistance(box, px12, py12, pz12);
				rbmd::Real f_ij,e_ij;

				lj126(cut_off, px12, py12, pz12, eps_ij, sigma_ij,f_ij,e_ij);
				sum_fx += f_ij * px12;
				sum_fy += f_ij * py12;
				sum_fz += f_ij * pz12;

				sum_eij += e_ij;
				printf("--------test-10-------\n");

			}
			printf("--------test-11-------\n");

			// ʹ�� atomicAdd �����������ܣ��������ݾ���
			//atomicAdd(&force_x[tid1],sum_fx);
			//atomicAdd(&force_y[tid1],sum_fy);
			//atomicAdd(&force_z[tid1],sum_fz);
			//atomicAdd(&evdwl[tid1], sum_eij);

			force_x[tid1] += sum_fx;
			force_y[tid1] += sum_fy;
			force_z[tid1] += sum_fz;

		    evdwl[tid1] += sum_eij;
			printf("--------test-12-------\n");

		}
	}



	__device__
		rbmd::Real CoulForce(
			 rbmd::Real cut_off,
			 rbmd::Real alpha,
			 rbmd::Real charge_pi,
			 rbmd::Real charge_pj,
			 rbmd::Real px12,
			 rbmd::Real py12,
			 rbmd::Real pz12)
	{
		rbmd::Real MY_pi = 3.14159265358979323846; //pi
		rbmd::Real MY_pis = 1.77245385090551602729; // sqrt(pi)

		rbmd::Real CoulForce = 0;
		const rbmd::Real  small_value = 0.0001;
		const rbmd::Real dis_2 = px12 * px12 + py12 * py12 + pz12 * pz12;
		const rbmd::Real dis = sqrt(dis_2);
		const rbmd::Real cut_off_2 = cut_off * cut_off;

		if (dis_2 < cut_off_2 && dis_2 > small_value)
		{
			rbmd::Real erfcx = sqrt(alpha) * dis;
			rbmd::Real expx = -alpha * dis * dis;
			rbmd::Real Gnearvalue = (1.0 - erf(erfcx)) / (dis * dis) +
				2 * sqrt(alpha) * exp(expx) / (MY_pis * dis);

			CoulForce = -charge_pi * charge_pj * Gnearvalue / dis;


		}
		return CoulForce;
	}


	//LJVirial
	//__global__ void ComputeJVirial(
	//		Box* box,
	//	    const rbmd::Real cut_off,
	//		const rbmd::Id& num_atoms,
	//		const rbmd::Id* atoms_type,
	//		const rbmd::Id* molecular_type,
	//		const rbmd::Real* eps,
	//		const rbmd::Real* sigma,
	//		const rbmd::Id* start_id,
	//		const rbmd::Id* end_id,
	//		const rbmd::Id* id_verletlist,
	//		const rbmd::Real* px,
	//		const rbmd::Real* py,
	//		const rbmd::Real* pz,
	//		rbmd::Real* virial_xx,
	//		rbmd::Real* virial_yy,
	//		rbmd::Real* virial_zz,
	//		rbmd::Real* virial_xy,
	//		rbmd::Real* virial_xz,
	//		rbmd::Real* virial_yz)
	//{

	//	rbmd::Real sum_virial_xx = 0;
	//	rbmd::Real sum_virial_yy= 0;
	//	rbmd::Real sum_virial_zz = 0;
	//	rbmd::Real sum_virial_xy = 0;
	//	rbmd::Real sum_virial_xz = 0;
	//	rbmd::Real sum_virial_yz = 0;

	//	int tid1 = threadIdx.x + blockIdx.x * blockDim.x;
	//	if (tid1 < num_atoms)
	//	{
	//		rbmd::Id typei = atoms_type[tid1];
	//		rbmd::Id molecular_id_i = molecular_type[tid1];
	//		rbmd::Real eps_i = eps[typei];
	//		rbmd::Real sigma_i = sigma[typei];

	//		rbmd::Real x1 = px[tid1];
	//		rbmd::Real y1 = py[tid1];
	//		rbmd::Real z1 = pz[tid1];


	//		for (int j = start_id[tid1]; j < end_id[tid1]; ++j)
	//		{
	//			rbmd::Id tid2 = id_verletlist[j];

	//			rbmd::Id typej = atoms_type[tid2];
	//			rbmd::Id molecular_id_j = molecular_type[tid2];
	//			rbmd::Real eps_j = eps[typej];
	//			rbmd::Real sigma_j = sigma[typej];

	//			//mix
	//			rbmd::Real eps_ij = sqrt(eps_i * eps_j);
	//			rbmd::Real sigma_ij = (sigma_i + sigma_j) / 2;

	//			rbmd::Real x2 = px[tid2];
	//			rbmd::Real y2 = py[tid2];
	//			rbmd::Real z2 = pz[tid2];
	//			rbmd::Real px12 = x2 - x1;
	//			rbmd::Real py12 = y2 - y1;
	//			rbmd::Real pz12 = z2 - z1;

	//			if (molecular_id_i == molecular_id_j)
	//				continue;

	//			//MinMirror(box, px12, py12, pz12);
	//			rbmd::Real Virial_f;
	//			Virial_f = LJVirial(cut_off, px12, py12, pz12, eps_ij, sigma_ij);

	//			rbmd::Real Virial_fx = Virial_f * px12;
	//			rbmd::Real Virial_fy = Virial_f * py12;
	//			rbmd::Real Virial_fz = Virial_f * pz12;

	//			sum_virial_xx +=  px12 * Virial_fx;
	//			sum_virial_yy +=  py12 * Virial_fy;
	//			sum_virial_zz +=  pz12 * Virial_fz;
	//			sum_virial_xy +=  px12 * Virial_fy;
	//			sum_virial_xz +=  px12 * Virial_fz;
	//			sum_virial_yz +=  py12 * Virial_fz;
	//		}

	//		// save virial
	//		virial_xx[tid1] += sum_virial_xx;
	//		virial_yy[tid1] += sum_virial_yy;
	//		virial_zz[tid1] += sum_virial_zz;
	//		virial_xy[tid1] += sum_virial_xy;
	//		virial_xz[tid1] += sum_virial_xz;
	//		virial_yz[tid1] += sum_virial_yz;

	//	}
	//}

	//__device__ rbmd::Real LJVirial(
	//		const rbmd::Real cut_off,
	//		const rbmd::Real px12,
	//		const rbmd::Real py12,
	//		const rbmd::Real pz12,
	//		const rbmd::Real eps_ij,
	//		const rbmd::Real sigma_ij)
	//{
	//	rbmd::Real virial_f = 0;;

	//	const rbmd::Real  small_value = 0.0001;
	//	const rbmd::Real dis_2 = px12 * px12 + py12 * py12 + pz12 * pz12;
	//	const rbmd::Real cut_off_2 = cut_off * cut_off;

	//	if (dis_2 < cut_off_2 && dis_2 > small_value)
	//	{
	//		rbmd::Real sigmaij_6 = sigma_ij * sigma_ij * sigma_ij * sigma_ij * sigma_ij * sigma_ij;
	//		rbmd::Real dis_6 = dis_2 * dis_2 * dis_2;
	//		rbmd::Real sigmaij_dis_6 = sigmaij_6 / dis_6;

	//		 virial_f = 0.5 * 24 * eps_ij * ((2 * sigmaij_dis_6 - 1) * sigmaij_dis_6) / dis_2;

	//}
	//	return virial_f;
	//}


	//
	void LJForceOp<device::DEVICE_GPU>::operator()(
						Box* box,
			            const rbmd::Real cut_off,
			            const rbmd::Id num_atoms,
		                const rbmd::Id* atoms_type,
			            const rbmd::Id* molecular_type,
			            const rbmd::Real* sigma,
			            const rbmd::Real* eps,
			            const rbmd::Id* start_id,
		                const rbmd::Id* end_id,
		                const rbmd::Id* id_verletlist,
			            const rbmd::Real* px,
			            const rbmd::Real* py,
			            const rbmd::Real* pz,
			            rbmd::Real* force_x,
			            rbmd::Real* force_y,
			            rbmd::Real* force_z,
			            rbmd::Real* evdwl)
		{
		    unsigned int blocks_per_grid = (num_atoms + BLOCK_SIZE - 1) / BLOCK_SIZE;
			printf("value of a ：%f\n", cut_off);
			printf("value of num_atoms ：%d\n", num_atoms);
			printf("value of atoms_type ：%d\n", atoms_type[0]);
			printf("value of molecular_type ：%d\n", molecular_type[0]);
			printf("value of sigma ：%f\n", sigma[0]);
			printf("value of eps ：%f\n", eps[0]);
			printf("value of start_id ：%d\n", start_id[0]);
			printf("value of end_id ：%d\n", end_id[0]);
			printf("value of id_verletlist ：%d\n", id_verletlist[0]);
			printf("value of px ：%f\n", px[0]);
			printf("value of py ：%f\n", py[0]);
			printf("value of pz ：%f\n", pz[0]);
			printf("value of force_x ：%f\n", force_x[0]);
			printf("value of force_y ：%f\n", force_y[0]);
			printf("value of force_z ：%f\n", force_z[0]);
			printf("value of evdwl ：%f\n", evdwl[0]);



		    CHECK_KERNEL(ComputeLJForce <<<blocks_per_grid, BLOCK_SIZE, 0, 0 >>> (box, cut_off, num_atoms, atoms_type, molecular_type,
				sigma, eps, start_id, end_id, id_verletlist, px, py, pz, force_x, force_y, force_z, evdwl));
		}

}

