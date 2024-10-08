#include "ljforce_op.h"
#include "model/box.h"
#include <hip/hip_runtime.h>
#include "../common/rbmd_define.h"

namespace op
{

    //---------device---------//
    //lj126
    __device__ void lj126(
		    rbmd::Real cut_off,
		    rbmd::Real px12,
		    rbmd::Real py12,
		    rbmd::Real pz12,
		    rbmd::Real eps_ij,
		    rbmd::Real sigma_ij,
                    rbmd::Real& force_lj,
                    rbmd::Real& energy_lj)
	{
		const rbmd::Real  small_value = 0.0001;
		const rbmd::Real dis_2 = px12 * px12 + py12 * py12 + pz12 * pz12;
		const rbmd::Real cut_off_2 = cut_off * cut_off;

		if (dis_2 < cut_off_2 && dis_2 > small_value)
		{
		   rbmd::Real sigmaij_6 = POW(sigma_ij,6.0);
		   rbmd::Real dis_6 =  POW(dis_2,3.0);
		   rbmd::Real sigmaij_dis_6 = sigmaij_6 / dis_6;

                   force_lj = -24 * eps_ij * ((2 * sigmaij_dis_6 - 1) * sigmaij_dis_6) / dis_2 ;
                   energy_lj = 0.5 * (4 * eps_ij * (sigmaij_6 / dis_6 - 1) * (sigmaij_6 / dis_6));
		}
	}

	//lj126_rs
	__device__ void lj126_rs(
		rbmd::Real rs,
		rbmd::Real px12,
		rbmd::Real py12,
		rbmd::Real pz12,
		rbmd::Real eps_ij,
		rbmd::Real sigma_ij,
		rbmd::Real& fs_ij)
	{
		const rbmd::Real  small_value = 0.0001;
		const rbmd::Real dis_2 = px12 * px12 + py12 * py12 + pz12 * pz12;
		const rbmd::Real rs_2 = rs * rs;

		if (dis_2 < rs_2 && dis_2 > small_value)
		{
			rbmd::Real sigmaij_6 = POW(sigma_ij,6.0);
			rbmd::Real dis_6 = POW(dis_2,3.0);
			rbmd::Real sigmaij_dis_6 = sigmaij_6 / dis_6;

			fs_ij = -24 * eps_ij * ((2 * sigmaij_dis_6 - 1) * sigmaij_dis_6) / dis_2;
		}
	}

	//lj126_rcs
	__device__ void lj126_rcs(
		rbmd::Real rc,
		rbmd::Real rs,
		rbmd::Id  pice_num,
		rbmd::Real px12,
		rbmd::Real py12,
		rbmd::Real pz12,
		rbmd::Real eps_ij,
		rbmd::Real sigma_ij,
		rbmd::Real& fcs_ij)
	{
		const rbmd::Real  small_value = 0.0001;
		const rbmd::Real dis_2 = px12 * px12 + py12 * py12 + pz12 * pz12;
		const rbmd::Real rc_2 = rc * rc;
		const rbmd::Real rs_2 = rs * rs;

		if (dis_2 < rc_2 && dis_2 > rs_2)
		{
			rbmd::Real sigmaij_6 = POW(sigma_ij,6.0);
			rbmd::Real dis_6 = POW(dis_2,3.0);
			rbmd::Real sigmaij_dis_6 = sigmaij_6 / dis_6;

			fcs_ij = pice_num * (-24 * eps_ij * ((2 * sigmaij_dis_6 - 1) * sigmaij_dis_6) / dis_2);
		}
	}

	//CoulForce
	__device__
		void CoulCutForce(
			rbmd::Real cut_off,
			rbmd::Real alpha,
			rbmd::Real qqr2e,
			rbmd::Real charge_i,
			rbmd::Real charge_j,
			rbmd::Real px12,
			rbmd::Real py12,
			rbmd::Real pz12,
			rbmd::Real force_coul,
                        rbmd::Real energy_coul)
	{
		rbmd::Real MY_pi = 3.14159265358979323846; //pi
		rbmd::Real MY_pis = 1.77245385090551602729; // sqrt(pi)

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

                   force_coul = qqr2e * (-charge_i * charge_j * Gnearvalue / dis);
                   energy_coul = qqr2e * (0.5 * charge_i * charge_j * (1.0 - erf(sqrt(alpha) * dis)) / dis);
		}
	}

	//EwaldForce
	__device__ void EwaldForce(
		Box* box,
		const rbmd::Real alpha,
		Real3 M,
		const rbmd::Real rhok_real_i,
		const rbmd::Real rhok_imag_i,
		const rbmd::Real charge,
		const rbmd::Real px,
		const rbmd::Real py,
		const rbmd::Real pz,
		rbmd::Real& force_ewald_x,
		rbmd::Real& force_ewald_y,
		rbmd::Real& force_ewald_z)
	{
		rbmd::Real force_ewald;
		rbmd::Real volume = box->_length[0] * box->_length[1] * box->_length[2];
		Real3 K = make_Real3(2 * M_PI * M.x / box->_length[0],
		                       2 * M_PI * M.y / box->_length[1],
			                 2 * M_PI * M.z / box->_length[2]);


		rbmd::Real range_K_2 = K.x * K.x + K.y * K.y + K.z * K.z;
		rbmd::Real dot_product = K.x * px + K.y * py + K.z * pz;

		rbmd::Real factor_a = -4 * M_PI * charge;
		rbmd::Real factor_b = exp(-range_K_2 / (4 * alpha));
		rbmd::Real factor_c = cos(dot_product) * rhok_imag_i;
		rbmd::Real factor_d = sin(dot_product) * rhok_real_i;

		force_ewald = factor_a / (volume * range_K_2) * factor_b * (factor_c - factor_d);
		force_ewald_x = force_ewald * K.x;
		force_ewald_y = force_ewald * K.y;
		force_ewald_z = force_ewald * K.z;
	}


	//------global---------//
	//verlet-list: LJForce
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
            rbmd::Real* fx,
            rbmd::Real* fy,
            rbmd::Real* fz,
            rbmd::Real* total_evdwl)
        {
          rbmd::Real sum_fx = 0;
          rbmd::Real sum_fy = 0;
          rbmd::Real sum_fz = 0;

          rbmd::Real sum_elj = 0;

          unsigned int tid1 = blockIdx.x * blockDim.x + threadIdx.x;
          if (tid1 < num_atoms)
          {
            rbmd::Id typei = atoms_type[tid1];
            rbmd::Id molecular_id_i =  molecular_type[tid1];
            rbmd::Real eps_i = eps[typei];
            rbmd::Real sigma_i = sigma[typei];

            rbmd::Real x1 = px[tid1];
            rbmd::Real y1 = py[tid1];
            rbmd::Real z1 = pz[tid1];


            for (int j = start_id[tid1]; j < end_id[tid1]; ++j)
            {

              rbmd::Id tid2 = id_verletlist[j];
              rbmd::Id typej = atoms_type[tid2];
              rbmd::Id molecular_id_j = molecular_type[tid2];
              rbmd::Real eps_j = eps[typej];
              rbmd::Real sigma_j = sigma[typej];

              //mix
              rbmd::Real eps_ij = sqrt(eps_i * eps_j);
              rbmd::Real sigma_ij = (sigma_i + sigma_j) / 2;

              rbmd::Real x2 = px[tid2];
              rbmd::Real y2 = py[tid2];
              rbmd::Real z2 = pz[tid2];
              rbmd::Real px12 = x2 - x1;
              rbmd::Real py12 = y2 - y1;
              rbmd::Real pz12 = z2 - z1;
              //if (molecular_id_i == molecular_id_j)
              //continue;

              MinImageDistance(box, px12, py12, pz12);

              rbmd::Real force_lj,fpair;
              rbmd::Real energy_lj;

              lj126(cut_off, px12, py12, pz12, eps_ij, sigma_ij,force_lj,energy_lj);
              sum_fx += force_lj * px12;
              sum_fy += force_lj * py12;
              sum_fz += force_lj * pz12;

              sum_elj += energy_lj;
            }

            fx[tid1] = sum_fx;
            fy[tid1] = sum_fy;
            fz[tid1] = sum_fz;

            atomicAdd(total_evdwl, sum_elj);
            //printf("--------test---evdwl[tid1]:%f---\n",evdwl[tid1]);
          }
        }


        __global__ void ComputeLJForceViral(
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
			            rbmd::Real* fx,
			            rbmd::Real* fy,
			            rbmd::Real* fz,
				    rbmd::Real* virial_xx,
		                    rbmd::Real* virial_yy,
		                    rbmd::Real* virial_zz,
		                    rbmd::Real* virial_xy,
		                    rbmd::Real* virial_xz,
                                    rbmd::Real* virial_yz,
				    rbmd::Real* total_evdwl)
	{
		rbmd::Real sum_fx = 0;
		rbmd::Real sum_fy = 0;
	        rbmd::Real sum_fz = 0;

		rbmd::Real sum_virial_xx = 0;
		rbmd::Real sum_virial_yy = 0;
		rbmd::Real sum_virial_zz = 0;
		rbmd::Real sum_virial_xy = 0;
		rbmd::Real sum_virial_xz = 0;
		rbmd::Real sum_virial_yz = 0;

		rbmd::Real sum_elj = 0;

		unsigned int tid1 = blockIdx.x * blockDim.x + threadIdx.x;
		if (tid1 < num_atoms)
		{
			rbmd::Id typei = atoms_type[tid1]; 
			rbmd::Id molecular_id_i =  molecular_type[tid1];
			rbmd::Real eps_i = eps[typei];
			rbmd::Real sigma_i = sigma[typei];

			rbmd::Real x1 = px[tid1];
			rbmd::Real y1 = py[tid1];
			rbmd::Real z1 = pz[tid1];


			for (int j = start_id[tid1]; j < end_id[tid1]; ++j)
			{

				rbmd::Id tid2 = id_verletlist[j];
				rbmd::Id typej = atoms_type[tid2];
				rbmd::Id molecular_id_j = molecular_type[tid2];
				rbmd::Real eps_j = eps[typej];
				rbmd::Real sigma_j = sigma[typej];

				//mix
				rbmd::Real eps_ij = sqrt(eps_i * eps_j);
				rbmd::Real sigma_ij = (sigma_i + sigma_j) / 2;

				rbmd::Real x2 = px[tid2];
				rbmd::Real y2 = py[tid2];
				rbmd::Real z2 = pz[tid2];
				rbmd::Real px12 = x2 - x1;
				rbmd::Real py12 = y2 - y1;
				rbmd::Real pz12 = z2 - z1;
				//if (molecular_id_i == molecular_id_j)
					//continue; 

				MinImageDistance(box, px12, py12, pz12);

				rbmd::Real force_lj,fpair;
				rbmd::Real energy_lj;

				lj126(cut_off, px12, py12, pz12, eps_ij, sigma_ij,force_lj,energy_lj);
				sum_fx += force_lj * px12;
				sum_fy += force_lj * py12;
				sum_fz += force_lj * pz12;

				fpair = -0.5 * force_lj;
				sum_virial_xx +=  px12 * px12 * fpair;
				sum_virial_yy +=  py12 * py12 * fpair;
				sum_virial_zz +=  pz12 * pz12 * fpair;
				sum_virial_xy +=  px12 * py12 * fpair;
				sum_virial_xz +=  px12 * pz12 * fpair;
				sum_virial_yz +=  py12 * pz12 * fpair;

				sum_elj += energy_lj;
			}

			fx[tid1] = sum_fx;
			fy[tid1] = sum_fy;
			fz[tid1] = sum_fz;

			virial_xx[tid1] = sum_virial_xx;
			virial_yy[tid1] = sum_virial_yy;
			virial_zz[tid1] = sum_virial_zz;
			virial_xy[tid1] = sum_virial_xy;
			virial_xz[tid1] = sum_virial_xz;
			virial_yz[tid1] = sum_virial_yz;

                        atomicAdd(total_evdwl, sum_elj);
			//printf("--------test---evdwl[tid1]:%f---\n",evdwl[tid1]);
		}
		
	}

	//RBL: LJForce
	__global__ void ComputeLJRBLForce(
		Box* box,
		const rbmd::Real rs,
		const rbmd::Real rc,
		const rbmd::Id num_atoms,
		const rbmd::Id neighbor_sample_num,
		const rbmd::Id pice_num,
		const rbmd::Id* atoms_type,
		const rbmd::Id* molecular_type,
		const rbmd::Real* sigma,
		const rbmd::Real* eps,
		const rbmd::Id* start_id,
		const rbmd::Id* end_id,
		const rbmd::Id* id_verletlist,
		const rbmd::Id* id_random_neighbor,
		const rbmd::Id* random_neighbor_num,
		const rbmd::Real* px,
		const rbmd::Real* py,
		const rbmd::Real* pz,
		rbmd::Real* fx,
                rbmd::Real* fy,
                rbmd::Real* fz)
	{
		rbmd::Real sum_fx = 0;
		rbmd::Real sum_fy = 0;
		rbmd::Real sum_fz = 0;

		rbmd::Real sum_fsx = 0;
		rbmd::Real sum_fsy = 0;
		rbmd::Real sum_fsz = 0;

		rbmd::Real sum_fcsx = 0;
		rbmd::Real sum_fcsy = 0;
		rbmd::Real sum_fcsz = 0;



		unsigned int tid1 = blockIdx.x * blockDim.x + threadIdx.x;
		if (tid1 < num_atoms)
		{
			rbmd::Id typei = atoms_type[tid1];
			rbmd::Id molecular_id_i = molecular_type[tid1];
			rbmd::Real eps_i = eps[typei];
			rbmd::Real sigma_i = sigma[typei];

			rbmd::Real x1 = px[tid1];
			rbmd::Real y1 = py[tid1];
			rbmd::Real z1 = pz[tid1];

			rbmd::Real fs_ij, fcs_ij;


			//rs 
			for (rbmd::Id j = start_id[tid1]; j < end_id[tid1]; ++j)
			{

				rbmd::Id tid2 = id_verletlist[j];
				rbmd::Id typej = atoms_type[tid2];
				rbmd::Id molecular_id_j = molecular_type[tid2];
				rbmd::Real eps_j = eps[typej];
				rbmd::Real sigma_j = sigma[typej];

				//mix
				rbmd::Real eps_ij = sqrt(eps_i * eps_j);
				rbmd::Real sigma_ij = (sigma_i + sigma_j) / 2;

				rbmd::Real x2 = px[tid2];
				rbmd::Real y2 = py[tid2];
				rbmd::Real z2 = pz[tid2];
				rbmd::Real px12 = x2 - x1;
				rbmd::Real py12 = y2 - y1;
				rbmd::Real pz12 = z2 - z1;
				//if (molecular_id_i == molecular_id_j)
					//continue; 

				MinImageDistance(box, px12, py12, pz12);

				//compute the force_rs
				lj126_rs(rs, px12, py12, pz12, eps_ij, sigma_ij, fs_ij);

				sum_fsx += fs_ij * px12;
				sum_fsy += fs_ij * py12;
				sum_fsz += fs_ij * pz12;
			}

			//rcs
			rbmd::Id real_random_num = random_neighbor_num[tid1];
			for (rbmd::Id jj = 0; jj < real_random_num; ++jj)
			{

				rbmd::Id tid2 = id_random_neighbor[tid1 * neighbor_sample_num + jj];
				rbmd::Id typej = atoms_type[tid2];
				rbmd::Id molecular_id_j = molecular_type[tid2];
				rbmd::Real eps_j = eps[typej];
				rbmd::Real sigma_j = sigma[typej];

				//mix
				rbmd::Real eps_ij = sqrt(eps_i * eps_j);
				rbmd::Real sigma_ij = (sigma_i + sigma_j) / 2;

				rbmd::Real x2 = px[tid2];
				rbmd::Real y2 = py[tid2];
				rbmd::Real z2 = pz[tid2];
				rbmd::Real px12 = x2 - x1;
				rbmd::Real py12 = y2 - y1;
				rbmd::Real pz12 = z2 - z1;
				//if (molecular_id_i == molecular_id_j)
					//continue; 

				MinImageDistance(box, px12, py12, pz12);

				//compute the force_rcs
				lj126_rcs(rc, rs, pice_num, px12, py12, pz12, eps_ij, sigma_ij, fcs_ij); 

				sum_fcsx += fcs_ij * px12;
				sum_fcsy += fcs_ij * py12;
				sum_fcsz += fcs_ij * pz12;
			}

			//total force = fs + fcs 
			sum_fx = sum_fsx + sum_fcsx;
			sum_fy = sum_fsy + sum_fcsy;
			sum_fz = sum_fsz + sum_fcsz;

			fx[tid1] = sum_fx;
			fy[tid1] = sum_fy;
			fz[tid1] = sum_fz;

			//printf("--------------corr_force_x[tid1]:%f--,corr_force_y[tid1]:%f--,corr_force_z[tid1]:%f-\n",corr_force_x[tid1],corr_force_y[tid1],corr_force_z[tid1]);
		}


	}

	//RBL: Fix LJForce
	__global__ void FixLJRBLForce(
		const rbmd::Id num_atoms,
		const rbmd::Real corr_value_x,
		const rbmd::Real corr_value_y,
		const rbmd::Real corr_value_z,
		rbmd::Real* fx,
		rbmd::Real* fy,
		rbmd::Real* fz)
	{
		unsigned int tid1 = blockIdx.x * blockDim.x + threadIdx.x;
		if (tid1 < num_atoms)
		{
			fx[tid1] = fx[tid1] - corr_value_x;
			fy[tid1] = fy[tid1] - corr_value_y;
			fz[tid1] = fz[tid1] - corr_value_z;
		}
	}

	//verlet-list: LJEnergy
	__global__ void ComputeLJEnergy(
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
		rbmd::Real* total_evdwl)
	{
		rbmd::Real sum_elj = 0;

		unsigned int tid1 = blockIdx.x * blockDim.x + threadIdx.x;
		if (tid1 < num_atoms)
		{
			rbmd::Id typei = atoms_type[tid1];
			rbmd::Id molecular_id_i = molecular_type[tid1];
			rbmd::Real eps_i = eps[typei];
			rbmd::Real sigma_i = sigma[typei];

			rbmd::Real x1 = px[tid1];
			rbmd::Real y1 = py[tid1];
			rbmd::Real z1 = pz[tid1];

			for (int j = start_id[tid1]; j < end_id[tid1]; ++j)
			{

				rbmd::Id tid2 = id_verletlist[j];
				rbmd::Id typej = atoms_type[tid2];
				rbmd::Id molecular_id_j = molecular_type[tid2];
				rbmd::Real eps_j = eps[typej];
				rbmd::Real sigma_j = sigma[typej];

				//mix
				rbmd::Real eps_ij = sqrt(eps_i * eps_j);
				rbmd::Real sigma_ij = (sigma_i + sigma_j) / 2;

				rbmd::Real x2 = px[tid2];
				rbmd::Real y2 = py[tid2];
				rbmd::Real z2 = pz[tid2];
				rbmd::Real px12 = x2 - x1;
				rbmd::Real py12 = y2 - y1;
				rbmd::Real pz12 = z2 - z1;
				//if (molecular_id_i == molecular_id_j)
					//continue; 

				MinImageDistance(box, px12, py12, pz12);

				rbmd::Real force_lj;
				rbmd::Real energy_ij;

				lj126(cut_off, px12, py12, pz12, eps_ij, sigma_ij, force_lj, energy_ij);
				sum_elj += energy_ij;
			}

			atomicAdd(total_evdwl, sum_elj);
			//printf("--------test---evdwl[tid1]:%f---\n",evdwl[tid1]);
		}

	}

//LJCoulCutForce
        __global__ void ComputeLJCoulCutForce(
            Box* box,
            const rbmd::Real cut_off,
            const rbmd::Id num_atoms,
            const rbmd::Real alpha,
            const rbmd::Real qqr2e,
            const rbmd::Id* atoms_type,
            const rbmd::Id* molecular_type,
            const rbmd::Real* sigma,
            const rbmd::Real* eps,
            const rbmd::Id* start_id,
            const rbmd::Id* end_id,
            const rbmd::Id* id_verletlist,
            const rbmd::Real* charge,
            const rbmd::Real* px,
            const rbmd::Real* py,
            const rbmd::Real* pz,
            rbmd::Real* fx,
            rbmd::Real* fy,
            rbmd::Real* fz,
            rbmd::Real* total_evdwl,
            rbmd::Real* total_ecoul)
        {
          rbmd::Real sum_fx = 0;
          rbmd::Real sum_fy = 0;
          rbmd::Real sum_fz = 0;

          rbmd::Real sum_elj = 0;
          rbmd::Real sum_ecoul = 0;

          unsigned int tid1 = blockIdx.x * blockDim.x + threadIdx.x;
          if (tid1 < num_atoms)
          {
            rbmd::Id typei = atoms_type[tid1];
            rbmd::Id molecular_id_i =  molecular_type[tid1];
            rbmd::Real eps_i = eps[typei];
            rbmd::Real sigma_i = sigma[typei];
            rbmd::Real charge_i = charge[tid1];

            rbmd::Real x1 = px[tid1];
            rbmd::Real y1 = py[tid1];
            rbmd::Real z1 = pz[tid1];


            for (int j = start_id[tid1]; j < end_id[tid1]; ++j)
            {

              rbmd::Id tid2 = id_verletlist[j];
              rbmd::Id typej = atoms_type[tid2];
              rbmd::Id molecular_id_j = molecular_type[tid2];
              rbmd::Real eps_j = eps[typej];
              rbmd::Real sigma_j = sigma[typej];
              rbmd::Real charge_j = charge[tid2];

              //mix
              rbmd::Real eps_ij = sqrt(eps_i * eps_j);
              rbmd::Real sigma_ij = (sigma_i + sigma_j) / 2;

              rbmd::Real x2 = px[tid2];
              rbmd::Real y2 = py[tid2];
              rbmd::Real z2 = pz[tid2];
              rbmd::Real px12 = x2 - x1;
              rbmd::Real py12 = y2 - y1;
              rbmd::Real pz12 = z2 - z1;
              //if (molecular_id_i == molecular_id_j)
              //continue;

              MinImageDistance(box, px12, py12, pz12);

              rbmd::Real force_lj,force_coul,force_pair;
              rbmd::Real energy_lj,energy_coul;

              //lj cut
              lj126(cut_off, px12, py12, pz12, eps_ij, sigma_ij,force_lj,energy_lj);

              //Coul cut
              CoulCutForce(cut_off, alpha,qqr2e, charge_i, charge_j, px12, py12, pz12, force_coul,energy_coul);

              force_pair = force_lj + force_coul;

              sum_fx += force_pair * px12;
              sum_fy += force_pair * py12;
              sum_fz += force_pair * pz12;

              sum_elj += energy_lj;
              sum_ecoul += energy_coul;
            }

            fx[tid1] = sum_fx;
            fy[tid1] = sum_fy;
            fz[tid1] = sum_fz;

            atomicAdd(total_evdwl, sum_elj);
            atomicAdd(total_ecoul, sum_ecoul);
            //printf("--------test---evdwl[tid1]:%f---\n",evdwl[tid1]);
          }
        }


	//StructureFactor
	__global__ void ComputeChargeStructureFactorComponent(
		const rbmd::Id num_atoms,
		Real3 K,
		const rbmd::Real* charge,
		const rbmd::Real* px,
		const rbmd::Real* py,
		const rbmd::Real* pz,
		rbmd::Real* density_real,
		rbmd::Real* density_imag)
	{
		unsigned int tid1 = blockIdx.x * blockDim.x + threadIdx.x;
		if (tid1 < num_atoms)
		{
			rbmd::Real local_charge = charge[tid1];
			rbmd::Real dot_product = K.x * px[tid1] + K.y * py[tid1] + K.z * pz[tid1];

			density_real[tid1] = local_charge * cos(dot_product);
			density_imag[tid1] = local_charge * sin(dot_product);
		}
	}

	//EwaldForce
	__global__ void ComputeEwaldForce(
		Box* box,
		const rbmd::Id num_atoms,
	        const rbmd::Id  Kmax,
		const rbmd::Real alpha,
	        const rbmd::Real* real_array,
		const rbmd::Real* imag_array,
		const rbmd::Real* charge,
		const rbmd::Real* px,
		const rbmd::Real* py,
		const rbmd::Real* pz,
		rbmd::Real* fx,
		rbmd::Real* fy,
		rbmd::Real* fz)
	{
		rbmd::Real sum_fx = 0;
		rbmd::Real sum_fy = 0;
		rbmd::Real sum_fz = 0;

		unsigned int tid1 = blockIdx.x * blockDim.x + threadIdx.x;

		if (tid1 < num_atoms)
		{
			rbmd::Real p_x = px[tid1];
			rbmd::Real p_y = py[tid1];
			rbmd::Real p_z = pz[tid1];
			rbmd::Real charge_i = charge[tid1];


			rbmd::Id indexEwald = 0;
			for (rbmd::Id i = -Kmax; i <= Kmax; ++i)
			{
				for (rbmd::Id j = -Kmax; j <= Kmax; ++j)
				{
					for (rbmd::Id k = -Kmax; k <= Kmax; ++k)
					{
						if (i != 0 || j != 0 || k != 0)
						{
							indexEwald++;            //   total = (2*Kmax+1) * (2*Kmax+1)* (2*Kmax+1) -1
							Real3 M = { (rbmd::Real)i, (rbmd::Real)j, (rbmd::Real)k };

							rbmd::Real rhok_real_i = real_array[indexEwald -1];
							rbmd::Real rhok_imag_i = imag_array[indexEwald -1];

							 rbmd::Real force_Ewald_x, force_Ewald_y, force_Ewald_z;
							 EwaldForce(box, alpha,M, rhok_real_i, rhok_imag_i,charge_i, p_x, p_y, p_z,
										force_Ewald_x, force_Ewald_y, force_Ewald_z);

							 sum_fx += force_Ewald_x;
							 sum_fy += force_Ewald_y;
							 sum_fz += force_Ewald_z;
						}
					}
				}
			}

			fx[tid1] = sum_fx;
			fy[tid1] = sum_fy;
			fz[tid1] = sum_fz;
	    }
	}

        __global__  void ComputeBondOrder(
            Box* box,
            const rbmd::Id num_atoms,
            const rbmd::Id* atoms_type,
            const rbmd::Id* start_id,
            const rbmd::Id* end_id,
            const rbmd::Id* id_verletlist,
            const rbmd::Real* px,
            const rbmd::Real* py,
            const rbmd::Real* pz,
            rbmd::Real* fx,
            rbmd::Real* fy,
            rbmd::Real* fz)
        {
          unsigned int tid1 = blockIdx.x * blockDim.x + threadIdx.x;
          if(tid1<num_atoms)
          {
            rbmd::Id  type1 = atoms_type[tid1];
            rbmd::Real x1 = px[tid1];
            rbmd::Real y1 = py[tid1];
            rbmd::Real z1 = pz[tid1];
            for(int j1=start_id[tid1];j1<end_id[tid1];++j1)
            {
                rbmd::Id tid2 = id_verletlist[j1];
                rbmd::Id type2 = atoms_type[tid2];
                rbmd::Real x2 = px[tid2];
                rbmd::Real y2 = py[tid2];
                rbmd::Real z2 = pz[tid2];

                rbmd::Real x12 = x2-x1;
                rbmd::Real y12 = y2-y1;
                rbmd::Real z12 = z2-z1;

                MinImageDistance(box,x12,y12,z12);
                rbmd::Real  d12 = sqrt(x12*x12+y12*y12+z12*z12);
                for(int j2=start_id[tid1];j2<end_id[tid1];++j2)
                {
                  rbmd::Id tid3 = id_verletlist[j2];
                  rbmd::Id type3 = atoms_type[tid3];
                  rbmd::Real x3 = px[tid3];
                  rbmd::Real y3 = py[tid3];
                  rbmd::Real z3 = pz[tid3];

                  rbmd::Real x13  = x3 -x1;
                  rbmd::Real y13  = y3 -y1;
                  rbmd::Real z13  = z3 -z1;
                  MinImageDistance(box,x13,y13,z13);
                  rbmd::Real  d13 = sqrt(x13*x13+y13*y13+z13*z13);
                }

            }


          }

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





    //verlet-list: LJForce
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
            rbmd::Real* fx,
            rbmd::Real* fy,
            rbmd::Real* fz,
            rbmd::Real* total_evdwl)
        {
          unsigned int blocks_per_grid = (num_atoms + BLOCK_SIZE - 1) / BLOCK_SIZE;

          CHECK_KERNEL(ComputeLJForce <<<blocks_per_grid, BLOCK_SIZE, 0, 0 >>> (box, cut_off, num_atoms, atoms_type, molecular_type,
                                                                             sigma, eps, start_id, end_id, id_verletlist, px, py, pz,
                                                                             fx, fy, fz,total_evdwl));
        }

	void LJForceVirialOp<device::DEVICE_GPU>::operator()(
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
			            rbmd::Real* fx,
			            rbmd::Real* fy,
			            rbmd::Real* fz,
		                    rbmd::Real* virial_xx,
		                    rbmd::Real* virial_yy,
		                    rbmd::Real* virial_zz,
		                    rbmd::Real* virial_xy,
		                    rbmd::Real* virial_xz,
		                    rbmd::Real* virial_yz,
				    rbmd::Real* total_evdwl)
		{
		    unsigned int blocks_per_grid = (num_atoms + BLOCK_SIZE - 1) / BLOCK_SIZE;

		    CHECK_KERNEL(ComputeLJForceViral <<<blocks_per_grid, BLOCK_SIZE, 0, 0 >>> (box, cut_off, num_atoms, atoms_type, molecular_type,
				sigma, eps, start_id, end_id, id_verletlist, px, py, pz,
				fx, fy, fz, virial_xx,virial_yy, virial_zz, virial_xy, virial_xz, virial_yz,total_evdwl));
		}

	//RBL:  LJLForce
	void LJRBLForceOp<device::DEVICE_GPU>::operator()(
		Box* box,
		const rbmd::Real rs,
		const rbmd::Real rc,
		const rbmd::Id num_atoms,
		const rbmd::Id neighbor_sample_num,
		const rbmd::Id pice_num,
		const rbmd::Id* atoms_type,
		const rbmd::Id* molecular_type,
		const rbmd::Real* sigma,
		const rbmd::Real* eps,
		const rbmd::Id* start_id,
		const rbmd::Id* end_id,
		const rbmd::Id* id_verletlist,
		const rbmd::Id* id_random_neighbor,
		const rbmd::Id* random_neighbor_num,
		const rbmd::Real* px,
		const rbmd::Real* py,
		const rbmd::Real* pz,
		rbmd::Real* fx,
		rbmd::Real* fy,
		rbmd::Real* fz)
	{
		unsigned int blocks_per_grid = (num_atoms + BLOCK_SIZE - 1) / BLOCK_SIZE;

		CHECK_KERNEL(ComputeLJRBLForce <<<blocks_per_grid, BLOCK_SIZE, 0, 0 >>> 
			(box, rs, rc, num_atoms, neighbor_sample_num,pice_num,
			atoms_type, molecular_type,sigma, eps, 
			start_id, end_id, id_verletlist, id_random_neighbor, random_neighbor_num,
			px, py, pz,fx, fy, fz));
	}

	//RBL: Fix LJForce
	void FixLJRBLForceOp<device::DEVICE_GPU>::operator()(
		const rbmd::Id num_atoms,
		const rbmd::Real corr_value_x,
		const rbmd::Real corr_value_y,
		const rbmd::Real corr_value_z,
		rbmd::Real* fx,
		rbmd::Real* fy,
		rbmd::Real* fz)
	{
		unsigned int blocks_per_grid = (num_atoms + BLOCK_SIZE - 1) / BLOCK_SIZE;

		CHECK_KERNEL(FixLJRBLForce <<<blocks_per_grid, BLOCK_SIZE, 0, 0 >>>
			(num_atoms, corr_value_x, corr_value_y, corr_value_z, fx, fy, fz));
	}

	//verlet-list: LJEnergy
	void LJEnergyOp<device::DEVICE_GPU>::operator()(
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
		rbmd::Real* total_evdwl)
	{
		unsigned int blocks_per_grid = (num_atoms + BLOCK_SIZE - 1) / BLOCK_SIZE;

		CHECK_KERNEL(ComputeLJEnergy <<<blocks_per_grid, BLOCK_SIZE, 0, 0 >>> (box, cut_off, num_atoms, atoms_type, molecular_type,
			sigma, eps, start_id, end_id, id_verletlist, px, py, pz, total_evdwl));
	}

        void LJCoulCutForceOp<device::DEVICE_GPU>::operator()(
            Box* box,
            const rbmd::Real cut_off,
            const rbmd::Id num_atoms,
            const rbmd::Real alpha,
            const rbmd::Real qqr2e,
            const rbmd::Id* atoms_type,
            const rbmd::Id* molecular_type,
            const rbmd::Real* sigma,
            const rbmd::Real* eps,
            const rbmd::Id* start_id,
            const rbmd::Id* end_id,
            const rbmd::Id* id_verletlist,
            const rbmd::Real* charge,
            const rbmd::Real* px,
            const rbmd::Real* py,
            const rbmd::Real* pz,
            rbmd::Real* fx,
            rbmd::Real* fy,
            rbmd::Real* fz,
            rbmd::Real* total_evdwl,
            rbmd::Real* total_ecoul)
        {
          unsigned int blocks_per_grid = (num_atoms + BLOCK_SIZE - 1) / BLOCK_SIZE;

          CHECK_KERNEL(ComputeLJCoulCutForce <<<blocks_per_grid, BLOCK_SIZE, 0, 0 >>> (box, cut_off, num_atoms, alpha,qqr2e,atoms_type, molecular_type,
                                                                             sigma, eps, start_id, end_id, id_verletlist, charge,px, py, pz,
                                                                             fx, fy, fz,total_evdwl,total_ecoul));
        }

	//StructureFactor
	void ComputeChargeStructureFactorComponentOp<device::DEVICE_GPU>::operator()(
		const rbmd::Id num_atoms,
		const Real3 K,
		const rbmd::Real* charge,
		const rbmd::Real* px,
		const rbmd::Real* py,
		const rbmd::Real* pz,
		rbmd::Real* density_real,
		rbmd::Real* density_imag)
	{
		unsigned int blocks_per_grid = (num_atoms + BLOCK_SIZE - 1) / BLOCK_SIZE;
		
		CHECK_KERNEL(ComputeChargeStructureFactorComponent <<<blocks_per_grid, BLOCK_SIZE, 0, 0 >>>
			(num_atoms, K, charge, px, py ,pz , density_real, density_imag));

	}

	//EwaldForce
	void ComputeEwaldForceOp<device::DEVICE_GPU>::operator()(
		Box* box,
		const rbmd::Id num_atoms,
		const rbmd::Id  Kmax,
		const rbmd::Real alpha,
		const rbmd::Real* real_array,
		const rbmd::Real* imag_array,
		const rbmd::Real* charge,
		const rbmd::Real* px,
		const rbmd::Real* py,
		const rbmd::Real* pz,
		rbmd::Real* fx,
                rbmd::Real* fy,
                rbmd::Real* fz)
	{
		unsigned int blocks_per_grid = (num_atoms + BLOCK_SIZE - 1) / BLOCK_SIZE;

		CHECK_KERNEL(ComputeEwaldForce <<<blocks_per_grid, BLOCK_SIZE, 0, 0 >>>
			(box,num_atoms, Kmax, alpha, real_array, imag_array, charge,
		     px, py, pz, fx, fy, fz));

	}

}

