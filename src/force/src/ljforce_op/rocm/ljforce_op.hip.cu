#include <hip/hip_runtime.h>

#include "../common/rbmd_define.h"
#include "ljforce_op.h"
#include "model/box.h"

namespace op
{

	//---------device---------//
        __device__ rbmd::Id index2D(rbmd::Id row, rbmd::Id col, rbmd::Id num_columns) {
          return row * num_columns + col;
        }

        __device__ rbmd::Id index3D(
          rbmd::Id depth,
          rbmd::Id row,
          rbmd::Id col,
          rbmd::Id num_rows,
          rbmd::Id num_columns) {
          return depth * (num_rows * num_columns) + row * num_columns + col;
        }

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
		const rbmd::Real dis_2 = px12 * px12 + py12 * py12 + pz12 * pz12;
		const rbmd::Real cut_off_2 = cut_off * cut_off;

		if (dis_2 < cut_off_2 && dis_2 > EPSILON)
		{
			rbmd::Real sigmaij_6 = POW(sigma_ij, 6.0);
			rbmd::Real dis_6 = POW(dis_2, 3.0);
			rbmd::Real sigmaij_dis_6 = sigmaij_6 / dis_6;

			force_lj = -24 * eps_ij * ((2 * sigmaij_dis_6 - 1) * sigmaij_dis_6) / dis_2;
			energy_lj = 0.5 * (4 * eps_ij * (sigmaij_6 / dis_6 - 1) * (sigmaij_6 / dis_6));
		}
		else
		{
		  force_lj  = 0.0;
		  energy_lj = 0.0;
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
		const rbmd::Real dis_2 = px12 * px12 + py12 * py12 + pz12 * pz12;
		const rbmd::Real rs_2 = rs * rs;

		if (dis_2 < rs_2 && dis_2 > EPSILON)
		{
		   rbmd::Real sigmaij_6 = POW(sigma_ij, 6.0);
		   rbmd::Real dis_6 = POW(dis_2, 3.0);
		   rbmd::Real sigmaij_dis_6 = sigmaij_6 / dis_6;
		   fs_ij = -24 * eps_ij * ((2 * sigmaij_dis_6 - 1) * sigmaij_dis_6) / dis_2;
		}
                else fs_ij = 0.0;

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
		const rbmd::Real dis_2 = px12 * px12 + py12 * py12 + pz12 * pz12;
		const rbmd::Real rc_2 = rc * rc;
		const rbmd::Real rs_2 = rs * rs;

		if (dis_2 < rc_2 && dis_2 > rs_2)
		{
			rbmd::Real sigmaij_6 = POW(sigma_ij, 6.0);
			rbmd::Real dis_6 = POW(dis_2, 3.0);
			rbmd::Real sigmaij_dis_6 = sigmaij_6 / dis_6;

			fcs_ij = pice_num * (-24 * eps_ij * ((2 * sigmaij_dis_6 - 1) * sigmaij_dis_6) / dis_2);
		}
		else fcs_ij = 0.0;
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
			rbmd::Real& force_coul,
			rbmd::Real& energy_coul)
	{
		const rbmd::Real dis_2 = px12 * px12 + py12 * py12 + pz12 * pz12;
		const rbmd::Real dis = SQRT(dis_2);
		const rbmd::Real cut_off_2 = cut_off * cut_off;

		if (dis_2 < cut_off_2 && dis_2 > EPSILON)
		{
			rbmd::Real erfcx = SQRT(alpha) * dis;
			rbmd::Real expx = -alpha * dis_2;
			rbmd::Real gnear_value = (1.0 - ERF(erfcx)) / dis_2 +
				2 * SQRT(alpha) * EXP(expx) / (SQRT(M_PI) * dis);

			force_coul = qqr2e * (-charge_i * charge_j * gnear_value / dis);
			energy_coul = qqr2e * (0.5 * charge_i * charge_j * (1.0 - ERF(SQRT(alpha) * dis)) / dis);
		}
		else
		{
		  force_coul  = 0.0;
                  energy_coul = 0.0;
		}
	}

        __device__
        void CoulCutForce_erf(
                rbmd::Real cut_off,
                rbmd::Real alpha,
                rbmd::Real qqr2e,
                rbmd::Real table_pij,
                rbmd::Real charge_i,
                rbmd::Real charge_j,
                rbmd::Real px12,
                rbmd::Real py12,
                rbmd::Real pz12,
                rbmd::Real& force_coul,
                rbmd::Real& energy_coul)
	{
	  const rbmd::Real dis_2 = px12 * px12 + py12 * py12 + pz12 * pz12;
	  const rbmd::Real dis = SQRT(dis_2);
	  const rbmd::Real cut_off_2 = cut_off * cut_off;

	  if (dis_2 < cut_off_2 && dis_2 > EPSILON)
	  {
	    force_coul = qqr2e * (-charge_i * charge_j * table_pij / dis);
	    energy_coul = qqr2e * (0.5 * charge_i * charge_j * (1.0 - ERF(SQRT(alpha) * dis)) / dis);
	  }
	  else
	  {
	    force_coul = 0.0;
            energy_coul =0.0;
	  }
	}

        __device__
        void CoulCutForce_rs(
                rbmd::Real rs,
                rbmd::Real alpha,
                rbmd::Real qqr2e,
                rbmd::Real charge_i,
                rbmd::Real charge_j,
                rbmd::Real px12,
                rbmd::Real py12,
                rbmd::Real pz12,
                rbmd::Real& force_coul)
	{
	  const rbmd::Real dis_2 = px12 * px12 + py12 * py12 + pz12 * pz12;
	  const rbmd::Real dis = SQRT(dis_2);
	  const rbmd::Real rs_2 = rs * rs;

	  if (dis_2 < rs_2 && dis_2 > EPSILON)
	  {
	    rbmd::Real erfcx = SQRT(alpha) * dis;
	    rbmd::Real expx = -alpha * dis_2;
	    rbmd::Real gnear_value = (1.0 - ERF(erfcx)) / dis_2 +
                    2 * SQRT(alpha) * EXP(expx) / (SQRT(M_PI) * dis);

	    force_coul = qqr2e * (-charge_i * charge_j * gnear_value / dis);
	  }
	  else force_coul = 0.0;
	}

__device__
void CoulCutForce_rs_erf(
        rbmd::Real rs,
        rbmd::Real alpha,
        rbmd::Real qqr2e,
        rbmd::Real table_pij,
        rbmd::Real charge_i,
        rbmd::Real charge_j,
        rbmd::Real px12,
        rbmd::Real py12,
        rbmd::Real pz12,
        rbmd::Real& force_coul)
	{
	  const rbmd::Real dis_2 = px12 * px12 + py12 * py12 + pz12 * pz12;
	  const rbmd::Real dis = SQRT(dis_2);
	  const rbmd::Real rs_2 = rs * rs;

	  if (dis_2 < rs_2 && dis_2 > EPSILON)
	  {
	    force_coul = qqr2e * (-charge_i * charge_j * table_pij / dis);
	  }
	  else force_coul = 0.0;
	}

__device__
void CoulCutForce_rcs(
        rbmd::Real rc,
        rbmd::Real rs,
        rbmd::Id  pice_num,
        rbmd::Real alpha,
        rbmd::Real qqr2e,
        rbmd::Real charge_i,
        rbmd::Real charge_j,
        rbmd::Real px12,
        rbmd::Real py12,
        rbmd::Real pz12,
        rbmd::Real& force_coul)
	{
	  const rbmd::Real dis_2 = px12 * px12 + py12 * py12 + pz12 * pz12;
	  const rbmd::Real dis = SQRT(dis_2);
	  const rbmd::Real rc_2 = rc * rc;
	  const rbmd::Real rs_2 = rs * rs;

	  if (dis_2 < rc_2 && dis_2 > rs_2)
	  {
	    rbmd::Real erfcx = SQRT(alpha) * dis;
	    rbmd::Real expx = -alpha * dis_2;
	    rbmd::Real gnear_value = (1.0 - ERF(erfcx)) / dis_2 +
                    2 * SQRT(alpha) * EXP(expx) / (SQRT(M_PI) * dis);

	    force_coul = pice_num * qqr2e * (-charge_i * charge_j * gnear_value / dis);
	  }
	  else force_coul = 0.0;
	}

__device__
void CoulCutForce_rcs_erf(
        rbmd::Real rc,
        rbmd::Real rs,
        rbmd::Id  pice_num,
        rbmd::Real alpha,
        rbmd::Real qqr2e,
        rbmd::Real table_pij,
        rbmd::Real charge_i,
        rbmd::Real charge_j,
        rbmd::Real px12,
        rbmd::Real py12,
        rbmd::Real pz12,
        rbmd::Real& force_coul)
	{
	  const rbmd::Real dis_2 = px12 * px12 + py12 * py12 + pz12 * pz12;
	  const rbmd::Real dis = SQRT(dis_2);
	  const rbmd::Real rc_2 = rc * rc;
	  const rbmd::Real rs_2 = rs * rs;

	  if (dis_2 < rc_2 && dis_2 > rs_2)
	  {
	    force_coul = pice_num * qqr2e * (-charge_i * charge_j * table_pij / dis);
	  }
	  else force_coul = 0.0;
	}


	//EwaldForce
	__device__ void EwaldForce(
		Box* box,
		const rbmd::Real alpha,
		const int3 M,
		const rbmd::Real qqr2e,
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
		rbmd::Real alpha_inv = 1 / alpha;

		rbmd::Real factor_a = -4 * M_PI * charge;
		rbmd::Real factor_b = EXP(-0.25 * range_K_2 * alpha_inv);
		rbmd::Real factor_c = COS(dot_product) * rhok_imag_i;
		rbmd::Real factor_d = SIN(dot_product) * rhok_real_i;

		force_ewald = factor_a / (volume * range_K_2) * factor_b * (factor_c - factor_d);
		force_ewald *= qqr2e;

		force_ewald_x = force_ewald * K.x;
		force_ewald_y = force_ewald * K.y;
		force_ewald_z = force_ewald * K.z;
	}

       //RBEForce
       __device__ void RBEForce(
           Box* box,
           const Real3 M,
           const rbmd::Real qqr2e,
           const rbmd::Real rhok_real_i,
           const rbmd::Real rhok_imag_i,
           const rbmd::Real charge,
           const rbmd::Real px,
           const rbmd::Real py,
           const rbmd::Real pz,
           rbmd::Real& force_rbe_x,
           rbmd::Real& force_rbe_y,
           rbmd::Real& force_rbe_z)
	{
	  rbmd::Real force_rbe;
	  rbmd::Real volume = box->_length[0] * box->_length[1] * box->_length[2];
	  Real3 K = make_Real3(2 * M_PI * M.x / box->_length[0],
                  2 * M_PI * M.y / box->_length[1],
                  2 * M_PI * M.z / box->_length[2]);


	  rbmd::Real range_K_2 = K.x * K.x + K.y * K.y + K.z * K.z;
	  rbmd::Real dot_product = K.x * px + K.y * py + K.z * pz;

	  rbmd::Real factor_a = -4 * M_PI * charge;
	  rbmd::Real factor_b = COS(dot_product) * rhok_imag_i;
	  rbmd::Real factor_c = SIN(dot_product) * rhok_real_i;

	  force_rbe= (factor_a / (volume * range_K_2))* (factor_b- factor_c);
	  force_rbe *= qqr2e;

	  force_rbe_x = force_rbe * K.x;
	  force_rbe_y = force_rbe * K.y;
	  force_rbe_z = force_rbe * K.z;
	}

        __device__  void ComputeS(
            Box* box,
            const rbmd::Real alpha,
            rbmd::Real& S)
	{
	  Real3 H{0.0,0.0,0.0};
	  for (rbmd::Id i = 0; i < 3; ++i)
	  {
	    const rbmd::Real factor = -(alpha * box->_length[i] * box->_length[i]);

	    for (rbmd::Id m = -10; m <= 10; m++)
	    {
	      rbmd::Real expx = m * m * factor;
	      H.data[i] += EXP(expx);
	    }
	    H.data[i] *= SQRT(-(factor) / M_PI);
	  }

	   rbmd::Real factor_3 = H.data[0] * H.data[1] * H.data[2];
	   S = factor_3 - 1;
	}

	template<typename Func>
	__device__ void ExecuteOnKmax(const rbmd::Id& k_maxconst, Func& function)
	{
		rbmd::Id indexEwald = 0;
		for (rbmd::Id i = -k_maxconst; i <= k_maxconst; i++)
		{
			for (rbmd::Id j = -k_maxconst; j <= k_maxconst; j++)
			{
				for (rbmd::Id k = -k_maxconst; k <= k_maxconst; k++)
				{
					if (i != 0 || j != 0 || k != 0)
					{
						indexEwald++;
						int3 M = make_Int3(i, j, k);
						function(M, indexEwald);
					}
				}
			}
		}
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
		__shared__ typename hipcub::BlockReduce<rbmd::Real, BLOCK_SIZE>::TempStorage
			temp_storage;

		rbmd::Real sum_fx = 0;
		rbmd::Real sum_fy = 0;
		rbmd::Real sum_fz = 0;

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
				rbmd::Real eps_ij = SQRT(eps_i * eps_j);
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

				rbmd::Real force_lj, fpair;
				rbmd::Real energy_lj;

				lj126(cut_off, px12, py12, pz12, eps_ij, sigma_ij, force_lj, energy_lj);
				sum_fx += force_lj * px12;
				sum_fy += force_lj * py12;
				sum_fz += force_lj * pz12;

				sum_elj += energy_lj;
			}

			fx[tid1] = sum_fx;
			fy[tid1] = sum_fy;
			fz[tid1] = sum_fz;

			rbmd::Real block_sum =
				hipcub::BlockReduce<rbmd::Real, BLOCK_SIZE>(temp_storage).Sum(sum_elj);
			if (threadIdx.x == 0) {
				atomicAdd(total_evdwl, block_sum);
			}
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
				rbmd::Real eps_ij = SQRT(eps_i * eps_j);
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

				rbmd::Real force_lj, fpair;
				rbmd::Real energy_lj;

				lj126(cut_off, px12, py12, pz12, eps_ij, sigma_ij, force_lj, energy_lj);
				sum_fx += force_lj * px12;
				sum_fy += force_lj * py12;
				sum_fz += force_lj * pz12;

				fpair = -0.5 * force_lj;
				sum_virial_xx += px12 * px12 * fpair;
				sum_virial_yy += py12 * py12 * fpair;
				sum_virial_zz += pz12 * pz12 * fpair;
				sum_virial_xy += px12 * py12 * fpair;
				sum_virial_xz += px12 * pz12 * fpair;
				sum_virial_yz += py12 * pz12 * fpair;

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
				rbmd::Real eps_ij = SQRT(eps_i * eps_j);
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
				rbmd::Real eps_ij = SQRT(eps_i * eps_j);
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
		}


	}

	//RBL: Fix RBL Force
	__global__ void FixRBLForce(
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
	  __shared__ typename hipcub::BlockReduce<rbmd::Real, BLOCK_SIZE>::TempStorage
        temp_storage;
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
				rbmd::Real eps_ij = SQRT(eps_i * eps_j);
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

		  rbmd::Real block_sum = hipcub::BlockReduce<rbmd::Real, BLOCK_SIZE>
		  (temp_storage).Sum(sum_elj);
		  if (threadIdx.x == 0) {
		    atomicAdd(total_evdwl, block_sum);
		  }
		}

	}

	//LJCoulCutForce
	__global__ void ComputeLJCutCoulForce(
		Box* box,ERFTable* erf_table,
		const rbmd::Real cut_off,
		const rbmd::Id num_atoms,
		const rbmd::Real alpha,
		const rbmd::Real qqr2e,
		const rbmd::Id* atoms_type,
		const rbmd::Id* molecular_id,
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
	  __shared__ typename hipcub::BlockReduce<rbmd::Real, BLOCK_SIZE>::TempStorage
        temp_storage_elj;
	  __shared__ typename hipcub::BlockReduce<rbmd::Real, BLOCK_SIZE>::TempStorage
        temp_storage_ecoul;
		rbmd::Real sum_fx = 0;
		rbmd::Real sum_fy = 0;
		rbmd::Real sum_fz = 0;

		rbmd::Real sum_elj = 0;
		rbmd::Real sum_ecoul = 0;

		unsigned int tid1 = blockIdx.x * blockDim.x + threadIdx.x;
		if (tid1 < num_atoms)
		{
			rbmd::Id typei = atoms_type[tid1];
			rbmd::Id molecular_id_i = molecular_id[tid1];
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
				rbmd::Id molecular_id_j = molecular_id[tid2];
				rbmd::Real eps_j = eps[typej];
				rbmd::Real sigma_j = sigma[typej];
				rbmd::Real charge_j = charge[tid2];

				//mix
				rbmd::Real eps_ij = SQRT(eps_i * eps_j);
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

			        //erf
			        rbmd::Real  dis = SQRT(px12*px12+py12*py12+pz12*pz12);
			        rbmd::Id index_table_pij = erf_table->Extract(dis);
			        rbmd::Real  table_pij = erf_table->TableGnearValue(dis,index_table_pij);

				rbmd::Real force_lj, force_coul, force_pair;
				rbmd::Real energy_lj, energy_coul;

				//lj cut
				lj126(cut_off, px12, py12, pz12, eps_ij, sigma_ij, force_lj, energy_lj);

				//Coul cut
				CoulCutForce(cut_off, alpha, qqr2e, charge_i, charge_j, px12, py12, pz12, force_coul, energy_coul);

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
			//printf("--------test---fx[tid1]:%f---\n",fx[tid1]);
		  rbmd::Real block_sum_elj =
                          hipcub::BlockReduce<rbmd::Real, BLOCK_SIZE>(temp_storage_elj).Sum(sum_elj);
		  rbmd::Real block_sum_ecoul =
                         hipcub::BlockReduce<rbmd::Real, BLOCK_SIZE>(temp_storage_ecoul).Sum(sum_ecoul);

		  if (threadIdx.x == 0) {
		    atomicAdd(total_evdwl, block_sum_elj);
		    atomicAdd(total_ecoul, block_sum_ecoul);
		  }
		}
	}

       //LJCutCoulEnergy
	__global__ void ComputeLJCutCoulEnergy(
		Box* box,ERFTable* erf_table,
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
		rbmd::Real* total_evdwl,
		rbmd::Real* total_ecoul)
	{
	  __shared__ typename hipcub::BlockReduce<rbmd::Real, BLOCK_SIZE>::TempStorage
        temp_storage_elj;
	  __shared__ typename hipcub::BlockReduce<rbmd::Real, BLOCK_SIZE>::TempStorage
        temp_storage_ecoul;

		rbmd::Real sum_elj = 0;
		rbmd::Real sum_ecoul = 0;

		unsigned int tid1 = blockIdx.x * blockDim.x + threadIdx.x;
		if (tid1 < num_atoms)
		{
			rbmd::Id typei = atoms_type[tid1];
			rbmd::Id molecular_id_i = molecular_type[tid1];
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
				rbmd::Real eps_ij = SQRT(eps_i * eps_j);
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

			        //erf
			        rbmd::Real  dis = SQRT(px12*px12+py12*py12+pz12*pz12);
			        rbmd::Id index_table_pij = erf_table->Extract(dis);
			        rbmd::Real  table_pij = erf_table->TableGnearValue(dis,index_table_pij);

			        rbmd::Real force_lj, force_coul, force_pair;
				rbmd::Real energy_lj, energy_coul;

				//lj cut
				lj126(cut_off, px12, py12, pz12, eps_ij, sigma_ij,
				  force_lj, energy_lj);

				//Coul cut
				CoulCutForce_erf(cut_off, alpha, qqr2e,table_pij, charge_i, charge_j,
				  px12, py12, pz12, force_coul, energy_coul);

				sum_elj += energy_lj;
				sum_ecoul += energy_coul;
			}

		  rbmd::Real block_sum_elj =
                          hipcub::BlockReduce<rbmd::Real, BLOCK_SIZE>(temp_storage_elj).Sum(sum_elj);
		  rbmd::Real block_sum_ecoul =
                         hipcub::BlockReduce<rbmd::Real, BLOCK_SIZE>(temp_storage_ecoul).Sum(sum_ecoul);

		  if (threadIdx.x == 0) {
		    atomicAdd(total_evdwl, block_sum_elj);
		    atomicAdd(total_ecoul, block_sum_ecoul);
		  }
		}
	}

       //LJCutCoul RBL
	__global__ void ComputeLJCutCoulRBLForce(
		Box* box,
		ERFTable* erf_table,
		const rbmd::Real rs,
                const rbmd::Real rc,
		const rbmd::Id num_atoms,
		const rbmd::Id neighbor_sample_num,
                const rbmd::Id pice_num,
		const rbmd::Real alpha,
		const rbmd::Real qqr2e,
		const rbmd::Id* atoms_type,
		const rbmd::Id* molecular_type,
		const rbmd::Real* sigma,
		const rbmd::Real* eps,
		const rbmd::Id* start_id,
		const rbmd::Id* end_id,
		const rbmd::Id* id_verletlist,
		const rbmd::Id* id_random_neighbor,
                const rbmd::Id* random_neighbor_num,
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
		        rbmd::Real charge_i = charge[tid1];

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
			  rbmd::Real charge_j = charge[tid2];

			  //mix
			  rbmd::Real eps_ij = SQRT(eps_i * eps_j);
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

			  //erf
			  rbmd::Real  dis = SQRT(px12*px12+py12*py12+pz12*pz12);
			  rbmd::Id index_table_pij = erf_table->Extract(dis);
			  rbmd::Real  table_pij = erf_table->TableGnearValue(dis,index_table_pij);

			  //compute the force_rs
			  rbmd::Real force_lj_rs, force_coul_rs, force_pair;
			  lj126_rs(rs, px12, py12, pz12, eps_ij, sigma_ij,
			    force_lj_rs);
			  CoulCutForce_rs_erf(rs, alpha, qqr2e,table_pij,charge_i, charge_j,
			    px12, py12, pz12, force_coul_rs);

			  fs_ij = force_lj_rs + force_coul_rs;
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
			  rbmd::Real charge_j = charge[tid2];

			  //mix
			  rbmd::Real eps_ij = SQRT(eps_i * eps_j);
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

			  //erf
			  rbmd::Real  dis = SQRT(px12*px12+py12*py12+pz12*pz12);
			  rbmd::Id index_table_pij = erf_table->Extract(dis);
			  rbmd::Real  table_pij = erf_table->TableGnearValue(dis,index_table_pij);

			  //compute the force_rcs
			  rbmd::Real force_lj_rcs, force_coul_rcs, force_pair;
			  lj126_rcs(rc, rs, pice_num, px12, py12, pz12,
			    eps_ij, sigma_ij, force_lj_rcs);
			  CoulCutForce_rcs_erf(rc,rs, pice_num, alpha, qqr2e,table_pij,
			    charge_i, charge_j, px12, py12, pz12, force_coul_rcs);
			  fcs_ij = force_lj_rcs+force_coul_rcs;

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
		}
	}

	//StructureFactor
	__global__ void ComputeChargeStructureFactor(
		const rbmd::Id num_atoms,
		const Real3 K,
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

			density_real[tid1] = local_charge * COS(dot_product);
			density_imag[tid1] = local_charge * SIN(dot_product);
		}
	}

        //Charge Structure  Factor on Pnumber
      __global__ void  ComputePnumberChargeStructureFactor(
                        Box* box,
        		const rbmd::Id num_atoms,
         		const rbmd::Id p_number,
         		const rbmd::Real* charge,
         		const rbmd::Real* p_sample_x,
         		const rbmd::Real* p_sample_y,
         		const rbmd::Real* p_sample_z,
        		const rbmd::Real* px,
                        const rbmd::Real* py,
                        const rbmd::Real* pz,
                        rbmd::Real* density_real,
                        rbmd::Real* density_imag)
        {
	  unsigned int tid1 = blockIdx.x * blockDim.x + threadIdx.x;
	  if (tid1 < num_atoms)
	  {
	    rbmd::Real chargei = charge[tid1];
	    rbmd::Real p_x = px[tid1];
	    rbmd::Real p_y = py[tid1];
	    rbmd::Real p_z = pz[tid1];

	    for (rbmd::Id i = 0; i < p_number; i++)
	    {
	      rbmd::Id index = tid1 + i * num_atoms;

	      rbmd::Real  k_x = p_sample_x[i];
	      rbmd::Real  k_y = p_sample_y[i];
	      rbmd::Real  k_z = p_sample_z[i];
	      k_x = 2 * M_PI * k_x  / box->_length[0];
	      k_y = 2 * M_PI * k_y  / box->_length[1];
	      k_z = 2 * M_PI * k_z  / box->_length[2];

	      rbmd::Real dot_product = k_x*p_x + k_y*p_y + k_z*p_z;
	      density_real[index] = chargei * COS(dot_product);
	      density_imag[index] = chargei * SIN(dot_product);
	    }
	  }
	}

	//EwaldForce
	__global__ void ComputeEwaldForce(
		Box* box,
		const rbmd::Id num_atoms,
		const rbmd::Id  Kmax,
		const rbmd::Real alpha,
		const rbmd::Real qqr2e,
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

			auto function = [&](const int3& M, const rbmd::Id& indexEwald)
			{
				const rbmd::Real rhok_real_i = real_array[indexEwald - 1];
				const rbmd::Real rhok_imag_i = imag_array[indexEwald - 1];

				rbmd::Real force_Ewald_x, force_Ewald_y, force_Ewald_z;
				EwaldForce(box, alpha, M, qqr2e, rhok_real_i, rhok_imag_i, charge_i,
					p_x, p_y, p_z, force_Ewald_x, force_Ewald_y, force_Ewald_z);

				sum_fx += force_Ewald_x;
				sum_fy += force_Ewald_y;
				sum_fz += force_Ewald_z;
			};
			ExecuteOnKmax(Kmax, function);

			fx[tid1] = sum_fx;
			fy[tid1] = sum_fy;
			fz[tid1] = sum_fz;
		}
	}

	__global__ void ComputeSqCharge(
		const rbmd::Id num_atoms,
		const rbmd::Real* charge,
		rbmd::Real* sq_charge)
	{
		unsigned int tid1 = blockIdx.x * blockDim.x + threadIdx.x;
		if (tid1 < num_atoms)
		{
			rbmd::Real chargei = charge[tid1];
			sq_charge[tid1] = chargei * chargei;
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
		if (tid1 < num_atoms)
		{
			rbmd::Id  type1 = atoms_type[tid1];
			rbmd::Real x1 = px[tid1];
			rbmd::Real y1 = py[tid1];
			rbmd::Real z1 = pz[tid1];
			for (int j1 = start_id[tid1]; j1 < end_id[tid1]; ++j1)
			{
				rbmd::Id tid2 = id_verletlist[j1];
				rbmd::Id type2 = atoms_type[tid2];
				rbmd::Real x2 = px[tid2];
				rbmd::Real y2 = py[tid2];
				rbmd::Real z2 = pz[tid2];

				rbmd::Real x12 = x2 - x1;
				rbmd::Real y12 = y2 - y1;
				rbmd::Real z12 = z2 - z1;

				MinImageDistance(box, x12, y12, z12);
				rbmd::Real  d12 = SQRT(x12 * x12 + y12 * y12 + z12 * z12);
				for (int j2 = start_id[tid1]; j2 < end_id[tid1]; ++j2)
				{
					rbmd::Id tid3 = id_verletlist[j2];
					rbmd::Id type3 = atoms_type[tid3];
					rbmd::Real x3 = px[tid3];
					rbmd::Real y3 = py[tid3];
					rbmd::Real z3 = pz[tid3];

					rbmd::Real x13 = x3 - x1;
					rbmd::Real y13 = y3 - y1;
					rbmd::Real z13 = z3 - z1;
					MinImageDistance(box, x13, y13, z13);
					rbmd::Real  d13 = SQRT(x13 * x13 + y13 * y13 + z13 * z13);
				}

			}


		}

	}

           //index
         __global__ void GenerateIndexArray(
           const rbmd::Id num_atoms,
           const rbmd::Id RBE_P,
           rbmd::Id* psample_key)
	{
	  unsigned int tid1 = blockIdx.x * blockDim.x + threadIdx.x;
	  if (tid1 < num_atoms * RBE_P)
	  {
	    psample_key[tid1] = tid1 / num_atoms;
	  }
	}

        //RBEForce
        __global__ void ComputeRBEForce(
        Box* box,
        const rbmd::Id num_atoms,
        const rbmd::Id  p_number,
        const rbmd::Real alpha,
        const rbmd::Real qqr2e,
        const rbmd::Real* real_array,
        const rbmd::Real* imag_array,
        const rbmd::Real* charge,
        const rbmd::Real* p_sample_x,
        const rbmd::Real* p_sample_y,
        const rbmd::Real* p_sample_z,
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
	    //
	    for (rbmd::Id i = 0; i < p_number; i++)
	    {
	      const Real3 M = make_Real3(p_sample_x[i],
	                          p_sample_y[i],
	                           p_sample_z[i]);

	      const rbmd::Real rhok_real_i = real_array[i];
	      const rbmd::Real rhok_imag_i = imag_array[i];

	      rbmd::Real force_rbe_x, force_rbe_y, force_rbe_z;
	      RBEForce(box, M, qqr2e, rhok_real_i, rhok_imag_i, charge_i,
                      p_x, p_y, p_z, force_rbe_x, force_rbe_y, force_rbe_z);

	      sum_fx += force_rbe_x;
	      sum_fy += force_rbe_y;
	      sum_fz += force_rbe_z;
	    }
            //
            rbmd::Real sum_gauss;
	    ComputeS(box,alpha,sum_gauss);
	    //printf("--------test---sum_gauss:%f\n",sum_gauss);

	    sum_fx = sum_fx * sum_gauss / p_number;
	    sum_fy = sum_fy * sum_gauss / p_number;
	    sum_fz = sum_fz * sum_gauss / p_number;

	    fx[tid1] = sum_fx;
	    fy[tid1] = sum_fy;
	    fz[tid1] = sum_fz;
	    //printf("--------test---force0:%f---,force1:%f---force2:%f\n",fx[tid1],fy[tid1],fz[tid1]);
	  }
	}

        //
      __global__ void ComputeSpecialCoulForce(
       Box* box,
       const rbmd::Id num_atoms,
       const rbmd::Real qqr2e,
       const rbmd::Id* atoms_id,
       const rbmd::Id*  atoms_vec,
       const rbmd::Id*  atoms_offset,
       const rbmd::Id*  atom_count,
       const rbmd::Id*  special_ids,
       const rbmd::Real*  special_weights,
       const rbmd::Id*  special_offset,
       const rbmd::Id*  special_count,
       const rbmd::Real* charge,
       const rbmd::Real* px,
       const rbmd::Real* py,
       const rbmd::Real* pz,
       rbmd::Real* fx,
       rbmd::Real* fy,
       rbmd::Real* fz,
       rbmd::Real* total_especial_coul)
       {
          __shared__ typename hipcub::BlockReduce<rbmd::Real, BLOCK_SIZE>::TempStorage
    temp_storage;

          rbmd::Real  sum_fx=0.0;
          rbmd::Real  sum_fy=0.0;
          rbmd::Real  sum_fz=0.0;
          rbmd::Real  sum_energy_special_coul=0.0;

          unsigned int tid1 = blockIdx.x * blockDim.x + threadIdx.x;
          if (tid1 < num_atoms)
          {
            rbmd::Id id = atoms_id[tid1];
            if(atom_count[id] < 3 )
            {
              sum_fx=sum_fx=sum_fx=0.0;
              sum_energy_special_coul=0.0;
            }
            else
            {
              rbmd::Real charge_i = charge[id];
              rbmd::Real x1 = px[id];
              rbmd::Real y1 = py[id];
              rbmd::Real z1 = pz[id];
              rbmd::Id  num_offset = atoms_offset[id];
              rbmd::Id  num_components = special_offset[id];
              //printf("tid1 %i id %i\n",tid1,id);
              //printf("tid1 %i id %i num_offset %i\n",tid1 ,id,num_offset);
              //printf("id %i num_offset %i\n",id ,num_offset);
              //printf("id %i num_components %i\n",id ,num_components);
              //printf(" id %i  px[id] %f  tid %i  px[tid1] %f\n", id , px[id] ,
              //tid1 ,px[tid1]);
              for (rbmd::Id j = 0; j < atom_count[id]; ++j)
              {
                rbmd::Id id2 = atoms_vec[num_offset + j];
                //printf("atom_count %i tid2 %i\n",atom_count[id] ,tid2);

                if (id == id2)
                  continue;

                rbmd::Real charge_j = charge[id2];
                rbmd::Real x2 = px[id2];
                rbmd::Real y2 = py[id2];
                rbmd::Real z2 = pz[id2];

                rbmd::Real x12 = x2 - x1;
                rbmd::Real y12 = y2 - y1;
                rbmd::Real z12 = z2 - z1;
                MinImageDistance(box, x12, y12, z12);

                rbmd::Real dis_ij = SQRT(x12*x12+y12*y12+z12*z12);
                rbmd::Real dis_ij3 = POW(dis_ij, 3.0);
                rbmd::Real force_component = -qqr2e * charge_i * charge_j / dis_ij3;
                rbmd::Real energy_atom = 0.5 * qqr2e * charge_i * charge_j / dis_ij;

                rbmd::Real weight = 1.0;
                for (rbmd::Id k = 0; k < special_count[id]; ++k)
                {
                  //printf("special_count %i tid22 %i\n",special_count[id] ,special_ids[num_components+k]);
                  if (special_ids[num_components+k] == id2)
                  {
                    weight = special_weights[num_components+k];
                    //printf("weight %f\n", weight);
                  }
                }
                sum_fx += (1.0 - weight) * force_component * x12;
                sum_fy += (1.0 - weight) * force_component * y12;
                sum_fz += (1.0 - weight) * force_component * z12;
                sum_energy_special_coul += (1.0 - weight) * energy_atom;
                //printf(" energy_special_coul %f\n" ,sum_energy_special_coul);
              }
            }

            //printf(" force %f %f %f\n" ,sum_fx,sum_fy,sum_fz);
            fx[tid1] = sum_fx;
            fy[tid1] = sum_fy;
            fz[tid1] = sum_fz;

            rbmd::Real block_sum_especial_coul = hipcub::BlockReduce<rbmd::Real, BLOCK_SIZE>
            (temp_storage).Sum(sum_energy_special_coul);

            if (threadIdx.x == 0) {
              atomicAdd(total_especial_coul, block_sum_especial_coul);
            }
         }
       }

       __global__ void ComputeBondForce(
       Box* box,
       const rbmd::Id num_bonds,
       const rbmd::Id* atom_id_to_idx,
       const rbmd::Real* bond_coeffs_k,
       const rbmd::Real* bond_coeffs_equilibrium,
       const rbmd::Id* bond_type,
       const rbmd::Id* bondlisti,
       const rbmd::Id* bondlistj,
       const rbmd::Real* px,
       const rbmd::Real* py,
       const rbmd::Real* pz,
       rbmd::Real* fx,
       rbmd::Real* fy,
       rbmd::Real* fz,
       rbmd::Id* atom_ids_out,
       rbmd::Real* energy_bond)
        {
	  __shared__ typename hipcub::BlockReduce<rbmd::Real, BLOCK_SIZE>::TempStorage
    temp_storage;

	  rbmd::Real local_energy_bond  = 0;

	  unsigned int tid1 = blockIdx.x * blockDim.x + threadIdx.x;
	  if (tid1 < num_bonds)
	  {
	    rbmd::Id bondi = bondlisti[tid1];
	    rbmd::Id bondj = bondlistj[tid1];
	    rbmd::Id bondii = atom_id_to_idx[bondi];
	    rbmd::Id bondjj = atom_id_to_idx[bondj];
	  // printf("--------test--tid1:%i---bondi:%i--bondj:%i\n",
	    // tid1,bondi,bondj);

	   //printf("--------test--tid1:%i---idx_bondi:%i---idx_bondj:%i\n"
	    // ,tid1,atom_id_to_idx[bondi],atom_id_to_idx[bondj]);
	    rbmd::Id bondtype = bond_type[tid1];
	    rbmd::Real k = bond_coeffs_k[bondtype];
	    rbmd::Real equilibrium_bond = bond_coeffs_equilibrium[bondtype];
	    //printf("--------test---bondtype:%i---bond_coeffs_k:%f---equilibrium_bond:%f\n",
             // bondtype,k,equilibrium_bond);

	    rbmd::Real x12 = px[bondjj]-px[bondii];
	    rbmd::Real y12 = py[bondjj]-py[bondii];
	    rbmd::Real z12 = pz[bondjj]-pz[bondii];
	    //printf("--------test---x12:%f---y12:%f---z12:%f\n",x12,y12,z12);
	    MinImageDistance(box,x12,y12,z12);
	    rbmd::Real  dis_12 = SQRT(x12 * x12 + y12 * y12 + z12 * z12);
	    rbmd::Real dr = dis_12 - equilibrium_bond;
	    rbmd::Real rk = k * dr;

	    //energy
	    local_energy_bond += rk * dr;
	    //printf("--------test--dr:%f--rk:%f----local_energy_bond:%f\n"
	     // ,dr,rk,local_energy_bond);

	    rbmd::Real forcebondij;
	    if (dis_12 > 0.01)
	      forcebondij =  -2.0 * rk / dis_12;
	    else
	      forcebondij = 0.0;

	    // // apply force to each of 2 atoms
	    // fx[bondi] += forcebondij * x12;
	    // fy[bondi] += forcebondij * y12;
	    // fz[bondi] += forcebondij * z12;
	    //
	    // fx[bondj] -= forcebondij * x12;
	    // fy[bondj] -= forcebondij * y12;
	    // fz[bondj] -= forcebondij * z12;

	    // 计算作用力分量
	    rbmd::Real fx_ij = forcebondij * x12;
	    rbmd::Real fy_ij = forcebondij * y12;
	    rbmd::Real fz_ij = forcebondij * z12;
	    //printf("force_bond, %i  %f %f %f\n",tid1,fx_ij,fy_ij,fz_ij);

	    atomicAdd(&fx[bondii], fx_ij); // 将作用力添加到原子bondi
	    atomicAdd(&fy[bondii], fy_ij);
	    atomicAdd(&fz[bondii], fz_ij);

	    atomicAdd(&fx[bondjj], -fx_ij); // 对于bondj施加相反的作用力
	    atomicAdd(&fy[bondjj], -fy_ij);
	    atomicAdd(&fz[bondjj], -fz_ij);

	    // // 保存bondi和bondj的力到输出数组中
	    // atom_ids_out[2 * tid1]     = bondii;   // bondi -> bondii
	    // fx[2 * tid1]           = fx_ij;
	    // fy[2 * tid1]           = fy_ij;
	    // fz[2 * tid1]           = fz_ij;
	    //
	    // atom_ids_out[2 * tid1 + 1] = bondjj;   // bondj -> bondjj
	    // fx[2 * tid1 + 1]       = -fx_ij;   // 反向力
	    // fy[2 * tid1 + 1]       = -fy_ij;
	    // fz[2 * tid1 + 1]       = -fz_ij;

	  }

	  rbmd::Real block_sum =
    hipcub::BlockReduce<rbmd::Real, BLOCK_SIZE>(temp_storage).Sum(local_energy_bond);

	  if (threadIdx.x == 0){
	    atomicAdd(energy_bond, block_sum);
	  }
	}


      //
       __global__ void ComputeAngleForce(
       Box* box,
       const rbmd::Id num_anglels,
       const rbmd::Id* atom_id_to_idx,
       const rbmd::Real* anglel_coeffs_k,
       const rbmd::Real* anglel_coeffs_equilibrium,
       const rbmd::Id* anglel_type,
       const rbmd::Id* anglelisti,
       const rbmd::Id* anglelistj,
       const rbmd::Id* anglelistk,
       const rbmd::Real* px,
       const rbmd::Real* py,
       const rbmd::Real* pz,
       rbmd::Real* fx,
       rbmd::Real* fy,
       rbmd::Real* fz,
       rbmd::Real* energy_angle)
        {
          __shared__ typename hipcub::BlockReduce<rbmd::Real, BLOCK_SIZE>::TempStorage
temp_storage;
          rbmd::Real local_energy_angle  = 0;

	  unsigned int tid1 = blockIdx.x * blockDim.x + threadIdx.x;
	  if (tid1 < num_anglels)
	  {
	    rbmd::Id angleli = anglelisti[tid1];
	    rbmd::Id anglelj = anglelistj[tid1];
	    rbmd::Id anglelk = anglelistk[tid1];
	    rbmd::Id anglelii = atom_id_to_idx[angleli];
	    rbmd::Id angleljj = atom_id_to_idx[anglelj];
	    rbmd::Id anglelkk = atom_id_to_idx[anglelk];
	    //printf("--------test---angleli:%i---anglelj:%i---anglelk:%i\n",
	     // angleli,anglelj,anglelk);

	    rbmd::Id angleltype = anglel_type[tid1];
	    rbmd::Real k = anglel_coeffs_k[angleltype];
	    rbmd::Real equilibrium_angle = anglel_coeffs_equilibrium[angleltype];
	   // printf("--------test---angleltype:%i---k:%f---equilibrium_angle:%f\n",
	    // angleltype,k,equilibrium_angle);

	    rbmd::Real x12 = px[anglelii]-px[angleljj];  //i j
	    rbmd::Real y12 = py[anglelii]-py[angleljj];
	    rbmd::Real z12 = pz[anglelii]-pz[angleljj];
	    MinImageDistance(box,x12,y12,z12);

	    rbmd::Real x23 = px[anglelkk]-px[angleljj]; // k j
	    rbmd::Real y23 = py[anglelkk]-py[angleljj];
	    rbmd::Real z23 = pz[anglelkk]-pz[angleljj];
	    MinImageDistance(box,x23,y23,z23);

	    rbmd::Real dis_12_2 = x12 * x12 + y12 * y12 + z12 * z12;
	    rbmd::Real  dis_12 = SQRT(dis_12_2);

	    rbmd::Real dis_23_2 = x23 * x23 + y23 * y23 + z23 * z23;
	    rbmd::Real  dis_23 = SQRT(dis_23_2);

	    rbmd::Real cosangle = x12*x23 + y12*y23 + z12*z23;
	    cosangle /= dis_12 * dis_23;

	    if (cosangle > 1.0)
	      cosangle = 1.0;
	    if (cosangle < -1.0)
	      cosangle = -1.0;
	    rbmd::Real s = SQRT(1.0 - cosangle * cosangle);

	    const rbmd::Real SMALL = 0.001;
	    if (s < SMALL)
	      s = SMALL;
	    s = 1.0 / s;

	    rbmd::Real dtheta = ACOS(cosangle) - (equilibrium_angle * M_PI) / 180;
	    rbmd::Real tk = k * dtheta;

	    //energy
	    local_energy_angle += tk * dtheta;
	    //printf("--------test--dtheta:%f--tk:%f----local_energy_angle:%f\n"
	     // ,dtheta,tk ,local_energy_angle);

	    rbmd::Real a = -2.0 * tk * s;
	    rbmd::Real a11 = a * cosangle / dis_12_2;
	    rbmd::Real a12 = -a / (dis_12 * dis_23);
	    rbmd::Real a22 = a * cosangle / dis_23_2;

	    rbmd::Real force_anglei_x,force_anglei_y,force_anglei_z;
	    rbmd::Real force_anglek_x,force_anglek_y,force_anglek_z;
	    rbmd::Real force_anglej_x,force_anglej_y,force_anglej_z;

	    force_anglei_x = a11 * x12 + a12 * x23;
	    force_anglei_y = a11 * y12 + a12 * y23;
	    force_anglei_z = a11 * z12 + a12 * z23;

	    force_anglek_x = a22 * x23 + a12 * x12;
	    force_anglek_y = a22 * y23 + a12 * y12;
	    force_anglek_z = a22 * z23 + a12 * z12;

	    force_anglej_x = -(force_anglei_x+force_anglek_x);
	    force_anglej_y = -(force_anglei_y+force_anglek_y);
	    force_anglej_z = -(force_anglei_z+force_anglek_z);

// 	    printf("force_angle_i, %i  %f %f %f\n",tid1,
// 	    force_anglei_x,force_anglei_y,force_anglei_z);
//
// 	    printf("force_angle_j, %i  %f %f %f\n",tid1,
// force_anglej_x,force_anglej_y,force_anglej_z);
//
// 	    printf("force_angle_k, %i  %f %f %f\n",tid1,
// force_anglek_x,force_anglek_y,force_anglek_z);
	    atomicAdd(&fx[anglelii], force_anglei_x);
	    atomicAdd(&fy[anglelii], force_anglei_y);
	    atomicAdd(&fz[anglelii], force_anglei_z);

	    atomicAdd(&fx[anglelkk], force_anglek_x);
	    atomicAdd(&fy[anglelkk], force_anglek_y);
	    atomicAdd(&fz[anglelkk], force_anglek_z);

	    atomicAdd(&fx[angleljj], force_anglej_x);
	    atomicAdd(&fy[angleljj], force_anglej_y);
	    atomicAdd(&fz[angleljj], force_anglej_z);
	  }

          rbmd::Real block_sum =
hipcub::BlockReduce<rbmd::Real, BLOCK_SIZE>(temp_storage).Sum(local_energy_angle);

          if (threadIdx.x == 0){
            atomicAdd(energy_angle, block_sum);
          }
	}

       __global__ void ComputeDihedralForce(
         Box* box,
         const rbmd::Id num_dihedrals,
         const rbmd::Id num_atoms,
         const rbmd::Real* dihedral_coeffs_k,
         const rbmd::Id* dihedral_coeffs_sign ,
         const rbmd::Id* dihedral_coeffs_multiplicity ,
         const rbmd::Id* dihedral_type,
         const int4* dihedrallist,
         const rbmd::Real* px,
         const rbmd::Real* py,
         const rbmd::Real* pz,
         rbmd::Real* fx,
         rbmd::Real* fy,
         rbmd::Real* fz,
         rbmd::Real* energy_dihedral)
       {
          __shared__ typename hipcub::BlockReduce<rbmd::Real, BLOCK_SIZE>::TempStorage
temp_storage;
          rbmd::Real local_energy_dihedral  = 0;

          unsigned int tid1 = blockIdx.x * blockDim.x + threadIdx.x;
          if (tid1 < num_dihedrals)
          {

            rbmd::Id dihedrali = dihedrallist[tid1].x;
            rbmd::Id dihedralj = dihedrallist[tid1].y;
            rbmd::Id dihedralk = dihedrallist[tid1].z;
            rbmd::Id dihedralw = dihedrallist[tid1].w;
            rbmd::Real cos_shift, sin_shift;
            if (dihedral_coeffs_sign[tid1] == 1)
            {
              cos_shift = 1.0;
              sin_shift = 0.0;
            }
            else
            {
              cos_shift = -1.0;
              sin_shift = 0.0;
            }

            rbmd::Id bondtype = dihedral_type[tid1];
            rbmd::Real k = dihedral_coeffs_k[bondtype];

            rbmd::Real x12 = px[dihedrali]-px[dihedralj];  // i j =vb1
            rbmd::Real y12 = py[dihedrali]-py[dihedralj];
            rbmd::Real z12 = pz[dihedrali]-pz[dihedralj];
            MinImageDistance(box,x12,y12,z12);

            rbmd::Real x23 = px[dihedralk]-px[dihedralj];  //  k j=vb2
            rbmd::Real y23 = py[dihedralk]-py[dihedralj];
            rbmd::Real z23 = pz[dihedralk]-pz[dihedralj];
            MinImageDistance(box,x23,y23,z23);

            rbmd::Real x23m = -x23;  // =vb2m
            rbmd::Real y23m = -y23;
            rbmd::Real z23m = -z23;

            rbmd::Real x34 = px[dihedralw]-px[dihedralk];  // w k   =vb3
            rbmd::Real y34 = py[dihedralw]-py[dihedralk];
            rbmd::Real z34 = pz[dihedralw]-pz[dihedralk];
            MinImageDistance(box,x34,y34,z34);
            // c,s calculation

            rbmd::Real ax = y12 * z23m - z12 * y23m;
            rbmd::Real ay = z12 * x23m - x12 * z23m;
            rbmd::Real az = x12 * y23m - y12 * x23m;
            rbmd::Real bx = y34 * z23m - z34 * y23m;
            rbmd::Real by = z34 * x23m - x34 * z23m;
            rbmd::Real bz = x34 * y23m - z34 * x23m;

            rbmd::Real rasq = ax * ax + ay * ay + az * az;
            rbmd::Real rbsq = bx * bx + by * by + bz * bz;
            rbmd::Real rgsq = x23m * x23m + y23m * y23m + z23m * z23m;
            rbmd::Real rg = SQRT(rgsq);

            rbmd::Real  rginv, ra2inv, rb2inv;
            rginv = ra2inv = rb2inv = 0.0;
            if (rg > 0) rginv = 1.0 / rg;

            if (rasq > 0) ra2inv = 1.0 / rasq;

            if (rbsq > 0)  rb2inv = 1.0 / rbsq;

            rbmd::Real rabinv = SQRT(ra2inv * rb2inv);

            rbmd::Real c = (ax * bx + ay * by + az * bz) * rabinv;
            rbmd::Real s = rg * rabinv * (ax * x34 + ay * y34 + az * z34);

            if (c > 1.0) c = 1.0;
            if (c < -1.0)  c = -1.0;

            rbmd::Id m = dihedral_coeffs_multiplicity[tid1];
            rbmd::Real p = 1.0;
            rbmd::Real ddf1, df1;
            ddf1 = df1 = 0.0;
            for (rbmd::Id i = 0; i < m; i++)
            {
              ddf1 = p * c - df1 * s;
              df1 = p * s + df1 * c;
              p = ddf1;
            }

            p = p * cos_shift + df1 * sin_shift;
            df1 = df1 * cos_shift - ddf1 * sin_shift;
            df1 *= -m;
            p += 1.0;

            if (m == 0)
            {
              p = 1.0 + cos_shift;
              //p = 1.0 + cos_shift[type];
              df1 = 0.0;
            }

            //energy
            local_energy_dihedral += k * p;

            rbmd::Real fg = x12 * x23m + y12 * y23m + z12 * z23m;
            rbmd::Real hg = x34 * x23m + y34 * y23m + z34 * z23m;

            rbmd::Real  fga = fg * ra2inv * rginv;
            rbmd::Real  hgb = hg * rb2inv * rginv;
            rbmd::Real  gaa = -ra2inv * rg;
            rbmd::Real  gbb = rb2inv * rg;

            rbmd::Real  dtfx = gaa * ax;
            rbmd::Real  dtfy = gaa * ay;
            rbmd::Real  dtfz = gaa * az;
            rbmd::Real  dtgx = fga * ax - hgb * bx;
            rbmd::Real  dtgy = fga * ay - hgb * by;
            rbmd::Real  dtgz = fga * az - hgb * bz;
            rbmd::Real  dthx = gbb * bx;
            rbmd::Real  dthy = gbb * by;
            rbmd::Real  dthz = gbb * bz;

            rbmd::Real df = -k  * df1;
            //df = -k[type] * df1;
            rbmd::Real sx2 = df * dtgx;
            rbmd::Real sy2 = df * dtgy;
            rbmd::Real sz2 = df * dtgz;

            //force
            rbmd::Real force_dihedrali_x,force_dihedrali_y,force_dihedrali_z;
            rbmd::Real force_dihedralj_x,force_dihedralj_y,force_dihedralj_z;
            rbmd::Real force_dihedralk_x,force_dihedralk_y,force_dihedralk_z;
            rbmd::Real force_dihedralw_x,force_dihedralw_y,force_dihedralw_z;

            force_dihedrali_x = df * dtfx;
            force_dihedrali_y = df * dtfy;
            force_dihedrali_z = df * dtfz;

            force_dihedralj_x = sx2 - force_dihedrali_x;
            force_dihedralj_y = sy2 - force_dihedrali_y;
            force_dihedralj_z = sz2 - force_dihedrali_z;

            force_dihedralw_x = df * dthx;
            force_dihedralw_y = df * dthy;
            force_dihedralw_z = df * dthz;

            force_dihedralk_x = -sx2 - force_dihedralw_x;
            force_dihedralk_y = -sy2 - force_dihedralw_y;
            force_dihedralk_z = -sz2 - force_dihedralw_z;

            atomicAdd(&fx[dihedrali], force_dihedrali_x);
            atomicAdd(&fy[dihedrali], force_dihedrali_y);
            atomicAdd(&fz[dihedrali], force_dihedrali_z);

            atomicAdd(&fx[dihedralj], force_dihedralj_x);
            atomicAdd(&fy[dihedralj], force_dihedralj_y);
            atomicAdd(&fz[dihedralj], force_dihedralj_z);

            atomicAdd(&fx[dihedralw], force_dihedralw_x);
            atomicAdd(&fy[dihedralw], force_dihedralw_y);
            atomicAdd(&fz[dihedralw], force_dihedralw_z);

            atomicAdd(&fx[dihedralk], force_dihedralk_x);
            atomicAdd(&fy[dihedralk], force_dihedralk_y);
            atomicAdd(&fz[dihedralk], force_dihedralk_z);
          }

          rbmd::Real block_sum =
hipcub::BlockReduce<rbmd::Real, BLOCK_SIZE>(temp_storage).Sum(local_energy_dihedral);

          if (threadIdx.x == 0){
            atomicAdd(energy_dihedral, block_sum);
          }
       }


        /////////////////////
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
			fx, fy, fz, total_evdwl));
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
			fx, fy, fz, virial_xx, virial_yy, virial_zz, virial_xy, virial_xz, virial_yz, total_evdwl));
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
			(box, rs, rc, num_atoms, neighbor_sample_num, pice_num,
				atoms_type, molecular_type, sigma, eps,
				start_id, end_id, id_verletlist, id_random_neighbor, random_neighbor_num,
				px, py, pz, fx, fy, fz));
	}

	//RBL: Fix LJForce
	void FixRBLForceOp<device::DEVICE_GPU>::operator()(
		const rbmd::Id num_atoms,
		const rbmd::Real corr_value_x,
		const rbmd::Real corr_value_y,
		const rbmd::Real corr_value_z,
		rbmd::Real* fx,
		rbmd::Real* fy,
		rbmd::Real* fz)
	{
		unsigned int blocks_per_grid = (num_atoms + BLOCK_SIZE - 1) / BLOCK_SIZE;

		CHECK_KERNEL(FixRBLForce <<<blocks_per_grid, BLOCK_SIZE, 0, 0 >>>
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

	void LJCutCoulForceOp<device::DEVICE_GPU>::operator()(
		Box* box,ERFTable* erf_table,
		const rbmd::Real cut_off,
		const rbmd::Id num_atoms,
		const rbmd::Real alpha,
		const rbmd::Real qqr2e,
		const rbmd::Id* atoms_type,
		const rbmd::Id* molecular_id,
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

		CHECK_KERNEL(ComputeLJCutCoulForce <<<blocks_per_grid, BLOCK_SIZE, 0, 0 >>> (box,erf_table, cut_off, num_atoms, alpha, qqr2e, atoms_type, molecular_id,
			sigma, eps, start_id, end_id, id_verletlist, charge, px, py, pz,
			fx, fy, fz, total_evdwl, total_ecoul));
	}

        void LJCutCoulRBLForceOp<device::DEVICE_GPU>::operator()(
                Box* box,ERFTable* erf_table,
                const rbmd::Real rs,
                const rbmd::Real rc,
                const rbmd::Id num_atoms,
                const rbmd::Id neighbor_sample_num,
                const rbmd::Id pice_num,
                const rbmd::Real alpha,
                const rbmd::Real qqr2e,
                const rbmd::Id* atoms_type,
                const rbmd::Id* molecular_type,
                const rbmd::Real* sigma,
                const rbmd::Real* eps,
                const rbmd::Id* start_id,
                const rbmd::Id* end_id,
                const rbmd::Id* id_verletlist,
                const rbmd::Id* id_random_neighbor,
                const rbmd::Id* random_neighbor_num,
                const rbmd::Real* charge,
                const rbmd::Real* px,
                const rbmd::Real* py,
                const rbmd::Real* pz,
                rbmd::Real* fx,
                rbmd::Real* fy,
                rbmd::Real* fz)
        {
	  unsigned int blocks_per_grid = (num_atoms + BLOCK_SIZE - 1) / BLOCK_SIZE;

	  CHECK_KERNEL(ComputeLJCutCoulRBLForce <<<blocks_per_grid, BLOCK_SIZE, 0, 0 >>>
	    (box, erf_table,rs,rc, num_atoms,neighbor_sample_num,pice_num,alpha, qqr2e,
	      atoms_type, molecular_type,sigma, eps, start_id, end_id, id_verletlist,
	      id_random_neighbor, random_neighbor_num,charge, px, py, pz,fx,fy,fz));
	}

       void LJCutCoulEnergyOp<device::DEVICE_GPU>::operator()(
        Box* box,ERFTable* erf_table,
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
        rbmd::Real* total_evdwl,
        rbmd::Real* total_ecoul)
	{
	  unsigned int blocks_per_grid = (num_atoms + BLOCK_SIZE - 1) / BLOCK_SIZE;

	  CHECK_KERNEL(ComputeLJCutCoulEnergy <<<blocks_per_grid, BLOCK_SIZE, 0, 0 >>>
	    (box, erf_table,cut_off, num_atoms, alpha, qqr2e, atoms_type, molecular_type,
	    sigma, eps, start_id, end_id, id_verletlist, charge, px, py, pz,
	    total_evdwl, total_ecoul));
	}

	//Charge Structure Factor
	void ComputeChargeStructureFactorOp<device::DEVICE_GPU>::operator()(
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

		CHECK_KERNEL(ComputeChargeStructureFactor <<<blocks_per_grid, BLOCK_SIZE, 0, 0 >>>
			(num_atoms, K, charge, px, py, pz, density_real, density_imag));

	}

	//EwaldForce
	void ComputeEwaldForceOp<device::DEVICE_GPU>::operator()(
		Box* box,
		const rbmd::Id num_atoms,
		const rbmd::Id  Kmax,
		const rbmd::Real alpha,
		const rbmd::Real qqr2e,
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
			(box, num_atoms, Kmax, alpha, qqr2e, real_array, imag_array, charge,
				px, py, pz, fx, fy, fz));

	}
	//sq_charge
	void SqchargeOp<device::DEVICE_GPU>::operator()(
		const rbmd::Id num_atoms,
		const rbmd::Real* charge,
		rbmd::Real* sq_charge)
	{
		unsigned int blocks_per_grid = (num_atoms + BLOCK_SIZE - 1) / BLOCK_SIZE;

		CHECK_KERNEL(ComputeSqCharge <<<blocks_per_grid, BLOCK_SIZE, 0, 0 >>>
			(num_atoms, charge, sq_charge));
	}

        //
        void GenerateIndexArrayOp<device::DEVICE_GPU>::operator()(
                const rbmd::Id num_atoms,
                const rbmd::Id RBE_P,
                rbmd::Id* psample_key)
	{
	  unsigned int blocks_per_grid = ((num_atoms* RBE_P)+ BLOCK_SIZE - 1) / BLOCK_SIZE;

	  CHECK_KERNEL(GenerateIndexArray <<<blocks_per_grid, BLOCK_SIZE, 0, 0 >>>
                  (num_atoms, RBE_P,psample_key));
	}

     //RBE
       void ComputePnumberChargeStructureFactorOp<device::DEVICE_GPU>::operator()(
                Box* box,
                const rbmd::Id num_atoms,
                const rbmd::Id p_number,
                const rbmd::Real* charge,
                const rbmd::Real* p_sample_x,
                const rbmd::Real* p_sample_y,
                const rbmd::Real* p_sample_z,
                const rbmd::Real* px,
                const rbmd::Real* py,
                const rbmd::Real* pz,
                rbmd::Real* density_real,
                rbmd::Real* density_imag)
	{
	  unsigned int blocks_per_grid = (num_atoms + BLOCK_SIZE - 1) / BLOCK_SIZE;

	  CHECK_KERNEL(ComputePnumberChargeStructureFactor <<<blocks_per_grid, BLOCK_SIZE, 0, 0 >>>
                  (box,num_atoms, p_number, charge, p_sample_x,p_sample_y,p_sample_z,
                    px, py, pz, density_real, density_imag));

	}


        void ComputeRBEForceOp<device::DEVICE_GPU>::operator()(
           Box* box,
           const rbmd::Id num_atoms,
           const rbmd::Id  p_number,
           const rbmd::Real alpha,
           const rbmd::Real qqr2e,
           const rbmd::Real* real_array,
           const rbmd::Real* imag_array,
           const rbmd::Real* charge,
           const rbmd::Real* p_sample_x,
           const rbmd::Real* p_sample_y,
           const rbmd::Real* p_sample_z,
           const rbmd::Real* px,
           const rbmd::Real* py,
           const rbmd::Real* pz,
           rbmd::Real* fx,
           rbmd::Real* fy,
           rbmd::Real* fz)
	{
	  unsigned int blocks_per_grid = (num_atoms + BLOCK_SIZE - 1) / BLOCK_SIZE;

	  CHECK_KERNEL(ComputeRBEForce <<<blocks_per_grid, BLOCK_SIZE, 0, 0 >>>
                  (box, num_atoms, p_number, alpha,qqr2e, real_array, imag_array, charge,
                         p_sample_x,p_sample_y,p_sample_z, px, py, pz, fx, fy, fz));

	}

       void ComputeSpecialCoulForceOp<device::DEVICE_GPU>::operator()(
       Box* box,
       const rbmd::Id num_atoms,
       const rbmd::Real qqr2e,
       const rbmd::Id* atoms_id,
       const rbmd::Id*  atoms_vec,
       const rbmd::Id*  atoms_offset,
       const rbmd::Id*  atom_count,
       const rbmd::Id*  special_ids,
       const rbmd::Real*  special_weights,
       const rbmd::Id*  special_offset,
       const rbmd::Id*  special_count,
       const rbmd::Real* charge,
       const rbmd::Real* px,
       const rbmd::Real* py,
       const rbmd::Real* pz,
       rbmd::Real* fx,
       rbmd::Real* fy,
       rbmd::Real* fz,
       rbmd::Real* total_especial_coul)
      {
          unsigned int blocks_per_grid = (num_atoms + BLOCK_SIZE - 1) / BLOCK_SIZE;

          CHECK_KERNEL(ComputeSpecialCoulForce <<<blocks_per_grid, BLOCK_SIZE, 0, 0 >>>
                  (box,num_atoms,qqr2e,atoms_id,atoms_vec,atoms_offset,atom_count,
                    special_ids,special_weights,special_offset,special_count,
                    charge,px, py, pz, fx, fy, fz,total_especial_coul));
      }

//
      void ComputeBondForceOp<device::DEVICE_GPU>::operator()(
       Box* box,
       const rbmd::Id num_bonds,
       const rbmd::Id* atom_id_to_idx,
       const rbmd::Real* bond_coeffs_k,
       const rbmd::Real* bond_coeffs_equilibrium,
       const rbmd::Id* bond_type,
       const rbmd::Id* bondlisti,
       const rbmd::Id* bondlistj,
       const rbmd::Real* px,
       const rbmd::Real* py,
       const rbmd::Real* pz,
       rbmd::Real* fx,
       rbmd::Real* fy,
       rbmd::Real* fz,
       rbmd::Id* atom_ids_out,
       rbmd::Real* energy_bond)
       {
	  unsigned int blocks_per_grid = (num_bonds + BLOCK_SIZE - 1) / BLOCK_SIZE;

	  CHECK_KERNEL(ComputeBondForce <<<blocks_per_grid, BLOCK_SIZE, 0, 0 >>>
                  (box, num_bonds,atom_id_to_idx,bond_coeffs_k,bond_coeffs_equilibrium,bond_type,
                    bondlisti,bondlistj,px, py, pz, fx, fy, fz,atom_ids_out,energy_bond));
       }

//
       void ComputeAngleForceOp<device::DEVICE_GPU>::operator()(
        Box* box,
        const rbmd::Id num_anglels,
        const rbmd::Id* _atom_id_to_idx,
        const rbmd::Real* anglel_coeffs_k,
        const rbmd::Real* anglel_coeffs_equilibrium,
        const rbmd::Id* anglel_type,
        const rbmd::Id* anglelisti,
        const rbmd::Id* anglelistj,
        const rbmd::Id* anglelistk,
        const rbmd::Real* px,
        const rbmd::Real* py,
        const rbmd::Real* pz,
        rbmd::Real* fx,
        rbmd::Real* fy,
        rbmd::Real* fz,
        rbmd::Real* energy_angle)
        {
	  unsigned int blocks_per_grid = (num_anglels + BLOCK_SIZE - 1) / BLOCK_SIZE;

	  CHECK_KERNEL(ComputeAngleForce <<<blocks_per_grid, BLOCK_SIZE, 0, 0 >>>
                  (box, num_anglels,_atom_id_to_idx,anglel_coeffs_k,
                    anglel_coeffs_equilibrium,anglel_type,
                    anglelisti,anglelistj,anglelistk,px, py, pz, fx, fy, fz,energy_angle));
	}

//

        void ComputeDihedralForceOp<device::DEVICE_GPU>::operator()(
          Box* box,
          const rbmd::Id num_dihedrals,
          const rbmd::Id num_atoms,
          const rbmd::Real* dihedral_coeffs_k,
          const rbmd::Id* dihedral_coeffs_sign ,
          const rbmd::Id* dihedral_coeffs_multiplicity ,
          const rbmd::Id* dihedral_type,
          const int4* dihedrallist,
          const rbmd::Real* px,
          const rbmd::Real* py,
          const rbmd::Real* pz,
          rbmd::Real* fx,
          rbmd::Real* fy,
          rbmd::Real* fz,
          rbmd::Real* energy_dihedral)
       {
          unsigned int blocks_per_grid = (num_dihedrals + BLOCK_SIZE - 1) / BLOCK_SIZE;

          CHECK_KERNEL(ComputeDihedralForce <<<blocks_per_grid, BLOCK_SIZE, 0, 0 >>>
                  (box, num_dihedrals,num_atoms,dihedral_coeffs_k,
                    dihedral_coeffs_sign,dihedral_coeffs_multiplicity,dihedral_type,
                    dihedrallist,px, py, pz, fx, fy, fz,energy_dihedral));

       }

}

