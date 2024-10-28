#pragma once


#include <random>
#include <numeric>
#include "types.h"
#include "rbmd_define.h"
#include "model/box.h"

template<typename T>
T RandomValue(const rbmd::Real& Min, const rbmd::Real& Max){

  std::mt19937 gen(std::random_device{}());
  std::uniform_real_distribution<T> dis(Min, Max);

  return dis(gen);
};


//
struct RBEPSAMPLE
{
  rbmd::Real _alpha;
  Box* _box;
  rbmd::Id _P;
  bool _RBE_random;

  rbmd::Real Compute_S() const
  {
    const Real3& factor = Compute_H();
    rbmd::Real factor_3 = factor.data[0] * factor.data[1] * factor.data[2];
    rbmd::Real S = factor_3 - 1;
    return S;
  }

  Real3 Compute_H() const
  {
    Real3 H;
    for (rbmd::Id i = 0; i < 3; ++i)
    {
      const rbmd::Real factor = -(_alpha * _box->_length[i] * _box->_length[i]);
      H.data[i] = 0.0;

      for (rbmd::Id m = -10; m <= 10; m++)
      {
        rbmd::Real expx = m * m * factor;
        H.data[i] += EXP(expx);
      }
      H.data[i] *= SQRT(-(factor) / M_PI);
    }

    return H;
  }

  rbmd::Real MH_Algorithm(
    rbmd::Real m,
    rbmd::Real mu,
    const Real3& sigma,
    rbmd::Id dimension) const
  {
    rbmd::Real x_wait = FetchSample_1D(mu, sigma.data[dimension]);
    rbmd::Real m_wait = rbmd::Real(ROUND(x_wait));
    rbmd::Real Prob = (Distribution_P(m_wait, dimension) / Distribution_P(m, dimension)) *
      (Distribution_q(m, dimension) / Distribution_q(m_wait, dimension));
    Prob = MIN(Prob, rbmd::Real(1.0));

    if (_RBE_random)
    {
      rbmd::Real u = RandomValue<rbmd::Real>(0.0, 1.0); //random?
      if (u <= Prob)
        m = m_wait;
      return m;
    }
    else
    {
      rbmd::Real u = 0.5;
      if (u <= Prob)
        m = m_wait;
      return m;
    }
  }

  rbmd::Real Distribution_P(const rbmd::Real& x, const rbmd::Id dimension) const
  {
    rbmd::Real P_m = EXP(-POW(2 * M_PI * x / _box->_length[dimension], 2)
      / (4 * _alpha));
    Real3 H = Compute_H();
    P_m = P_m / H.data[dimension];
    return P_m;
  }

  rbmd::Real Distribution_q(const rbmd::Real& x, const rbmd::Id dimension) const
  {
    rbmd::Real q_m;
    if (x == 0)
    {
      q_m = ERF((1.0 / 2) /(SQRT(_alpha * POW(_box->_length[dimension], 2) /
        POW(M_PI, 2))));
    }
    else
      q_m = (ERF(((1.0 / 2) + ABS(x)) /
        (SQRT(_alpha * POW(_box->_length[dimension], 2)
        / POW(M_PI, 2)))) -ERF((ABS(x) - (1.0 / 2)) /
               (SQRT(_alpha * POW(_box->_length[dimension], 2) /
                 POW(M_PI, 2))))) /2;
    return q_m;
  }

   rbmd::Real FetchSample_1D(const rbmd::Real& mu,
                             const rbmd::Real& sigma) const // Fetch 1D sample from Gaussion contribution
  {
    rbmd::Real U1, U2, epsilon;
    epsilon = 1e-6;
    if (_RBE_random)
    {
      do
      {
        U1 = RandomValue< rbmd::Real>(0.0, 1.0);
      } while (U1 < epsilon);
      U2 = RandomValue< rbmd::Real>(0.0, 1.0);
      Real2 ChooseSample{ 0.0, 0.0 };
      ChooseSample.data[0] = sigma * SQRT(-2.0 * LOG(U1)) * COS(2 * M_PI * U2) + mu;
      return ChooseSample.data[0];
    }
    else
    {
      U1 = 0.5;
      U2 = 0.5;
      Real2 ChooseSample{ 0.0, 0.0 };
      ChooseSample.data[0] = sigma * SQRT(-2.0 * LOG(U1)) * COS(2 * M_PI * U2) + mu;
      return ChooseSample.data[0];
    }
  }

  void Fetch_P_Sample(
    const rbmd::Real& mu,
    const Real3& sigma,
    rbmd::Real* P_Sample_x,
    rbmd::Real* P_Sample_y,
    rbmd::Real* P_Sample_z) const
  {
    rbmd::Real epsilonx = 1e-6; // precision
    Real3 X_0;
    do
    {
      X_0 = { rbmd::Real(ROUND(FetchSample_1D(mu, sigma.x))),
              rbmd::Real(ROUND(FetchSample_1D(mu, sigma.y))),
              rbmd::Real(ROUND(FetchSample_1D(mu, sigma.z)))};
    } while (ABS(X_0.x) < epsilonx && ABS(X_0.y) < epsilonx &&
             ABS(X_0.z) < epsilonx);
    /// 记录第一个样本
    P_Sample_x[0] = X_0.x;
    P_Sample_y[0] = X_0.y;
    P_Sample_z[0] = X_0.z;

    //
    for (rbmd::Id i = 1; i < _P; i++)
    {
      Real3 X_1 = { MH_Algorithm(X_0.x, mu, sigma, 0),
                    MH_Algorithm(X_0.y, mu, sigma, 1),
                    MH_Algorithm(X_0.z, mu, sigma, 2) };
      P_Sample_x[i]= X_1.x;
      P_Sample_y[i]= X_1.y;
      P_Sample_z[i]= X_1.z;

      X_0 = X_1;
      if (ABS(X_1.x) < epsilonx && ABS(X_1.y) < epsilonx &&
          ABS(X_1.z) < epsilonx)
      {
        i = i - 1;  //// 维持样本数量
        continue;
      }
    }
  }
};