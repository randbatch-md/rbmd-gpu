#pragma once

#include <cmath>
#include <limits>
#include <algorithm> // For std::swap and std::fabs
#include "../common/rbmd_define.h"

namespace MathLib {

/**
 * @brief Calculate the dot product of two 3D vectors.
 */
inline rbmd::Real dot3(const rbmd::Real* a, const rbmd::Real* b) {
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

/**
 * @brief Calculate the cross product of two 3D vectors c = a x b.
 */
inline void cross3(const rbmd::Real* a, const rbmd::Real* b, rbmd::Real* c) {
    c[0] = a[1] * b[2] - a[2] * b[1];
    c[1] = a[2] * b[0] - a[0] * b[2];
    c[2] = a[0] * b[1] - a[1] * b[0];
}

/**
 * @brief Reverse the 3D vector: v = -v
 */
inline void negate3(rbmd::Real* v) {
    v[0] = -v[0];
    v[1] = -v[1];
    v[2] = -v[2];
}


/**
 * @brief 使用雅可比旋转法计算3x3对称矩阵的特征值和特征向量.
 * 此函数直接实现了算法，不依赖外部库.
 * @param a 输入的3x3对称矩阵 (内容会被修改)
 * @param d 输出的特征值数组 (大小为3)
 * @param v 输出的特征向量矩阵 (3x3), 特征向量是列向量
 * @return 0表示成功收敛, 1表示达到最大迭代次数仍未收敛
 */
inline int jacobi3(rbmd::Real a[3][3], rbmd::Real d[3], rbmd::Real v[3][3]) {
    constexpr int max_rotations = 50;

    // 初始化特征向量矩阵 v 为单位矩阵
    v[0][0] = v[1][1] = v[2][2] = 1.0;
    v[0][1] = v[0][2] = v[1][0] = v[1][2] = v[2][0] = v[2][1] = 0.0;

    // 初始化特征值 d 和中间向量 b
    d[0] = a[0][0];
    d[1] = a[1][1];
    d[2] = a[2][2];
    rbmd::Real b[3] = {d[0], d[1], d[2]};
    rbmd::Real z[3] = {0.0, 0.0, 0.0};

    for (int iter = 0; iter < max_rotations; ++iter) {
        // 检查是否收敛: 如果所有非对角线元素之和足够小
        rbmd::Real sum_off_diag = std::fabs(a[0][1]) + std::fabs(a[0][2]) + std::fabs(a[1][2]);
        if (sum_off_diag < 1.0e-12) {
            return 0; // 成功收敛
        }

        // 选择绝对值最大的非对角线元素进行旋转
        int p = 0, q = 1;
        if (std::fabs(a[0][2]) > std::fabs(a[p][q])) { p = 0; q = 2; }
        if (std::fabs(a[1][2]) > std::fabs(a[p][q])) { p = 1; q = 2; }

        rbmd::Real t, tau, c, s;
        if (std::fabs(a[p][q]) < 1.0e-12) {
            t = 0.0;
        } else {
            tau = (d[q] - d[p]) / (2.0 * a[p][q]);
            t = 1.0 / (std::fabs(tau) + std::sqrt(1.0 + tau * tau));
            if (tau < 0.0) t = -t;
        }

        c = 1.0 / std::sqrt(1.0 + t * t);
        s = t * c;

        // 更新矩阵 a, 特征值 d, 和特征向量 v
        rbmd::Real h = t * a[p][q];
        z[p] -= h;
        z[q] += h;
        d[p] -= h;
        d[q] += h;
        a[p][q] = 0.0;

        for (int j = 0; j < p; ++j) {
            rbmd::Real gj = a[j][p];
            rbmd::Real hj = a[j][q];
            a[j][p] = c * gj - s * hj;
            a[j][q] = s * gj + c * hj;
        }
        for (int j = p + 1; j < q; ++j) {
            rbmd::Real gpj = a[p][j];
            rbmd::Real hj = a[j][q];
            a[p][j] = c * gpj - s * hj;
            a[j][q] = s * gpj + c * hj;
        }
        for (int j = q + 1; j < 3; ++j) {
            rbmd::Real gpj = a[p][j];
            rbmd::Real gqj = a[q][j];
            a[p][j] = c * gpj - s * gqj;
            a[q][j] = s * gpj + c * gqj;
        }
        for (int j = 0; j < 3; ++j) {
            rbmd::Real vjp = v[j][p];
            rbmd::Real vjq = v[j][q];
            v[j][p] = c * vjp - s * vjq;
            v[j][q] = s * vjp + c * vjq;
        }
    }
    return 1; // 达到最大迭代次数，未收敛
}


/**
 * @brief 根据角动量和对角化的转动惯量张量计算角速度.
 * 忠实复现 LAMMPS 中 math_extra.cpp::angmom_to_omega 的逻辑.
 * @param angmom 角动量 L (space frame)
 * @param ex 第一个主轴 (特征向量)
 * @param ey 第二个主轴 (特征向量)
 * @param ez 第三个主轴 (特征向量)
 * @param idiag 主转动惯量 (特征值)
 * @param w 输出的角速度 omega (space frame)
 */
inline void angmom_to_omega(const rbmd::Real* angmom,
                            const rbmd::Real* ex, const rbmd::Real* ey, const rbmd::Real* ez,
                            const rbmd::Real* idiag, rbmd::Real* w) {
    rbmd::Real wbody[3];

    // 将空间坐标系下的角动量投影到物体坐标系（主轴坐标系）
    // w_body = I_body^-1 * L_body = (1/I_diag) * (P^T * L_space)
    if (std::fabs(idiag[0]) < 1.0e-12) {
        wbody[0] = 0.0;
    } else {
        wbody[0] = dot3(angmom, ex) / idiag[0];
    }

    if (std::fabs(idiag[1]) < 1.0e-12) {
        wbody[1] = 0.0;
    } else {
        wbody[1] = dot3(angmom, ey) / idiag[1];
    }

    if (std::fabs(idiag[2]) < 1.0e-12) {
        wbody[2] = 0.0;
    } else {
        wbody[2] = dot3(angmom, ez) / idiag[2];
    }

    // 将物体坐标系下的角速度转换回空间坐标系
    // w_space = P * w_body
    w[0] = wbody[0] * ex[0] + wbody[1] * ey[0] + wbody[2] * ez[0];
    w[1] = wbody[0] * ex[1] + wbody[1] * ey[1] + wbody[2] * ez[1];
    w[2] = wbody[0] * ex[2] + wbody[1] * ey[2] + wbody[2] * ez[2];
}

// /**
//  * @brief Jacobi :: General Jacobian Eigenvalue Solver
//  * * uses the Jacobi rotation method to  calculate all the eigenvalues and eigenvectors of
//  * an N-dimensional real symmetric matrix.
//  * * @tparam T  ( float, double)
//  * @tparam N The dimension of the matrix
//  */
// template <typename T, int N>
// class Jacobi {
// public:
//     /**
//      * @param max_rotations The maximum number of iterations for Jacobi method to prevent infinite loop
//      */
//     Jacobi(int max_rotations = 50) : m_max_rotations(max_rotations) {}
//
//     /**
//      * @brief Diagonalize the input symmetric matrix.
//      * * @param a 输入的 N x N 对称矩阵 (注意: 此矩阵内容在计算中会被修改)
//      * @param d 输出的特征值数组 (大小为 N)
//      * @param v 输出的特征向量矩阵 (N x N), 每个特征向量是矩阵的一列
//      * @return 0 表示成功收敛, 1 表示达到最大迭代次数仍未收敛
//      */
//     int Diagonalize(T a[N][N], T d[N], T v[N][N]) {
//         // --- 1. 初始化 ---
//         // 初始化特征向量矩阵 v 为单位矩阵
//         for (int i = 0; i < N; ++i) {
//             for (int j = 0; j < N; ++j) {
//                 v[i][j] = (i == j) ? 1.0 : 0.0;
//             }
//         }
//
//         // 初始化特征值 d 和中间向量 b
//         T b[N], z[N];
//         for (int i = 0; i < N; ++i) {
//             b[i] = d[i] = a[i][i];
//             z[i] = 0.0;
//         }
//
//         // --- 2. 主迭代循环 ---
//         for (int iter = 0; iter < m_max_rotations; ++iter) {
//             // 计算所有非对角线元素的绝对值之和，用于检查收敛
//             T sum_off_diag = 0.0;
//             for (int i = 0; i < N; ++i) {
//                 for (int j = i + 1; j < N; ++j) {
//                     sum_off_diag += std::fabs(a[i][j]);
//                 }
//             }
//
//             // 如果已经收敛，则成功返回
//             if (sum_off_diag < 1.0e-12) {
//                 return 0;
//             }
//
//             // --- 3. 选择旋转目标并计算旋转参数 ---
//             // 阈值，用于在早期迭代中跳过较小的非对角线元素
//             T threshold = (iter < 4) ? 0.2 * sum_off_diag / (N * N) : 0.0;
//
//             for (int p = 0; p < N; ++p) {
//                 for (int q = p + 1; q < N; ++q) {
//                     T g = 100.0 * std::fabs(a[p][q]);
//
//                     // 如果元素太小，跳过
//                     if (iter > 4 && (std::fabs(d[p]) + g == std::fabs(d[p]))
//                                  && (std::fabs(d[q]) + g == std::fabs(d[q]))) {
//                         a[p][q] = 0.0;
//                     } else if (std::fabs(a[p][q]) > threshold) {
//                         T h = d[q] - d[p];
//                         T t;
//                         if (std::fabs(h) + g == std::fabs(h)) {
//                             t = a[p][q] / h;
//                         } else {
//                             T theta = 0.5 * h / a[p][q];
//                             t = 1.0 / (std::fabs(theta) + std::sqrt(1.0 + theta * theta));
//                             if (theta < 0.0) t = -t;
//                         }
//
//                         T c = 1.0 / std::sqrt(1.0 + t * t);
//                         T s = t * c;
//                         T tau = s / (1.0 + c);
//                         h = t * a[p][q];
//
//                         // --- 4. 应用旋转 ---
//                         z[p] -= h;
//                         z[q] += h;
//                         d[p] -= h;
//                         d[q] += h;
//                         a[p][q] = 0.0;
//
//                         // 更新矩阵 a
//                         for (int j = 0; j < p; ++j)   { T gj = a[j][p]; T hj = a[j][q]; a[j][p] = gj - s * (hj + gj * tau); a[j][q] = hj + s * (gj - hj * tau); }
//                         for (int j = p + 1; j < q; ++j) { T gpj = a[p][j]; T hj = a[j][q]; a[p][j] = gpj - s * (hj + gpj * tau); a[j][q] = hj + s * (gpj - hj * tau); }
//                         for (int j = q + 1; j < N; ++j) { T gpj = a[p][j]; T gqj = a[q][j]; a[p][j] = gpj - s * (gqj + gpj * tau); a[q][j] = gqj + s * (gpj - gqj * tau); }
//
//                         // 更新特征向量矩阵 v
//                         for (int j = 0; j < N; ++j) { T vjp = v[j][p]; T vjq = v[j][q]; v[j][p] = vjp - s * (vjq + vjp * tau); v[j][q] = vjq + s * (vjp - vjq * tau); }
//                     }
//                 }
//             }
//
//             // 更新特征值 d
//             for (int i = 0; i < N; ++i) {
//                 b[i] += z[i];
//                 d[i] = b[i];
//                 z[i] = 0.0;
//             }
//         }
//         return 1; // 达到最大迭代次数，未收敛
//     }
//
// private:
//     int m_max_rotations;
// };
//
// /**
//  * @brief 便捷的调用接口，专门用于求解3x3对称矩阵的特征值问题.
//  * * 这是一个封装器 (Wrapper)，它实例化并调用通用的 Jacobi<T, 3> 类来完成工作。
//  * * @param a_const 输入的3x3对称矩阵 (const)
//  * @param d 输出的特征值数组 (大小为3)
//  * @param v 输出的特征向量矩阵 (3x3), 特征向量是列向量
//  * @return 0表示成功收敛, 1表示失败
//  */
// inline int jacobi3(const rbmd::Real a_const[3][3], rbmd::Real d[3], rbmd::Real v[3][3]) {
//   // 因为通用的 Diagonalize 方法会修改输入矩阵，所以我们先创建一个副本
//   rbmd::Real a_copy[3][3];
//   for(int i = 0; i < 3; ++i) {
//     for(int j = 0; j < 3; ++j) {
//       a_copy[i][j] = a_const[i][j];
//     }
//   }
//
//   // 实例化一个3维的 Jacobi 求解器
//   Jacobi<rbmd::Real, 3> eigensolver;
//
//   // 调用通用的对角化方法
//   int error_code = eigensolver.Diagonalize(a_copy, d, v);
//
//   // （可选，但与LAMMPS保持一致）: LAMMPS的约定是特征向量作为列向量，
//   // 我们的实现也是如此。LAMMPS的旧版jacobi3接口在最后做了一次转置，
//   // 可能是为了匹配某些代码的行向量约定。为保持清晰，我们这里不转置，
//   // 并明确注释特征向量是列向量。
//
//   // （可选）如果需要，可以对特征值进行排序
//   // ...
//
//   return error_code;
// }

} // namespace MathLib

