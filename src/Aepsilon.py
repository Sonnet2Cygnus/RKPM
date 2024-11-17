import numpy as np
import math
from src import window

# 计算弧长
def compute_delta_s_I(radius, Ne):

    central_angle = 2 * math.pi / Ne
    delta_s_I = radius * central_angle
    return delta_s_I

# 计算矩阵A
def compute_A_matrix(delta_s_I, all_modified_w, dx, dy, lagrangian_points, delta_I, eta_I, all_S_I):

    Delta_A = dx * dy
    Ne = len(all_modified_w)
    
    # 初始化 A 矩阵
    A = np.zeros((Ne, Ne))
    
    for I in range(Ne):
        S_I = all_S_I[I]
        delta_I_lag = delta_I[I]
        eta_I_lag = eta_I[I]

        for K in range(Ne):
            # 计算 a_IK = delta_s_I * sum(w_I * w_K) * Delta_A
            nearest_grid_point = lagrangian_points[K]
            M_I = window.compute_m_ab_matrix(S_I, nearest_grid_point, dx, dy)
            H_I = window.compute_H_I(delta_I_lag, eta_I_lag)
            d_I = window.compute_b_I(M_I, H_I)
            w_k = window.modified_window_function(S_I, nearest_grid_point, d_I, dx, dy)
            sum_w = np.sum(np.array(all_modified_w[I] * np.array(w_k)))
            A[I, K] = delta_s_I * sum_w * Delta_A
    
    return A


# 计算epsilon
def compute_epsilon(A, Ne, dx, max_iterations=100, tolerance=1e-6):

    # 初始化epsilon
    epsilon = np.full(A.shape[0], 2.0 * np.pi * dx / Ne)

    # 使用迭代的方法 计算 
    for _ in range(max_iterations):
        residual = np.ones(A.shape[0]) - A @ epsilon
        if np.linalg.norm(residual) < tolerance:
            break

        # 更新残差
        epsilon += np.linalg.solve(A, residual)

    return epsilon

