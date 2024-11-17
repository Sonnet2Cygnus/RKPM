import numpy as np

# 计算窗函数
def window_function_d(r):
    abs_r = np.abs(r)
    if 0.5 <= abs_r <= 1.5:
        return (1/6) * (5 - 3*abs_r - np.sqrt(-3*(1-abs_r)**2 + 1))
    elif abs_r <= 0.5:
        return (1/3) * (1 + np.sqrt(-3*r**2 + 1))
    else:
        return 0

# 计算矩阵m
def compute_m_ab_matrix(S_I, nearest_grid_point, dx, dy):
    x_i, y_j = nearest_grid_point
    delta_A_mn = dx * dy
    m_ab_matrix = np.zeros((6, 6))

    for x_mn, y_mn in S_I:
        delta_x = (x_mn - x_i) / dx
        delta_y = (y_mn - y_j) / dy

        w_total = window_function_d(delta_x) * window_function_d(delta_y) * delta_A_mn

        powers_x = [1, delta_x, delta_x**2, delta_x**3, delta_x**4]
        powers_y = [1, delta_y, delta_y**2, delta_y**3, delta_y**4]

        for i in range(6):
            for j in range(6):
                m_ab_matrix[i, j] += powers_x[min(i, 4)] * powers_y[min(j, 4)] * w_total

    return m_ab_matrix

# 计算HI
def compute_H_I(delta_I, eta_I):

    H_I = np.diag([
        1,
        1 / delta_I,
        1 / eta_I,
        1 / (delta_I * eta_I),
        1 / (delta_I ** 2),
        1 / (eta_I ** 2)
    ])
    return H_I

# 计算bI 和 dI
def compute_b_I(H_I, M_I):
    d_I = np.zeros(6)
    d_I[0] = 1
    return d_I


# 计算修正窗函数
def modified_window_function(S_I, nearest_grid_point, d_I, dx, dy):

    x_i, y_j = nearest_grid_point

    modified_w_values = []

    for x_mn, y_mn in S_I:
        delta_x = (x_mn - x_i) / dx
        delta_y = (y_mn - y_j) / dy

        # 计算窗函数
        w_total = window_function_d(delta_x) * window_function_d(delta_y)

        # 计算修正窗函数
        modified_w = d_I[0] * w_total + \
                     d_I[1] * delta_x * w_total + \
                     d_I[2] * delta_y * w_total + \
                     d_I[3] * delta_x * delta_y * w_total + \
                     d_I[4] * delta_x**2 * w_total + \
                     d_I[5] * delta_y**2 * w_total

        modified_w_values.append(modified_w)

    return modified_w_values

# 对所有的点计算修正窗函数
def compute_all_modified_window_functions(all_S_I, lagrangian_points, nearest_euler_points, delta_I, eta_I, dx, dy):

    all_modified_w = []
    for idx in range(len(lagrangian_points)):
        S_I = all_S_I[idx]
        nearest_grid_point = lagrangian_points[idx]
        delta_I_lag = delta_I[idx]
        eta_I_lag = eta_I[idx]

        # 计算 m_ab_matrix
        M_I = compute_m_ab_matrix(S_I, nearest_grid_point, dx, dy)

        # 计算 H_I
        H_I = compute_H_I(delta_I_lag, eta_I_lag)

        # 计算 b_I 和 d_I
        d_I = compute_b_I(H_I, M_I)

        # 计算修正窗函数
        modified_w = modified_window_function(S_I, nearest_grid_point, d_I, dx, dy)

        all_modified_w.append(modified_w)

    return all_modified_w