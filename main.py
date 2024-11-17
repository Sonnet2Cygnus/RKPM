import numpy as np
import matplotlib.pyplot as plt
from src import SI_generated, window, Aepsilon, error

def main():

    # 设置测试数据参数
    dx = 0.05063
    dy = 0.05063
    Ne = 200
    radius = 1.0
    center = (2.5, 2.5)
    grid_size = 5.0
    Delta_A = dx * dy
    # 测试
    eulerian_points, lagrangian_points, nearest_grid_points, delta_I, eta_I, all_S_I = SI_generated.generate_grid(
        dx, dy, Ne, radius, center, grid_size
    )

    all_modified_w = window.compute_all_modified_window_functions(
            all_S_I, lagrangian_points, nearest_grid_points, delta_I, eta_I, dx, dy
    )

    delta_s_I = Aepsilon.compute_delta_s_I(radius, Ne)
    print(delta_s_I)

    A_matrix = Aepsilon.compute_A_matrix(
            delta_s_I, all_modified_w, dx, dy, lagrangian_points, delta_I, eta_I, all_S_I
    )
    print("A Matrix:")
    print(A_matrix)

    epsilon = Aepsilon.compute_epsilon(A_matrix, Ne, dx)
    print(epsilon)
    print(f"... 共 {Ne} 个 ε 值")

    
    infinity_norm_error = error.compute_error(eulerian_points, lagrangian_points, all_S_I, all_modified_w, epsilon, delta_s_I, Delta_A)
    print(f"Error: {infinity_norm_error}")

if __name__ == "__main__":
    main()
