import numpy as np
import scipy


# FIXME
# Discard must be in ascending order and not contain duplicates
def filtered_index_array(size: int, discard: np.ndarray) -> np.ndarray:
    discard = discard.flatten()

    result = np.zeros(size - discard.size, dtype=np.int32)
    index_result = 0
    index_discard = 0

    for i in range(size):
        if index_discard < discard.size and i == discard[index_discard]:
            index_discard += 1
        else:
            result[index_result] = i
            index_result += 1

    return result


def main() -> None:

    num_elements_x = 200
    num_elements_y = 100
    volume_fraction = 0.2
    penalization = 3.0
    radius_min = 2.0
    move = 0.2

    young_modulus = 1.0
    young_modulus_min = 1e-9
    poisson_ratio = 0.3

    num_elements = num_elements_x * num_elements_y
    num_nodes_x = num_elements_x + 1
    num_nodes_y = num_elements_y + 1
    num_nodes = num_nodes_x * num_nodes_y
    num_dofs_per_node = 2
    num_dofs = num_nodes * num_dofs_per_node

    node_indices = np.arange(num_nodes, dtype=np.int32).reshape((num_nodes_y, num_nodes_x))

    # Represents the index of the first DOF of each element
    connectivity_vector = num_dofs_per_node * node_indices[
        :num_elements_y, :num_elements_x
    ].reshape((num_elements, 1))

    # Each row of the connectivity matrix indexes the 8 DOFs of the
    # corresponding element
    connectivity_matrix = connectivity_vector + np.array(
        [
            2,
            3,
            num_dofs_per_node * num_nodes_y + 2,
            num_dofs_per_node * num_nodes_y + 3,
            num_dofs_per_node * num_nodes_y + 0,
            num_dofs_per_node * num_nodes_y + 1,
            0,
            1,
        ],
        dtype=np.int32,
    ).reshape((1, 8))

    # fmt: off
    dof_connectivities_i = np.array(
        [0, 1, 2, 3, 4, 5, 6, 7, 1, 2, 3, 4, 5, 6, 7, 2, 3, 4,
         5, 6, 7, 3, 4, 5, 6, 7, 4, 5, 6, 7, 5, 6, 7, 6, 7, 7],
        dtype=np.int32,
    )
    dof_connectivities_j = np.array(
        [0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2,
         2, 2, 2, 3, 3, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 6, 6, 7],
        dtype=np.int32,
    )
    # fmt: on
    stiffness_matrix_indices_i = connectivity_matrix[:, dof_connectivities_i].transpose().flatten()
    stiffness_matrix_indices_j = connectivity_matrix[:, dof_connectivities_j].transpose().flatten()
    stiffness_matrix_indices = np.stack(
        [
            np.maximum(stiffness_matrix_indices_i, stiffness_matrix_indices_j),
            np.minimum(stiffness_matrix_indices_i, stiffness_matrix_indices_j),
        ],
        axis=1,
    )

    # fmt: off
    element_stiffness_matrix_values_1 = np.array(
        [12, 3, -6, -3, -6, -3, 0, 3, 12, 3, 0, -3, -6, -3, -6, 12, -3, 0,
         -3, -6, 3, 12, 3, -6, 3, -6, 12, 3, -6, -3, 12, 3, 0, 12, -3, 12])
    element_stiffness_matrix_values_2 = np.array(
        [-4, 3, -2, 9, 2, -3, 4, -9, -4, -9, 4, -3, 2, 9, -2, -4, -3, 4,
         9, 2, 3, -4, -9, -2, 3, 2, -4, 3, -2, 9, -4, -9, 4, -4, -3, -4])
    # fmt: on
    element_stiffness_matrix_values = (
        1.0
        / (24.0 * (1.0 - poisson_ratio * poisson_ratio))
        * (element_stiffness_matrix_values_1 + poisson_ratio * element_stiffness_matrix_values_2)
    ).astype(np.float32)
    element_stiffness_matrix = np.zeros((8, 8), dtype=np.float32)
    index = 0
    for j in range(8):
        for i in range(j, 8):
            element_stiffness_matrix[i, j] = element_stiffness_matrix_values[index]
            element_stiffness_matrix[j, i] = element_stiffness_matrix_values[index]
            index += 1

    young_moduli = np.zeros(num_elements, dtype=np.float32)

    passive_solid = np.array([], dtype=np.int32)
    passive_void = np.array([], dtype=np.int32)
    passive_elements = np.concatenate([passive_solid, passive_void], dtype=np.int32)
    active_elements = filtered_index_array(num_elements, passive_elements)

    fixed_dofs = np.concatenate(
        [
            np.linspace(0, num_dofs_per_node * num_nodes_y, num=num_nodes_y, dtype=np.int32),
            [num_dofs_per_node * node_indices[-1, -1] + 1],
        ],
        dtype=np.int32,
    )
    free_dofs = filtered_index_array(num_dofs, fixed_dofs)

    # Maps DOF indices to indices in the global stiffness matrix. The value -1
    # indicates that the corresponding DOF is fixed and hence not present in
    # the stiffness matrix
    all_to_free = np.full(num_dofs, -1, dtype=np.int32)
    all_to_free[free_dofs] = np.arange(free_dofs.size, dtype=np.int32)

    forces = np.zeros(free_dofs.size, dtype=np.float32)
    forces[all_to_free[num_dofs_per_node * node_indices[0, 0] + 1]] = -1.0

    displacements = np.zeros(num_dofs, dtype=np.float32)

    kernel_size = 2 * int(np.ceil(radius_min)) - 1
    kernel_min_coord = -np.ceil(radius_min) + 1.0
    kernel_iy, kernel_ix = np.indices((kernel_size, kernel_size), dtype=np.float32)
    kernel_y = kernel_min_coord + kernel_iy
    kernel_x = kernel_min_coord + kernel_ix
    kernel_r = np.hypot(kernel_x, kernel_y)
    filter_kernel = np.maximum(radius_min - kernel_r, 0.0).astype(np.float32)

    filter_weights = scipy.signal.convolve2d(
        np.ones((num_elements_y, num_elements_x), dtype=np.float32),
        filter_kernel,
        mode="same",
        boundary="fill",
        fillvalue=0,
    )

    design_variables = np.zeros(num_elements, dtype=np.float32)
    # TODO: check this
    design_variables[active_elements] = (
        volume_fraction * (num_elements - passive_solid.size) - passive_solid.size
    ) / active_elements.size
    design_variables[passive_solid] = 1.0
    design_variables_filtered = np.zeros(num_elements, dtype=np.float32)
    design_variables_physical = design_variables.copy()
    design_variables_old = np.ones(num_elements, dtype=np.float32)
    design_variables_indexed_temp = np.zeros(active_elements.size, dtype=np.float32)

    displacement_matrix = np.zeros(connectivity_matrix.shape, dtype=np.float32)

    stiffness_derivative = np.zeros(num_elements, dtype=np.float32)
    compliance_derivative = np.zeros(num_elements, dtype=np.float32)
    filtered_compliance_derivative = np.zeros(num_elements, dtype=np.float32)
    volume_derivative = np.zeros(num_elements, dtype=np.float32)
    # TODO: check this
    volume_derivative[active_elements] = 1.0 / (num_elements * volume_fraction)
    filtered_volume_derivative = np.zeros(num_elements, dtype=np.float32)

    active_design_variables = np.zeros(active_elements.size, dtype=np.float32)
    lower_bound = np.zeros(active_elements.size, dtype=np.float32)
    upper_bound = np.zeros(active_elements.size, dtype=np.float32)
    resizing_rule_constant = np.zeros(active_elements.size, dtype=np.float32)


if __name__ == "__main__":
    main()
