import numpy as np
import scipy
import scipy.sparse


def filter_array(arr: np.ndarray, discard: np.ndarray) -> np.ndarray:
    return arr[~np.isin(arr, discard)]


def convolve(arr: np.ndarray, kernel: np.ndarray) -> np.ndarray:
    return scipy.signal.convolve2d(arr, kernel, mode="same", boundary="fill", fillvalue=0)


class FEA:
    def __init__(
        self,
        num_elements_x,
        num_elements_y,
        volume_fraction,
        penalization,
        radius_min,
        move,
        young_modulus,
        young_modulus_min,
        poisson_ratio,
    ) -> None:

        self.num_elements_x = num_elements_x
        self.num_elements_y = num_elements_y
        self.num_elements = num_elements_x * num_elements_y
        self.num_nodes_x = num_elements_x + 1
        self.num_nodes_y = num_elements_y + 1
        self.num_nodes = self.num_nodes_x * self.num_nodes_y
        self.num_dofs_per_node = 2
        self.num_dofs = self.num_nodes * self.num_dofs_per_node

        self.volume_fraction = volume_fraction
        self.penalization = penalization
        self.move = move
        self.young_modulus = young_modulus
        self.young_modulus_min = young_modulus_min

        self.node_indices = np.arange(self.num_nodes, dtype=np.int32).reshape(
            (self.num_nodes_y, self.num_nodes_x)
        )

        # Represents the index of the first DOF of each element
        connectivity_vector = self.num_dofs_per_node * self.node_indices[
            :num_elements_y, :num_elements_x
        ].reshape((self.num_elements, 1))

        # Each row of the connectivity matrix indexes the 8 DOFs of the
        # corresponding element
        self.connectivity_matrix = connectivity_vector + np.array(
            [
                2,
                3,
                self.num_dofs_per_node * self.num_nodes_y + 2,
                self.num_dofs_per_node * self.num_nodes_y + 3,
                self.num_dofs_per_node * self.num_nodes_y + 0,
                self.num_dofs_per_node * self.num_nodes_y + 1,
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
        stiffness_matrix_indices_i = (
            self.connectivity_matrix[:, dof_connectivities_i].transpose().flatten()
        )
        stiffness_matrix_indices_j = (
            self.connectivity_matrix[:, dof_connectivities_j].transpose().flatten()
        )
        self.stiffness_matrix_indices = np.stack(
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
        self.element_stiffness_matrix_values = (
            1.0
            / (24.0 * (1.0 - poisson_ratio * poisson_ratio))
            * (
                element_stiffness_matrix_values_1
                + poisson_ratio * element_stiffness_matrix_values_2
            )
        ).astype(np.float32)
        self.element_stiffness_matrix = np.zeros((8, 8), dtype=np.float32)
        index = 0
        for j in range(8):
            for i in range(j, 8):
                self.element_stiffness_matrix[i, j] = self.element_stiffness_matrix_values[index]
                self.element_stiffness_matrix[j, i] = self.element_stiffness_matrix_values[index]
                index += 1

        self.young_moduli = np.zeros(self.num_elements, dtype=np.float32)

        passive_solid = np.array([], dtype=np.int32)
        passive_void = np.array([], dtype=np.int32)
        passive_elements = np.concatenate([passive_solid, passive_void], dtype=np.int32)
        self.active_elements = filter_array(
            np.arange(self.num_elements, dtype=np.int32), passive_elements
        )

        self.fixed_dofs = np.concatenate(
            [
                np.linspace(
                    0,
                    self.num_dofs_per_node * self.num_nodes_y,
                    num=self.num_nodes_y,
                    dtype=np.int32,
                ),
                [self.num_dofs_per_node * self.node_indices[-1, -1] + 1],
            ],
            dtype=np.int32,
        )
        self.free_dofs = filter_array(np.arange(self.num_dofs, dtype=np.int32), self.fixed_dofs)

        # Maps DOF indices to indices in the global stiffness matrix. The value -1
        # indicates that the corresponding DOF is fixed and hence not present in
        # the stiffness matrix
        self.all_to_free = np.full(self.num_dofs, -1, dtype=np.int32)
        self.all_to_free[self.free_dofs] = np.arange(self.free_dofs.size, dtype=np.int32)

        self.forces = np.zeros(self.free_dofs.size, dtype=np.float32)
        self.forces[self.all_to_free[self.num_dofs_per_node * self.node_indices[0, 0] + 1]] = -1.0

        self.displacements = np.zeros(self.num_dofs, dtype=np.float32)

        kernel_size = 2 * int(np.ceil(radius_min)) - 1
        kernel_min_coord = -np.ceil(radius_min) + 1.0
        kernel_iy, kernel_ix = np.indices((kernel_size, kernel_size), dtype=np.float32)
        kernel_y = kernel_min_coord + kernel_iy
        kernel_x = kernel_min_coord + kernel_ix
        kernel_r = np.hypot(kernel_x, kernel_y)
        self.filter_kernel = np.maximum(radius_min - kernel_r, 0.0).astype(np.float32)
        self.filter_weights = convolve(
            np.ones((num_elements_y, num_elements_x), dtype=np.float32), self.filter_kernel
        )

        self.design_variables = np.zeros(self.num_elements, dtype=np.float32)
        # TODO: check this
        self.design_variables[self.active_elements] = (
            volume_fraction * (self.num_elements - passive_solid.size) - passive_solid.size
        ) / self.active_elements.size
        self.design_variables[passive_solid] = 1.0
        self.design_variables_filtered = np.zeros(self.num_elements, dtype=np.float32)
        self.design_variables_physical = self.design_variables.copy()
        self.design_variables_old = np.ones(self.num_elements, dtype=np.float32)

        self.displacement_matrix = np.zeros(self.connectivity_matrix.shape, dtype=np.float32)

        self.stiffness_derivative = np.zeros(self.num_elements, dtype=np.float32)
        self.compliance_derivative = np.zeros(self.num_elements, dtype=np.float32)
        self.filtered_compliance_derivative = np.zeros(self.num_elements, dtype=np.float32)
        self.volume_derivative = np.zeros(self.num_elements, dtype=np.float32)
        # TODO: check this
        self.volume_derivative[self.active_elements] = 1.0 / (self.num_elements * volume_fraction)
        self.filtered_volume_derivative = np.zeros(self.num_elements, dtype=np.float32)

        self.active_design_variables = np.zeros(self.active_elements.size, dtype=np.float32)
        self.lower_bound = np.zeros(self.active_elements.size, dtype=np.float32)
        self.upper_bound = np.zeros(self.active_elements.size, dtype=np.float32)
        self.resizing_rule_constant = np.zeros(self.active_elements.size, dtype=np.float32)

    def optimize(self, steps: int) -> None:
        for _ in range(steps):
            self._optimization_step()

    def _optimization_step(self) -> None:

        self.design_variables_filtered[:] = (
            convolve(
                self.design_variables.reshape((self.num_elements_y, self.num_elements_x)),
                self.filter_kernel,
            )
            / self.filter_weights
        ).flatten()

        self.design_variables_physical[self.active_elements] = self.design_variables_filtered[
            self.active_elements
        ]

        self.design_variables_old[:] = self.design_variables_physical

        self.young_moduli[:] = self.young_modulus_min + np.power(
            self.design_variables_physical, self.penalization
        ) * (self.young_modulus - self.young_modulus_min)
        self.stiffness_derivative[self.active_elements] = (
            -self.penalization
            * (self.young_modulus - self.young_modulus_min)
            * np.power(
                self.design_variables_physical[self.active_elements], self.penalization - 1.0
            )
        )
        # FIXME: pre-allocate
        self.stiffness_matrix_values = (
            self.element_stiffness_matrix_values[:, None] * self.young_moduli[None, :]
        ).flatten()

        self._solve_equilibrium_system()

        self.displacement_matrix = self.displacements[self.connectivity_matrix.flatten()].reshape(
            self.connectivity_matrix.shape
        )
        self.compliance_derivative = self.stiffness_derivative @ (
            (self.displacement_matrix @ self.element_stiffness_matrix) * self.displacement_matrix
        ).sum(axis=1)

        self.filtered_compliance_derivative = convolve(
            self.compliance_derivative.reshape((self.num_elements_y, self.num_elements_x))
            / self.filter_weights,
            self.filter_kernel,
        ).flatten()
        self.filtered_volume_derivative = convolve(
            self.volume_derivative.reshape((self.num_elements_y, self.num_elements_x))
            / self.filter_weights,
            self.filter_kernel,
        ).flatten()

        self.active_design_variables = self.design_variables[self.active_elements]

        self.lower_bound = self.active_design_variables - self.move
        self.upper_bound = self.active_design_variables + self.move
        self.resizing_rule_constant = self.active_design_variables * np.sqrt(
            np.maximum(
                -self.filtered_compliance_derivative[self.active_elements]
                / self.filtered_volume_derivative[self.active_elements],
                0.0,
            )
        )

        # Initial estimate for LM
        l1 = 0.0
        l2 = self.resizing_rule_constant.mean() / self.volume_fraction
        # OC resizing rule
        while (l2 - l1) / (l2 + l1) > 1e-4:
            l_middle = 0.5 * (l1 + l2)
            self.design_variables[self.active_elements] = np.clip(
                np.clip(self.resizing_rule_constant / l_middle, self.lower_bound, self.upper_bound),
                0.0,
                1.0,
            )
            if self.design_variables.mean() > self.volume_fraction:
                l1 = l_middle
            else:
                l2 = l_middle

    def _solve_equilibrium_system(self) -> None:

        # FIXME
        index_i = self.all_to_free[self.stiffness_matrix_indices[:, 0]]
        index_j = self.all_to_free[self.stiffness_matrix_indices[:, 1]]
        mask = (index_i >= 0) & (index_j >= 0)
        index_i = index_i[mask]
        index_j = index_j[mask]
        stiffness_matrix_values = self.stiffness_matrix_values[mask]

        stiffness_matrix = scipy.sparse.csc_matrix(
            (stiffness_matrix_values, (index_i, index_j)),
            shape=(self.free_dofs.size, self.free_dofs.size),
        )
        stiffness_matrix.prune()

        free_displacements = scipy.sparse.linalg.spsolve(stiffness_matrix, self.forces)

        self.displacements[:] = 0
        self.displacements[self.free_dofs] = free_displacements


def main() -> None:
    fea = FEA(
        num_elements_x=200,
        num_elements_y=100,
        volume_fraction=0.2,
        penalization=3.0,
        radius_min=2.0,
        move=0.2,
        young_modulus=1.0,
        young_modulus_min=1e-9,
        poisson_ratio=0.3,
    )
    fea.optimize(steps=100)


if __name__ == "__main__":
    main()
