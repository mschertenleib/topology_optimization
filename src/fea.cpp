#include "fea.hpp"

#include <Eigen/SparseCholesky>

#include <iostream>

namespace
{

// FIXME: this is temporary and not very safe
// Assumes all elements of b are contained in a
Eigen::VectorXi set_difference(const Eigen::VectorXi &a,
                               const Eigen::VectorXi &b)
{
    Eigen::VectorXi result(a.size() - b.size());
    Eigen::Index index_result {0};

    for (const auto value_a : a)
    {
        if (std::none_of(b.cbegin(),
                         b.cend(),
                         [value_a](Eigen::Index value_b)
                         { return value_a == value_b; }))
        {
            result(index_result) = value_a;
            ++index_result;
        }
    }

    return result;
}

} // namespace

FEA_problem fea_init(int num_elements_x, int num_elements_y)
{
    constexpr float young_modulus {1.0f};
    constexpr float young_modulus_min {1e-9f};
    constexpr float poisson_ratio {0.3f};

    const int num_elements {num_elements_x * num_elements_y};
    const int num_nodes_x {num_elements_x + 1};
    const int num_nodes_y {num_elements_y + 1};
    const int num_nodes {num_nodes_x * num_nodes_y};
    const int num_dofs_per_node {2};
    const int num_dofs {num_nodes * num_dofs_per_node};

    const Eigen::MatrixXi node_indices {
        Eigen::VectorXi::LinSpaced(num_nodes, 0, num_nodes - 1)
            .reshaped(num_nodes_y, num_nodes_x)};

    // Represents the index of the first DOF of each element
    const Eigen::VectorXi connectivity_vector {
        num_dofs_per_node *
        node_indices.topLeftCorner(num_elements_y, num_elements_x)
            .reshaped(num_elements, 1)};

    // Each row of the connectivity matrix indexes the 8 DOFs of the
    // corresponding element
    const Eigen::Matrix<int, Eigen::Dynamic, 8> connectivity_matrix {
        connectivity_vector.replicate(1, 8).rowwise() +
        Eigen::RowVector<int, 8> {2,
                                  3,
                                  num_dofs_per_node * num_nodes_y + 2,
                                  num_dofs_per_node * num_nodes_y + 3,
                                  num_dofs_per_node * num_nodes_y + 0,
                                  num_dofs_per_node * num_nodes_y + 1,
                                  0,
                                  1}};

    const Eigen::Vector<int, 36> dof_connectivities_i {
        0, 1, 2, 3, 4, 5, 6, 7, 1, 2, 3, 4, 5, 6, 7, 2, 3, 4,
        5, 6, 7, 3, 4, 5, 6, 7, 4, 5, 6, 7, 5, 6, 7, 6, 7, 7};
    const Eigen::Vector<int, 36> dof_connectivities_j {
        0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2,
        2, 2, 2, 3, 3, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 6, 6, 7};
    const Eigen::VectorXi stiffness_matrix_indices_i {
        connectivity_matrix(Eigen::all, dof_connectivities_i)
            .transpose()
            .reshaped()};
    const Eigen::VectorXi stiffness_matrix_indices_j {
        connectivity_matrix(Eigen::all, dof_connectivities_j)
            .transpose()
            .reshaped()};

    FEA_problem problem {};
    problem.num_elements_x = num_elements_x;
    problem.num_elements_y = num_elements_y;
    problem.num_elements = num_elements;
    problem.num_nodes_x = num_nodes_x;
    problem.num_nodes_y = num_nodes_y;
    problem.num_nodes = num_nodes;
    problem.num_dofs = num_dofs;

    problem.stiffness_matrix_indices.resize(stiffness_matrix_indices_i.rows(),
                                            2);
    problem.stiffness_matrix_indices
        << stiffness_matrix_indices_i.cwiseMax(stiffness_matrix_indices_j),
        stiffness_matrix_indices_i.cwiseMin(stiffness_matrix_indices_j);

    const Eigen::Vector<float, 36> element_stiffness_matrix_coefficients_1 {
        12, 3,  -6, -3, -6, -3, 0, 3,  12, 3, 0,  -3, -6, -3, -6, 12, -3, 0,
        -3, -6, 3,  12, 3,  -6, 3, -6, 12, 3, -6, -3, 12, 3,  0,  12, -3, 12};
    const Eigen::Vector<float, 36> element_stiffness_matrix_coefficients_2 {
        -4, 3, -2, 9,  2,  -3, 4, -9, -4, -9, 4,  -3, 2,  9,  -2, -4, -3, 4,
        9,  2, 3,  -4, -9, -2, 3, 2,  -4, 3,  -2, 9,  -4, -9, 4,  -4, -3, -4};
    problem.element_stiffness_matrix_values =
        1.0f / (1.0f - poisson_ratio * poisson_ratio) / 24.0f *
        (element_stiffness_matrix_coefficients_1 +
         poisson_ratio * element_stiffness_matrix_coefficients_2);

    problem.young_moduli.setConstant(num_elements, young_modulus);

    Eigen::VectorXi fixed_dofs(num_nodes_y + 1);
    fixed_dofs << Eigen::VectorXi::LinSpaced(
        num_nodes_y, 0, num_dofs_per_node * num_nodes_y - 1),
        num_dofs_per_node * node_indices(Eigen::last, Eigen::last) + 1;
    const Eigen::VectorXi all_dofs {
        Eigen::VectorXi::LinSpaced(num_dofs, 0, num_dofs - 1)};
    problem.free_dofs = set_difference(all_dofs, fixed_dofs);

    problem.forces.resize(num_dofs);
    problem.forces.setZero();
    problem.forces(num_dofs_per_node * node_indices(0, 0) + 1) = -1.0f;
    // problem.forces.insert(num_dofs_per_node * node_indices(0, 0) + 1) =
    // -1.0f;

    return problem;
}

Eigen::VectorXf fea_solve(const FEA_problem &problem)
{
    const Eigen::VectorXf stiffness_matrix_values {
        (problem.element_stiffness_matrix_values *
         problem.young_moduli.transpose())
            .reshaped()};

    const Eigen::Index num_values {stiffness_matrix_values.size()};
    std::vector<Eigen::Triplet<float>> triplets(
        static_cast<std::size_t>(num_values));
    for (Eigen::Index i {0}; i < num_values; ++i)
    {
        triplets[static_cast<std::size_t>(i)] = {
            problem.stiffness_matrix_indices(i, 0),
            problem.stiffness_matrix_indices(i, 1),
            stiffness_matrix_values(i)};
    }

    Eigen::SparseMatrix<float> stiffness_matrix(problem.num_dofs,
                                                problem.num_dofs);
    stiffness_matrix.setFromTriplets(triplets.cbegin(), triplets.cend());
    // TODO: the stiffness matrix is correct, however some coefficients are zero
    // as a result of summing coefficients from triplets. We could probably gain
    // some performance by removing them

    Eigen::SimplicialLDLT<Eigen::SparseMatrix<float>, Eigen::Lower> solver;
    solver.compute(stiffness_matrix);
    if (solver.info() != Eigen::Success)
    {
        throw std::runtime_error("Decomposition failed");
    }

    // FIXME: we might have solved our uninitialized/out of bounds memory
    // problem by making the force vector dense, but we are still not filtering
    // fixed DOFs. If we do that, we might actually get a working system, and it
    // would then just be a case of making it work with a sparse force vector...
    Eigen::VectorXf displacements {solver.solve(problem.forces)};
    if (solver.info() != Eigen::Success)
    {
        throw std::runtime_error("Solving failed");
    }

    if (displacements.maxCoeff() > 15.0f || displacements.minCoeff() < -50.0f)
    {
        std::cout << "Unexpected displacements in range ["
                  << displacements.minCoeff() << ", "
                  << displacements.maxCoeff() << "]\n";
    }

    return displacements;
}
