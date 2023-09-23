#include "fea.hpp"

#include <Eigen/SparseCholesky>

#include <chrono>
#include <iostream>
#include <sstream>

namespace
{

[[nodiscard]] inline std::chrono::steady_clock::time_point timer_start()
{
    return std::chrono::steady_clock::now();
}

inline void timer_stop(std::chrono::steady_clock::time_point start,
                       const char *label)
{
    const auto end = std::chrono::steady_clock::now();
    std::cout << label << ": "
              << std::chrono::duration<double>(end - start).count() * 1'000.0
              << " ms\n";
}

[[nodiscard]] constexpr const char *
to_string(Eigen::ComputationInfo computation_info) noexcept
{
    switch (computation_info)
    {
    case Eigen::Success: return "Success";
    case Eigen::NumericalIssue: return "NumericalIssue";
    case Eigen::NoConvergence: return "NoConvergence";
    case Eigen::InvalidInput: return "InvalidInput";
    }
    return "UNDEFINED";
}

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
    const auto t = timer_start();

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
            .reshaped(num_elements * 36, 1)};
    const Eigen::VectorXi stiffness_matrix_indices_j {
        connectivity_matrix(Eigen::all, dof_connectivities_j)
            .transpose()
            .reshaped(num_elements * 36, 1)};

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

    const Eigen::Vector<float, 36> element_stiffness_matrix_values_1 {
        12, 3,  -6, -3, -6, -3, 0, 3,  12, 3, 0,  -3, -6, -3, -6, 12, -3, 0,
        -3, -6, 3,  12, 3,  -6, 3, -6, 12, 3, -6, -3, 12, 3,  0,  12, -3, 12};
    const Eigen::Vector<float, 36> element_stiffness_matrix_values_2 {
        -4, 3, -2, 9,  2,  -3, 4, -9, -4, -9, 4,  -3, 2,  9,  -2, -4, -3, 4,
        9,  2, 3,  -4, -9, -2, 3, 2,  -4, 3,  -2, 9,  -4, -9, 4,  -4, -3, -4};
    problem.element_stiffness_matrix_values =
        1.0f / (1.0f - poisson_ratio * poisson_ratio) / 24.0f *
        (element_stiffness_matrix_values_1 +
         poisson_ratio * element_stiffness_matrix_values_2);

    problem.young_moduli.setConstant(num_elements, young_modulus);

    Eigen::VectorXi fixed_dofs(num_nodes_y + 1);
    fixed_dofs << Eigen::VectorXi::LinSpaced(
        num_nodes_y, 0, num_dofs_per_node * num_nodes_y - 1),
        num_dofs_per_node * node_indices(Eigen::last, Eigen::last) + 1;
    const Eigen::VectorXi all_dofs {
        Eigen::VectorXi::LinSpaced(num_dofs, 0, num_dofs - 1)};
    problem.free_dofs = set_difference(all_dofs, fixed_dofs);

    problem.forces.setZero(num_dofs);
    problem.forces(num_dofs_per_node * node_indices(0, 0) + 1) = -1.0f;

    timer_stop(t, "fea_init");

    return problem;
}

Eigen::VectorXf fea_solve(const FEA_problem &problem)
{
    const auto global_t = timer_start();
    auto t = timer_start();

    const Eigen::VectorXf stiffness_matrix_values {
        (problem.element_stiffness_matrix_values *
         problem.young_moduli.transpose())
            .reshaped()};

    // Maps DOF indices to indices in the global stiffness matrix. The value -1
    // indicates that the corresponding DOF is fixed and hence not present in
    // the stiffness matrix
    Eigen::VectorXi all_to_free(problem.num_dofs);
    int current_stiffness_matrix_index {0};
    for (Eigen::Index i {0}; i < problem.num_dofs; ++i)
    {
        if (std::binary_search(
                problem.free_dofs.cbegin(), problem.free_dofs.cend(), i))
        {
            all_to_free(i) = current_stiffness_matrix_index;
            ++current_stiffness_matrix_index;
        }
        else
        {
            all_to_free(i) = -1;
        }
    }

    const Eigen::Index num_values {stiffness_matrix_values.size()};
    std::vector<Eigen::Triplet<float>> triplets;
    triplets.reserve(static_cast<std::size_t>(num_values));
    for (Eigen::Index i {0}; i < num_values; ++i)
    {
        const auto mapped_i =
            all_to_free(problem.stiffness_matrix_indices(i, 0));
        const auto mapped_j =
            all_to_free(problem.stiffness_matrix_indices(i, 1));
        if (mapped_i != -1 && mapped_j != -1)
        {
            triplets.emplace_back(
                mapped_i, mapped_j, stiffness_matrix_values(i));
        }
    }

    Eigen::SparseMatrix<float> stiffness_matrix(problem.free_dofs.size(),
                                                problem.free_dofs.size());
    stiffness_matrix.setFromTriplets(triplets.cbegin(), triplets.cend());
    stiffness_matrix.prune(0.0f, 0.0f);

    timer_stop(t, "stiffness matrix assembly");
    t = timer_start();

    Eigen::SimplicialLDLT<Eigen::SparseMatrix<float>, Eigen::Lower> solver;
    solver.compute(stiffness_matrix);
    if (const auto result = solver.info(); result != Eigen::Success)
    {
        std::ostringstream message;
        message << "Decomposition failed: " << to_string(result);
        throw std::runtime_error(message.str());
    }

    timer_stop(t, "stiffness matrix decomposition");
    t = timer_start();

    const Eigen::VectorXf free_displacements {
        solver.solve(problem.forces(problem.free_dofs))};
    Eigen::VectorXf displacements(problem.num_dofs);
    displacements.setZero();
    displacements(problem.free_dofs) = free_displacements;
    if (const auto result = solver.info(); result != Eigen::Success)
    {
        std::ostringstream message;
        message << "Solving failed: " << to_string(result);
        throw std::runtime_error(message.str());
    }

    timer_stop(t, "solving system");
    timer_stop(global_t, "fea_solve");

    return displacements;
}
