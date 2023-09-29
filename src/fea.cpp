#include "fea.hpp"

#include <Eigen/SparseCholesky>
#include <Eigen/SparseCore>

#include <chrono>
#include <fstream>
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

// discard must be in ascending order and not contain duplicates
[[nodiscard]] Eigen::VectorXi
filtered_index_vector(int size, const Eigen::VectorXi &discard)
{
    Eigen::VectorXi result(size - discard.size());
    Eigen::Index index_result {0};
    Eigen::Index index_discard {0};

    for (int i {0}; i < size; ++i)
    {
        if (index_discard < discard.size() && i == discard(index_discard))
        {
            ++index_discard;
        }
        else
        {
            result(index_result) = i;
            ++index_result;
        }
    }

    return result;
}

void load_densities(Eigen::VectorXf &densities, const char *file_name)
{
    assert(densities.size() == 100 * 200);

    std::ifstream file(file_name);
    if (!file)
    {
        std::cerr << "Failed to open \"" << file_name << "\"\n";
        return;
    }

    for (auto &value : densities)
    {
        file >> value;
    }
}

} // namespace

FEA_state fea_init(int num_elements_x, int num_elements_y)
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

    FEA_state fea {};
    fea.num_elements_x = num_elements_x;
    fea.num_elements_y = num_elements_y;
    fea.num_elements = num_elements;
    fea.num_nodes_x = num_nodes_x;
    fea.num_nodes_y = num_nodes_y;
    fea.num_nodes = num_nodes;
    fea.num_dofs = num_dofs;
    fea.num_dofs_per_node = num_dofs_per_node;

    fea.stiffness_matrix_indices.resize(stiffness_matrix_indices_i.rows(), 2);
    fea.stiffness_matrix_indices
        << stiffness_matrix_indices_i.cwiseMax(stiffness_matrix_indices_j),
        stiffness_matrix_indices_i.cwiseMin(stiffness_matrix_indices_j);

    const Eigen::Vector<float, 36> element_stiffness_matrix_values_1 {
        12, 3,  -6, -3, -6, -3, 0, 3,  12, 3, 0,  -3, -6, -3, -6, 12, -3, 0,
        -3, -6, 3,  12, 3,  -6, 3, -6, 12, 3, -6, -3, 12, 3,  0,  12, -3, 12};
    const Eigen::Vector<float, 36> element_stiffness_matrix_values_2 {
        -4, 3, -2, 9,  2,  -3, 4, -9, -4, -9, 4,  -3, 2,  9,  -2, -4, -3, 4,
        9,  2, 3,  -4, -9, -2, 3, 2,  -4, 3,  -2, 9,  -4, -9, 4,  -4, -3, -4};
    fea.element_stiffness_matrix_values =
        1.0f / 24.0f / (1.0f - poisson_ratio * poisson_ratio) *
        (element_stiffness_matrix_values_1 +
         poisson_ratio * element_stiffness_matrix_values_2);
    Eigen::Index index {0};
    for (Eigen::Index j {0}; j < 8; ++j)
    {
        for (Eigen::Index i {j}; i < 8; ++i)
        {
            fea.element_stiffness_matrix(i, j) =
                fea.element_stiffness_matrix_values(index);
            fea.element_stiffness_matrix(j, i) =
                fea.element_stiffness_matrix_values(index);
            ++index;
        }
    }

    fea.young_moduli.setConstant(num_elements, young_modulus);
    //load_densities(fea.young_moduli, "../densities.txt");

    fea.passive_solid = {};
    fea.passive_void = {};
    Eigen::VectorXi passive_elements(fea.passive_solid.size() +
                                     fea.passive_void.size());
    passive_elements << fea.passive_solid, fea.passive_void;
    fea.active_elements = filtered_index_vector(num_elements, passive_elements);

    Eigen::VectorXi fixed_dofs(num_nodes_y + 1);
    fixed_dofs << Eigen::VectorXi::LinSpaced(
        num_nodes_y, 0, num_dofs_per_node * num_nodes_y - 1),
        num_dofs_per_node * node_indices(Eigen::last, Eigen::last) + 1;
    fea.free_dofs = filtered_index_vector(num_dofs, fixed_dofs);

    // Maps DOF indices to indices in the global stiffness matrix. The value -1
    // indicates that the corresponding DOF is fixed and hence not present in
    // the stiffness matrix
    fea.all_to_free.resize(fea.num_dofs);
    int current_stiffness_matrix_index {0};
    for (Eigen::Index i {0}; i < fea.num_dofs; ++i)
    {
        // TODO: this can probably be optimized just like
        // filtered_index_vector
        if (std::binary_search(fea.free_dofs.cbegin(), fea.free_dofs.cend(), i))
        {
            fea.all_to_free(i) = current_stiffness_matrix_index;
            ++current_stiffness_matrix_index;
        }
        else
        {
            fea.all_to_free(i) = -1;
        }
    }

    fea.forces.setZero(fea.free_dofs.size());
    fea.forces(fea.all_to_free(num_dofs_per_node * node_indices(0, 0) + 1)) =
        -1.0f;

    fea.displacements.setZero(num_dofs);

    timer_stop(t, "fea_init");

    return fea;
}

void fea_solve(FEA_state &fea)
{
    const auto global_t = timer_start();
    auto t = timer_start();

    const Eigen::VectorXf stiffness_matrix_values {
        (fea.element_stiffness_matrix_values * fea.young_moduli.transpose())
            .reshaped()};

    const Eigen::Index num_values {stiffness_matrix_values.size()};
    std::vector<Eigen::Triplet<float>> triplets;
    triplets.reserve(static_cast<std::size_t>(num_values));
    for (Eigen::Index i {0}; i < num_values; ++i)
    {
        const auto mapped_i =
            fea.all_to_free(fea.stiffness_matrix_indices(i, 0));
        const auto mapped_j =
            fea.all_to_free(fea.stiffness_matrix_indices(i, 1));
        if (mapped_i != -1 && mapped_j != -1)
        {
            triplets.emplace_back(
                mapped_i, mapped_j, stiffness_matrix_values(i));
        }
    }

    Eigen::SparseMatrix<float> stiffness_matrix(fea.free_dofs.size(),
                                                fea.free_dofs.size());
    stiffness_matrix.setFromTriplets(triplets.cbegin(), triplets.cend());
    stiffness_matrix.prune(0.0f);

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

    Eigen::VectorXf free_displacements(fea.free_dofs.size());
    free_displacements = solver.solve(fea.forces);
    fea.displacements.setZero();
    fea.displacements(fea.free_dofs) = free_displacements;

    if (const auto result = solver.info(); result != Eigen::Success)
    {
        std::ostringstream message;
        message << "Solving failed: " << to_string(result);
        throw std::runtime_error(message.str());
    }

    timer_stop(t, "solving system");
    timer_stop(global_t, "fea_solve");
}

void fea_init_optimization(FEA_state &fea,
                           float volume_fraction,
                           float penalization,
                           float radius_min,
                           float move)
{
}

void fea_optimization_step(FEA_state &fea)
{
}