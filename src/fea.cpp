#include "fea.hpp"
#include "utility.hpp"

#include <Eigen/SparseCholesky>
#include <Eigen/SparseCore>

#include <cassert>
#include <iostream>
#include <sstream>

namespace
{

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
[[nodiscard]] Eigen::ArrayXi filtered_index_array(int size,
                                                  const Eigen::ArrayXi &discard)
{
    Eigen::ArrayXi result(size - discard.size());
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

// Returns a [m.size() by neighboring_rows * neighboring_cols] matrix, for which
// each row is the flattened block of elements surrounding the corresponding
// element in m, padding with 0 for out-of-bounds elements
[[nodiscard]] Eigen::MatrixXf neighborhood(const Eigen::MatrixXf &m,
                                           Eigen::Index neighboring_rows,
                                           Eigen::Index neighboring_cols)
{
    Eigen::MatrixXf result(m.size(), neighboring_rows * neighboring_cols);
    result.setZero();

    for (Eigen::Index j {0}; j < neighboring_cols; ++j)
    {
        for (Eigen::Index i {0}; i < neighboring_rows; ++i)
        {
            const auto shift = (j - neighboring_cols / 2) * m.rows() +
                               (i - neighboring_rows / 2);
            const auto result_col = j * neighboring_rows + i;
            const auto result_min_row = std::max(shift, Eigen::Index {0});
            const auto result_max_row =
                std::min(m.size() - 1 + shift, m.size() - 1);
            const auto m_min_index = std::max(-shift, Eigen::Index {0});
            const auto m_max_index =
                std::min(m.size() - 1 - shift, m.size() - 1);
            result(Eigen::seq(result_min_row, result_max_row), result_col) =
                m.reshaped()(Eigen::seq(m_min_index, m_max_index));
        }
    }

    return result;
}

[[nodiscard]] Eigen::ArrayXXf filter(const Eigen::ArrayXXf &m,
                                     const Eigen::ArrayXXf &kernel)
{
#if 0

// FIXME: this is broken and does not give the right result

    return (neighborhood(m, kernel.rows(), kernel.cols()) *
            kernel.matrix().reshaped())
        .reshaped(m.rows(), m.cols())
        .array();

#else

    Eigen::ArrayXXf result;
    result.resizeLike(m);

    for (Eigen::Index i {0}; i < result.rows(); ++i)
    {
        for (Eigen::Index j {0}; j < result.cols(); ++j)
        {
            float sum {0.0f};
            for (Eigen::Index k {0}; k < kernel.rows(); ++k)
            {
                const auto m_i = i + k - kernel.rows() / 2;
                if (m_i >= 0 && m_i < m.rows())
                {
                    for (Eigen::Index l {0}; l < kernel.cols(); ++l)
                    {
                        const auto m_j = j + l - kernel.cols() / 2;
                        if (m_j >= 0 && m_j < m.cols())
                        {
                            sum += m(m_i, m_j) * kernel(k, l);
                        }
                    }
                }
            }
            result(i, j) = sum;
        }
    }

    return result;

#endif
}

void solve_equilibrium_system(FEA_state &fea)
{
    Scope_profiler prof_solve_system("Solve equilibrium system");
    Scope_profiler prof_assembly("Stiffness matrix assembly");

    const Eigen::Index num_values {fea.stiffness_matrix_values.size()};
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
                mapped_i, mapped_j, fea.stiffness_matrix_values(i));
        }
    }

    Eigen::SparseMatrix<float> stiffness_matrix(fea.free_dofs.size(),
                                                fea.free_dofs.size());
    stiffness_matrix.setFromTriplets(triplets.cbegin(), triplets.cend());
    stiffness_matrix.prune(0.0f);

    prof_assembly.stop();

    Eigen::SimplicialLDLT<Eigen::SparseMatrix<float>, Eigen::Lower> solver;
    {
        Scope_profiler prof_analyse_pattern("Pattern analysis");
        solver.compute(stiffness_matrix);
    }
    {
        Scope_profiler prof_factorize("Factorization");
        solver.factorize(stiffness_matrix);
    }

    if (const auto result = solver.info(); result != Eigen::Success)
    {
        std::ostringstream message;
        message << "Decomposition failed: " << to_string(result);
        throw std::runtime_error(message.str());
    }

    Eigen::VectorXf free_displacements(fea.free_dofs.size());
    {
        Scope_profiler prof_solve("Solving system");
        free_displacements = solver.solve(fea.forces);
    }
    fea.displacements.setZero();
    fea.displacements(fea.free_dofs) = free_displacements;

    if (const auto result = solver.info(); result != Eigen::Success)
    {
        std::ostringstream message;
        message << "Solving failed: " << to_string(result);
        throw std::runtime_error(message.str());
    }
}

} // namespace

FEA_state fea_init(int num_elements_x,
                   int num_elements_y,
                   float volume_fraction,
                   float penalization,
                   float radius_min,
                   float move,
                   Problem problem)
{
    FEA_state fea {};
    fea.num_elements_x = num_elements_x;
    fea.num_elements_y = num_elements_y;
    fea.num_elements = num_elements_x * num_elements_y;
    fea.num_nodes_x = num_elements_x + 1;
    fea.num_nodes_y = num_elements_y + 1;
    fea.num_nodes = fea.num_nodes_x * fea.num_nodes_y;
    fea.num_dofs_per_node = 2;
    fea.num_dofs = fea.num_nodes * fea.num_dofs_per_node;
    fea.young_modulus = 1.0f;
    fea.young_modulus_min = 1e-9f;
    fea.poisson_ratio = 0.3f;
    fea.volume_fraction = volume_fraction;
    fea.penalization = penalization;
    fea.radius_min = radius_min;
    fea.move = move;

    const Eigen::ArrayXXi node_indices {
        Eigen::ArrayXi::LinSpaced(fea.num_nodes, 0, fea.num_nodes - 1)
            .reshaped(fea.num_nodes_y, fea.num_nodes_x)};

    // Represents the index of the first DOF of each element
    const Eigen::VectorXi connectivity_vector {
        fea.num_dofs_per_node *
        node_indices.topLeftCorner(num_elements_y, num_elements_x)
            .reshaped(fea.num_elements, 1)};

    // Each row of the connectivity matrix indexes the 8 DOFs of the
    // corresponding element
    fea.connectivity_matrix =
        connectivity_vector.replicate(1, 8).rowwise() +
        Eigen::RowVector<int, 8> {2,
                                  3,
                                  fea.num_dofs_per_node * fea.num_nodes_y + 2,
                                  fea.num_dofs_per_node * fea.num_nodes_y + 3,
                                  fea.num_dofs_per_node * fea.num_nodes_y + 0,
                                  fea.num_dofs_per_node * fea.num_nodes_y + 1,
                                  0,
                                  1};

    const Eigen::Vector<int, 36> dof_connectivities_i {
        0, 1, 2, 3, 4, 5, 6, 7, 1, 2, 3, 4, 5, 6, 7, 2, 3, 4,
        5, 6, 7, 3, 4, 5, 6, 7, 4, 5, 6, 7, 5, 6, 7, 6, 7, 7};
    const Eigen::Vector<int, 36> dof_connectivities_j {
        0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2,
        2, 2, 2, 3, 3, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 6, 6, 7};
    const Eigen::VectorXi stiffness_matrix_indices_i {
        fea.connectivity_matrix(Eigen::placeholders::all, dof_connectivities_i)
            .transpose()
            .reshaped(fea.num_elements * 36, 1)};
    const Eigen::VectorXi stiffness_matrix_indices_j {
        fea.connectivity_matrix(Eigen::placeholders::all, dof_connectivities_j)
            .transpose()
            .reshaped(fea.num_elements * 36, 1)};

    fea.stiffness_matrix_indices.resize(stiffness_matrix_indices_i.rows(), 2);
    fea.stiffness_matrix_indices
        << stiffness_matrix_indices_i.cwiseMax(stiffness_matrix_indices_j),
        stiffness_matrix_indices_i.cwiseMin(stiffness_matrix_indices_j);
    fea.stiffness_matrix_values.setZero(36 * fea.num_elements);

    const Eigen::Vector<float, 36> element_stiffness_matrix_values_1 {
        12, 3,  -6, -3, -6, -3, 0, 3,  12, 3, 0,  -3, -6, -3, -6, 12, -3, 0,
        -3, -6, 3,  12, 3,  -6, 3, -6, 12, 3, -6, -3, 12, 3,  0,  12, -3, 12};
    const Eigen::Vector<float, 36> element_stiffness_matrix_values_2 {
        -4, 3, -2, 9,  2,  -3, 4, -9, -4, -9, 4,  -3, 2,  9,  -2, -4, -3, 4,
        9,  2, 3,  -4, -9, -2, 3, 2,  -4, 3,  -2, 9,  -4, -9, 4,  -4, -3, -4};
    fea.element_stiffness_matrix_values =
        1.0f / 24.0f / (1.0f - fea.poisson_ratio * fea.poisson_ratio) *
        (element_stiffness_matrix_values_1 +
         fea.poisson_ratio * element_stiffness_matrix_values_2);
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

    fea.young_moduli.setZero(fea.num_elements);

    if (problem == Problem::MBB_beam)
    {
        fea.passive_solid = {};
        fea.passive_void = {};
        Eigen::ArrayXi passive_elements(fea.passive_solid.size() +
                                        fea.passive_void.size());
        passive_elements << fea.passive_solid, fea.passive_void;
        fea.active_elements =
            filtered_index_array(fea.num_elements, passive_elements);

        Eigen::ArrayXi fixed_dofs(fea.num_nodes_y + 1);
        fixed_dofs << Eigen::ArrayXi::LinSpaced(
            fea.num_nodes_y, 0, fea.num_dofs_per_node * fea.num_nodes_y - 1),
            fea.num_dofs_per_node * node_indices(Eigen::placeholders::last,
                                                 Eigen::placeholders::last) +
                1;

        fea.free_dofs = filtered_index_array(fea.num_dofs, fixed_dofs);
    }
    else if (problem == Problem::arch)
    {
        const Eigen::ArrayXXi element_indices {
            Eigen::ArrayXi::LinSpaced(fea.num_elements, 0, fea.num_elements - 1)
                .reshaped(fea.num_elements_y, fea.num_elements_x)};

        fea.passive_solid =
            element_indices(0, Eigen::placeholders::all).reshaped();

        const auto base_length = fea.num_elements_x / 5;
        fea.passive_void =
            element_indices(
                Eigen::seq(fea.num_elements_y / 4, Eigen::placeholders::last),
                Eigen::seq(base_length,
                           Eigen::placeholders::last + 1 - base_length))
                .reshaped();

        Eigen::ArrayXi passive_elements(fea.passive_solid.size() +
                                        fea.passive_void.size());
        passive_elements << fea.passive_solid, fea.passive_void;
        std::sort(passive_elements.begin(), passive_elements.end());
        fea.active_elements =
            filtered_index_array(fea.num_elements, passive_elements);

        const auto bottom_fixed_nodes_left =
            node_indices(Eigen::placeholders::last,
                         Eigen::seq(0, base_length - 1))
                .reshaped();
        const auto bottom_fixed_nodes_right =
            node_indices(Eigen::placeholders::last,
                         Eigen::placeholders::lastN(base_length))
                .reshaped();

        Eigen::ArrayXi fixed_dofs(base_length * 4);
        fixed_dofs << fea.num_dofs_per_node * bottom_fixed_nodes_left,
            fea.num_dofs_per_node * bottom_fixed_nodes_left + 1,
            fea.num_dofs_per_node * bottom_fixed_nodes_right,
            fea.num_dofs_per_node * bottom_fixed_nodes_right + 1;
        std::sort(fixed_dofs.begin(), fixed_dofs.end());
        fea.free_dofs = filtered_index_array(fea.num_dofs, fixed_dofs);
    }
    else
    {
        assert(false);
    }

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
    if (problem == Problem::MBB_beam)
    {
        fea.forces(fea.all_to_free(fea.num_dofs_per_node * node_indices(0, 0) +
                                   1)) = -1.0f;
    }
    else if (problem == Problem::arch)
    {
        const Eigen::ArrayXi top_dofs {
            fea.num_dofs_per_node *
                node_indices(0, Eigen::placeholders::all).reshaped() +
            1};
        fea.forces(fea.all_to_free(top_dofs)).array() =
            -1.0f / static_cast<float>(fea.num_elements_x);
    }
    else
    {
        assert(false);
    }

    fea.displacements.setZero(fea.num_dofs);

    const auto kernel_size =
        2 * static_cast<Eigen::Index>(std::ceil(radius_min)) - 1;
    const float kernel_min_coord {-std::ceil(radius_min) + 1.0f};
    fea.filter_kernel = Eigen::ArrayXXf::NullaryExpr(
        kernel_size,
        kernel_size,
        [radius_min, kernel_min_coord](Eigen::Index i, Eigen::Index j)
        {
            const auto y = kernel_min_coord + static_cast<float>(i);
            const auto x = kernel_min_coord + static_cast<float>(j);
            return std::max(radius_min - std::hypot(x, y), 0.0f);
        });
    fea.filter_weights =
        filter(Eigen::ArrayXXf::Ones(num_elements_y, num_elements_x),
               fea.filter_kernel);

    fea.design_variables.setZero(fea.num_elements);
    fea.design_variables(fea.active_elements) =
        (volume_fraction *
             static_cast<float>(fea.num_elements - fea.passive_solid.size()) -
         static_cast<float>(fea.passive_solid.size())) /
        static_cast<float>(fea.active_elements.size());
    fea.design_variables(fea.passive_solid) = 1.0f;
    fea.design_variables_filtered.setZero(fea.num_elements);
    fea.design_variables_physical = fea.design_variables;
    fea.design_variables_old.setOnes(fea.num_elements);
    fea.design_variables_indexed_temp.setZero(fea.active_elements.size());

    fea.displacement_matrix.setZero(fea.connectivity_matrix.rows(),
                                    fea.connectivity_matrix.cols());

    fea.stiffness_derivative.setZero(fea.num_elements);
    fea.compliance_derivative.setZero(fea.num_elements);
    fea.filtered_compliance_derivative.setZero(fea.num_elements);
    fea.volume_derivative.setZero(fea.num_elements);
    fea.volume_derivative(fea.active_elements) =
        1.0f / static_cast<float>(fea.num_elements) / volume_fraction;
    fea.filtered_volume_derivative.setZero(fea.num_elements);

    fea.active_design_variables.setZero(fea.active_elements.size());
    fea.lower_bound.setZero(fea.active_elements.size());
    fea.upper_bound.setZero(fea.active_elements.size());
    fea.resizing_rule_constant.setZero(fea.active_elements.size());

    return fea;
}

void fea_optimization_step(FEA_state &fea)
{
    Scope_profiler prof_optimization("Optimization step");

    // TODO: design_variables could actually always have a 2D shape
    fea.design_variables_filtered =
        (filter(fea.design_variables.reshaped(fea.num_elements_y,
                                              fea.num_elements_x),
                fea.filter_kernel) /
         fea.filter_weights)
            .reshaped();

    // FIXME: even this expression requires the allocation of a temporary array.
    // I don't think there is anything we could do about it aside from making it
    // an explicit for-loop. This is probably where Eigen shows its limits when
    // an index array is necessary.
    // Note that for the current system we have to filter out inactive elements
    // because the mesh is always a rectangle (even when we explicitly disable
    // parts of it). Since we want to ultimately support arbitrary geometry,
    // this problem should go away (note that all the indexing used in here is
    // with active_elements). The indexing that will still be necessary is when
    // solving the equilibrium system, because we need to filter out fixed DOFs.
    // So, for the time being, just accept that we will have to allocate here.
    fea.design_variables_indexed_temp =
        fea.design_variables_filtered(fea.active_elements);
    fea.design_variables_physical(fea.active_elements) =
        fea.design_variables_indexed_temp;

    fea.design_variables_old = fea.design_variables_physical;

    fea.young_moduli = fea.young_modulus_min +
                       fea.design_variables_physical.pow(fea.penalization) *
                           (fea.young_modulus - fea.young_modulus_min);
    fea.stiffness_derivative(fea.active_elements) =
        -fea.penalization * (fea.young_modulus - fea.young_modulus_min) *
        fea.design_variables_physical(fea.active_elements)
            .pow(fea.penalization - 1.0f);
    fea.stiffness_matrix_values =
        (fea.element_stiffness_matrix_values * fea.young_moduli.transpose())
            .reshaped();

    solve_equilibrium_system(fea);

    fea.displacement_matrix =
        fea.displacements(fea.connectivity_matrix.reshaped())
            .reshaped(fea.connectivity_matrix.rows(),
                      fea.connectivity_matrix.cols());
    fea.compliance_derivative =
        fea.stiffness_derivative *
        (fea.displacement_matrix * fea.element_stiffness_matrix)
            .cwiseProduct(fea.displacement_matrix)
            .array()
            .rowwise()
            .sum();
    fea.filtered_compliance_derivative =
        filter(fea.compliance_derivative.reshaped(fea.num_elements_y,
                                                  fea.num_elements_x) /
                   fea.filter_weights,
               fea.filter_kernel)
            .reshaped();
    fea.filtered_volume_derivative =
        filter(fea.volume_derivative.reshaped(fea.num_elements_y,
                                              fea.num_elements_x) /
                   fea.filter_weights,
               fea.filter_kernel)
            .reshaped();

    fea.active_design_variables = fea.design_variables(fea.active_elements);

    fea.lower_bound = fea.active_design_variables - fea.move;
    fea.upper_bound = fea.active_design_variables + fea.move;
    fea.resizing_rule_constant =
        fea.active_design_variables *
        (-fea.filtered_compliance_derivative(fea.active_elements) /
         fea.filtered_volume_derivative(fea.active_elements))
            .cwiseMax(0.0f) // Drop negative coefficients for the sqrt
            .sqrt();
    // Initial estimate for LM
    float l1 {0.0f};
    float l2 {fea.resizing_rule_constant.mean() / fea.volume_fraction};
    // OC resizing rule
    while ((l2 - l1) / (l2 + l1) > 1e-4f)
    {
        const auto l_middle = 0.5f * (l1 + l2);
        fea.design_variables(fea.active_elements) =
            (fea.resizing_rule_constant / l_middle)
                .cwiseMin(fea.upper_bound)
                .cwiseMin(1.0f)
                .cwiseMax(fea.lower_bound)
                .cwiseMax(0.0f);
        if (fea.design_variables.mean() > fea.volume_fraction)
        {
            l1 = l_middle;
        }
        else
        {
            l2 = l_middle;
        }
    }
}
