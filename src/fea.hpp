#ifndef FEA_HPP
#define FEA_HPP

#include <Eigen/Core>
#include <Eigen/SparseCore>

struct FEA_problem
{
    int num_elements_x;
    int num_elements_y;
    int num_elements;
    int num_nodes_x;
    int num_nodes_y;
    int num_nodes;
    int num_dofs;
    Eigen::MatrixX2i stiffness_matrix_indices;
    Eigen::Vector<float, 36> element_stiffness_matrix_values;
    Eigen::VectorXf young_moduli;
    Eigen::VectorXi free_dofs;
    Eigen::VectorXf forces;
};

[[nodiscard]] FEA_problem fea_init(int num_elements_x, int num_elements_y);

[[nodiscard]] Eigen::VectorXf fea_solve(const FEA_problem &problem);

#endif // FEA_HPP
