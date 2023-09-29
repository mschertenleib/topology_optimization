#ifndef FEA_HPP
#define FEA_HPP

#include <Eigen/Core>

struct FEA_state
{
    int num_elements_x;
    int num_elements_y;
    int num_elements;
    int num_nodes_x;
    int num_nodes_y;
    int num_nodes;
    int num_dofs;
    int num_dofs_per_node;
    float young_modulus;
    float young_modulus_min;
    float poisson_ratio;
    float volume_fraction;
    float penalization;
    float radius_min;
    float move;
    Eigen::MatrixX2i stiffness_matrix_indices;
    Eigen::Vector<float, 36> element_stiffness_matrix_values;
    Eigen::Matrix<float, 8, 8> element_stiffness_matrix;
    Eigen::VectorXf young_moduli;
    Eigen::VectorXi passive_solid;
    Eigen::VectorXi passive_void;
    Eigen::VectorXi active_elements;
    Eigen::VectorXi free_dofs;
    Eigen::VectorXi all_to_free;
    Eigen::VectorXf forces;
    Eigen::VectorXf displacements;
    Eigen::VectorXf design_variables;
    Eigen::VectorXf design_variables_physical;
    Eigen::VectorXf design_variables_old;
    Eigen::VectorXf stiffness_derivative;
    Eigen::VectorXf volume_derivative;
    Eigen::MatrixXf filter;
};

[[nodiscard]] FEA_state fea_init(int num_elements_x,
                                 int num_elements_y,
                                 float volume_fraction,
                                 float penalization,
                                 float radius_min,
                                 float move);

void fea_solve(FEA_state &fea);

void fea_optimization_step(FEA_state &fea);

#endif // FEA_HPP
