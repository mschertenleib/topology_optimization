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
    int num_dofs_per_node;
    int num_dofs;
    float young_modulus;
    float young_modulus_min;
    float poisson_ratio;
    float volume_fraction;
    float penalization;
    float radius_min;
    float move;
    Eigen::Array<int, Eigen::Dynamic, 8> connectivity_matrix;
    Eigen::Vector<float, 36> element_stiffness_matrix_values;
    Eigen::Matrix<float, 8, 8> element_stiffness_matrix;
    Eigen::ArrayX2i stiffness_matrix_indices;
    Eigen::ArrayXf stiffness_matrix_values;
    Eigen::VectorXf young_moduli;
    Eigen::VectorXf forces;
    Eigen::VectorXf displacements;
    Eigen::ArrayXi passive_solid;
    Eigen::ArrayXi passive_void;
    Eigen::ArrayXi active_elements;
    Eigen::ArrayXi free_dofs;
    Eigen::ArrayXi all_to_free;
    Eigen::ArrayXf design_variables;
    Eigen::ArrayXf design_variables_filtered;
    Eigen::ArrayXf design_variables_physical;
    Eigen::ArrayXf design_variables_old;
    Eigen::ArrayXf design_variables_indexed_temp;
    Eigen::MatrixXf displacement_matrix;
    Eigen::ArrayXf stiffness_derivative;
    Eigen::ArrayXf compliance_derivative;
    Eigen::ArrayXf filtered_compliance_derivative;
    Eigen::ArrayXf volume_derivative;
    Eigen::ArrayXf filtered_volume_derivative;
    Eigen::ArrayXf active_design_variables;
    Eigen::ArrayXXf filter_weights;
    Eigen::ArrayXXf filter_kernel;
    Eigen::ArrayXf lower_bound;
    Eigen::ArrayXf upper_bound;
    Eigen::ArrayXf resizing_rule_constant;
};

enum struct Problem
{
    MBB_beam,
    arch
};

[[nodiscard]] FEA_state fea_init(int num_elements_x,
                                 int num_elements_y,
                                 float volume_fraction,
                                 float penalization,
                                 float radius_min,
                                 float move,
                                 Problem problem);

void fea_optimization_step(FEA_state &fea);

#endif // FEA_HPP
