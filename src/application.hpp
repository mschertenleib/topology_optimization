#ifndef APPLICATION_HPP
#define APPLICATION_HPP

#include "fea.hpp"

int application_main(const FEA_problem &fea_problem,
                     const Eigen::VectorXf &displacements);

#endif // APPLICATION_HPP
