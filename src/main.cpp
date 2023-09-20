#include "application.hpp"
#include "fea.hpp"

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <stdexcept>

namespace
{

template <typename T>
    requires std::is_fundamental_v<T>
void store(std::ofstream &file, const T &value)
{
    file.write(reinterpret_cast<const char *>(&value),
               static_cast<std::streamsize>(sizeof(value)));
    assert(file);
}

template <typename T>
    requires std::is_fundamental_v<T>
void load(std::ifstream &file, T &value)
{
    file.read(reinterpret_cast<char *>(&value),
              static_cast<std::streamsize>(sizeof(value)));
    assert(file);
}

template <typename Scalar,
          int Rows,
          int Cols,
          int Options,
          int MaxRows,
          int MaxCols>
void store(
    std::ofstream &file,
    const Eigen::Matrix<Scalar, Rows, Cols, Options, MaxRows, MaxCols> &m)
{
    store(file, m.rows());
    store(file, m.cols());
    file.write(reinterpret_cast<const char *>(m.data()),
               static_cast<std::streamsize>(m.size()) *
                   static_cast<std::streamsize>(sizeof(Scalar)));
    assert(file);
}

template <typename Scalar,
          int Rows,
          int Cols,
          int Options,
          int MaxRows,
          int MaxCols>
void load(std::ifstream &file,
          Eigen::Matrix<Scalar, Rows, Cols, Options, MaxRows, MaxCols> &m)
{
    Eigen::Index rows {};
    load(file, rows);
    Eigen::Index cols {};
    load(file, cols);
    m.resize(rows, cols);
    file.read(reinterpret_cast<char *>(m.data()),
              static_cast<std::streamsize>(m.size()) *
                  static_cast<std::streamsize>(sizeof(Scalar)));
    assert(file);
}

void store(const char *filename, const FEA_problem &problem)
{
    std::ofstream file(filename, std::ios::binary);
    assert(file);

    store(file, problem.num_elements_x);
    store(file, problem.num_elements_y);
    store(file, problem.num_elements);
    store(file, problem.num_nodes_x);
    store(file, problem.num_nodes_y);
    store(file, problem.num_nodes);
    store(file, problem.num_dofs);
    store(file, problem.stiffness_matrix_indices);
    store(file, problem.element_stiffness_matrix_values);
    store(file, problem.young_moduli);
    store(file, problem.free_dofs);
    store(file, problem.forces);
}

void load(const char *filename, FEA_problem &problem)
{
    std::ifstream file(filename, std::ios::binary);
    assert(file);

    load(file, problem.num_elements_x);
    load(file, problem.num_elements_y);
    load(file, problem.num_elements);
    load(file, problem.num_nodes_x);
    load(file, problem.num_nodes_y);
    load(file, problem.num_nodes);
    load(file, problem.num_dofs);
    load(file, problem.stiffness_matrix_indices);
    load(file, problem.element_stiffness_matrix_values);
    load(file, problem.young_moduli);
    load(file, problem.free_dofs);
    load(file, problem.forces);
}

void check_equal(const FEA_problem &a, const FEA_problem &b)
{
    assert(a.num_elements_x == b.num_elements_x);
    assert(a.num_elements_y == b.num_elements_y);
    assert(a.num_elements == b.num_elements);
    assert(a.num_nodes_x == b.num_nodes_x);
    assert(a.num_nodes_y == b.num_nodes_y);
    assert(a.num_nodes == b.num_nodes);
    assert(a.num_dofs == b.num_dofs);
    assert(a.stiffness_matrix_indices == b.stiffness_matrix_indices);
    assert(a.element_stiffness_matrix_values ==
           b.element_stiffness_matrix_values);
    assert(a.young_moduli == b.young_moduli);
    assert(a.free_dofs == b.free_dofs);
    assert(a.forces == b.forces);
}

} // namespace

int main()
{
    try
    {
        const auto fea_problem = fea_init(20, 10);

        /*store("dump.bin", fea_problem);
        FEA_problem fea_problem_2 {};
        load("dump.bin", fea_problem_2);
        check_equal(fea_problem, fea_problem_2);*/

        const auto displacements = fea_solve(fea_problem);
#if 0
        return EXIT_SUCCESS;
#else
        return application_main(fea_problem, displacements);
#endif
    }
    catch (const std::exception &e)
    {
        std::cerr << "Exception thrown: " << e.what() << '\n';
        return EXIT_FAILURE;
    }
    catch (...)
    {
        std::cerr << "Unknown exception thrown\n";
        return EXIT_FAILURE;
    }
}
