#include "application.hpp"

#include <cstdlib>
#include <iostream>
#include <stdexcept>

#include <Eigen/Dense>

int main()
{
    try
    {
        Eigen::MatrixXd m(2, 2);
        m(0, 0) = 3;
        m(1, 0) = 2.5;
        m(0, 1) = -1;
        m(1, 1) = m(1, 0) + m(0, 1);
        std::cout << m << std::endl;

        Application app;
        app.run();
        return EXIT_SUCCESS;
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
