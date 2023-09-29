#include "application.hpp"
#include "fea.hpp"

#include <cstdlib>
#include <iostream>
#include <stdexcept>

int main()
{
    try
    {
#if 1
        auto state = fea_init(20, 10);

        constexpr float volume_fraction {0.2f};
        constexpr float penalization {3.0f};
        constexpr float radius_min {2.0f};
        constexpr float move {0.2f};
        fea_init_optimization(
            state, volume_fraction, penalization, radius_min, move);

        // fea_solve(state);
        // return application_main(state);
#else
        Eigen::SparseVector<float> sparse(10);
        sparse.insert(4) = 3.1415f;
        sparse.insert(7) = -5.37f;

        Eigen::VectorXf dense1(sparse);
        std::cout << "dense1: " << dense1.transpose() << std::endl;

        Eigen::VectorXf dense2(10);
        dense2 = sparse;
        std::cout << "dense2: " << dense2.transpose() << std::endl;

        return EXIT_SUCCESS;
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
