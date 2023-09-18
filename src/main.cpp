#include "application.hpp"
#include "fea.hpp"

#include <cstdlib>
#include <iostream>
#include <stdexcept>

int main()
{
    try
    {
        const auto fea_problem = fea_init(20, 10);
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
