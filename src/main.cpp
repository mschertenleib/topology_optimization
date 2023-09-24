#include "application.hpp"
#include "fea.hpp"

#include <cstdlib>
#include <iostream>
#include <stdexcept>

int main()
{
    try
    {
        auto state = fea_init(20, 10);
        fea_solve(state);
#if 0
        return EXIT_SUCCESS;
#else
        return application_main(state);
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
