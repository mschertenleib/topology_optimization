#include "application.hpp"
#include "fea.hpp"

#include <cstdlib>
#include <iostream>
#include <stdexcept>

int main()
{
    try
    {
        constexpr float volume_fraction {0.3f};
        constexpr float penalization {3.0f};
        constexpr float radius_min {2.0f};
        constexpr float move {0.2f};

        auto fea =
            fea_init(60, 30, volume_fraction, penalization, radius_min, move);

        return application_main(fea);
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
