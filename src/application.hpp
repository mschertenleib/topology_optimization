#ifndef APPLICATION_HPP
#define APPLICATION_HPP

#include "context.hpp"

class Application
{
public:
    Application();

    void run();

private:
    void update();

    Context m_context {};
};

#endif // APPLICATION_HPP
