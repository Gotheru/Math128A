# cmake_minimum_required(VERSION <specify CMake version here>)
project(Math128A)

set(CMAKE_CXX_STANDARD 14)

include_directories("Chapter 6")

add_executable(Math128A
        "Chapter 3/Hermite_Interpolation.cpp"
        "Chapter 3/Polynomial.cpp"
        "Chapter 3/temp.cpp"
        "Chapter 4/AdaptiveQuadrature.cpp"
        "Chapter 4/Romberg.cpp"
        "Chapter 4/simpsonDoubleIntegral.cpp"
        "Chapter 4/temp.cpp"
        "Chapter 5/adams.cpp"
        "Chapter 5/rungekutta.cpp"
        "Chapter 5/rungekuttafehlberg.cpp"
        "Chapter 5/rungekuttasystems.cpp"
        "Chapter 5/temp.cpp"
        "Chapter 5/trapezoidal_newton.cpp"
        "Chapter 6/gaussian_elimination_with_substitution.cpp"
)
