#include <iostream>
#include "MatrixLib.Math.VectorN.h"

int main()
{
    std::cout << "Matrix Lib!\n";

    MatrixLib::VectorN vec(9, 1.0f);
    std::cout << vec << std::endl;

    return 0;
}
