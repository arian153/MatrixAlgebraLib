#include <iostream>
#include "MatrixLib.Math.MatrixMxN.h"
#include "MatrixLib.Math.VectorN.h"

int main()
{
    std::cout << "Matrix Lib!\n";

    MatrixLib::VectorN<7> v7_a { 2, 3, 4, 11, 12, 14, 3 };
    std::cout << v7_a << std::endl;

    MatrixLib::VectorN<7> v7_b = { 4, 5, 9, 3, 3, 2, 1 };
    std::cout << v7_b << std::endl;

    MatrixLib::VectorN<7> v7_c;
    v7_c = 4.0f;

    std::cout << v7_c << std::endl;
    v7_c = { 4, 5, 9, 3, 3, 2, 1 };
    std::cout << v7_c << std::endl;

    std::cout << v7_b - v7_a << std::endl;

    std::cout << MatrixLib::VectorN<7>::CrossProduct(v7_a, v7_b) << std::endl;

    v7_a.Normalize();
    auto v4_a = v7_a.Swizzle<4>({ 4, 4, 3 });
    std::cout << v4_a << std::endl;

    v4_a.ClearDigit(2);
    std::cout << v4_a << std::endl;

    MatrixLib::MatrixMxN<4, 3> mat_a = {
        v7_a.Swizzle<3>({ 5, 4, 3 })
        , v7_a.Swizzle<3>({ 2, 1, 6 })
        , v7_a.Swizzle<3>({ 0, 1, 2 })
        , v7_a.Swizzle<3>({ 4, 4, 5 }) };
    std::cout << "mat_a: " << mat_a << std::endl;

    mat_a.ClearDigit(3);

    std::cout << "mat_a: " << mat_a << std::endl;

    return 0;
}
