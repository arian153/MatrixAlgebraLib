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

    //v7_a.Normalize();
    auto v4_a = v7_a.Swizzle<4>({ 4, 4, 3 });
    std::cout << v4_a << std::endl;

    v4_a.ClearDigit(2);
    std::cout << v4_a << std::endl;

    MatrixLib::MatrixMxN<3, 3> mat_a = {
        v7_a.Swizzle<3>({ 5, 4, 3 }),
        v7_a.Swizzle<3>({ 2, 1, 6 }),
        v7_a.Swizzle<3>({ 0, 1, 2 })
    };
    std::cout << "mat_a: " << mat_a << std::endl;

    //mat_a.ClearDigit(3);
    std::cout << "mat_a: " << mat_a << std::endl;

    auto inv_a = MatrixLib::Matrix::Inverse(mat_a);
    std::cout << "inv_a: " << inv_a << std::endl;

    auto mul_a = MatrixLib::Matrix::Multiply(mat_a, inv_a);
    mul_a.ClearDigit(5);
    std::cout << "mul_a: " << mul_a << std::endl;

    auto mul_b = MatrixLib::Matrix::Multiply(inv_a, mat_a);
    mul_b.ClearDigit(5);
    std::cout << "mul_b: " << mul_b << std::endl;

    auto inv_aa = MatrixLib::Matrix::Inverse(inv_a);
    std::cout << "inv_aa: " << inv_aa << std::endl;
    auto mul_c = MatrixLib::Matrix::Multiply(inv_a, inv_aa);
    mul_c.ClearDigit(5);
    std::cout << "mul_c: " << mul_c << std::endl;

    MatrixLib::MatrixMxN<3, 5> mat_d = {
        v7_a.Swizzle<5>({ 5, 4, 3, 3, 1 }),
        v7_a.Swizzle<5>({ 2, 1, 6, 0, 3 }),
        v7_a.Swizzle<5>({ 0, 1, 2, 4, 5 })
    };

    std::cout << "mat_d: " << mat_d << std::endl;
    std::cout << "trs_d: " << MatrixLib::Matrix::Transpose(mat_d) << std::endl;

    auto inv_d = MatrixLib::Matrix::PseudoInverse(mat_d);
    std::cout << "inv_d: " << inv_d << std::endl;

    auto mul_d = MatrixLib::Matrix::Multiply(mat_d, inv_d);
    mul_d.ClearDigit(5);
    std::cout << "mul_d 1: " << mul_d << std::endl;

    MatrixLib::MatrixMxN<3, 1> mat_e = { 12.0f, 43.0f, 23.0f };

    auto x = MatrixLib::Matrix::Solve(mat_d, mat_e);

    std::cout << "Ax = b" << std::endl;
    std::cout << "A:" << mat_d << std::endl;

    std::cout << "b: " << mat_e << std::endl;

    std::cout << "x: " << x << std::endl;

    auto b = MatrixLib::Matrix::Multiply(mat_d, x);
    b.ClearDigit(5);
    std::cout << "b: " << b << std::endl;

    return 0;
}
