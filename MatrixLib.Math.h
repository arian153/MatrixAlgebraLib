#pragma once
#include <ostream>

#include "MatrixLib.Math.Constants.h"
#include "MatrixLib.Math.MatrixMxN.h"
#include "MatrixLib.Math.VectorN.h"

namespace MatrixLib
{
    class Matrix
    {
    public:
        // Add
        // Subtract
        // Multiply
        // Divide
        // Negate
        // Clear

        template <SizeT RowA, SizeT ColRow, SizeT ColB>
        static MatrixMxN<RowA, ColB> Multiply(MatrixMxN<RowA, ColRow> mat_a, MatrixMxN<ColRow, ColB> mat_b)
        {
            return MatrixMxN<RowA, ColRow>::MultiplyImpl(mat_a, mat_b);
        }

        template <SizeT Row, SizeT Col>
        static MatrixMxN<Col, Row> Transpose(const MatrixMxN<Row, Col>& mat)
        {
            return MatrixMxN<Row, Col>::TransposeImpl(mat);
        }

        template <SizeT Row, SizeT Col>
        static SizeT Rank(const MatrixMxN<Row, Col>& mat, Real tolerance = Constant::EPSILON)
        {
            return MatrixMxN<Row, Col>::RankImpl(mat, tolerance);
        }

        template <SizeT Row, SizeT Col>
        static Real Determinant(const MatrixMxN<Row, Col>& mat)
        {
            return MatrixMxN<Row, Col>::DeterminantImpl(mat);
        }

        template <SizeT Row, SizeT Col>
        static MatrixMxN<Row, Col> Inverse(const MatrixMxN<Row, Col>& mat, bool use_permute = true)
        {
            return MatrixMxN<Row, Col>::InverseImpl(mat, use_permute);
        }

        template <SizeT Row, SizeT Col>
        static MatrixMxN<Col, Row> PseudoInverse(const MatrixMxN<Row, Col>& mat, Real tolerance = Constant::EPSILON)
        {
            return MatrixMxN<Row, Col>::PseudoInverseImpl(mat, tolerance);
        }

        template <SizeT Row, SizeT ColA, SizeT ColB>
        static MatrixMxN<ColA, ColB> Solve(const MatrixMxN<Row, ColA>& mat_a, const MatrixMxN<Row, ColB>& mat_b)
        {
            return MatrixMxN<Row, ColA>::SolveImpl(mat_a, mat_b);
        }

        template <SizeT Row, SizeT Col>
        static MatrixMxN<Row, Col> HadamardProduct(const MatrixMxN<Row, Col>& mat_a, const MatrixMxN<Row, Col>& mat_b)
        {
            return MatrixMxN<Row, Col>::HadamardProductImpl(mat_a, mat_b);
        }

        template <SizeT Row, SizeT Col>
        static MatrixMxN<Row, Col> OuterProduct(const VectorN<Row>& vec_a, const VectorN<Col>& vec_b)
        {
            return MatrixMxN<Row, Col>::OuterProductImpl(vec_a, vec_b);
        }

        template <SizeT Row, SizeT Col>
        static MatrixMxN<Row - 1, Col - 1> MinorMatrix(const MatrixMxN<Row, Col>& mat, SizeT row, SizeT col)
        {
            return MatrixMxN<Row, Col>::MinorMatrixImpl(mat, row, col);
        }

        template <SizeT N>
        static MatrixMxN<N, N> AdjugateMatrix(const MatrixMxN<N, N>& mat)
        {
            // https://en.wikipedia.org/wiki/Adjugate_matrix
            Real det      = Matrix::Determinant(mat);
            auto adjugate = Matrix::Inverse(mat, true);
            adjugate *= det;

            return adjugate;
        }

        template <SizeT N>
        static MatrixMxN<N, N> CofactorMatrix(const MatrixMxN<N, N>& mat)
        {
            // https://en.wikipedia.org/wiki/Minor_(linear_algebra)
            Real det      = Matrix::Determinant(mat);
            auto cofactor = Matrix::Transpose(Matrix::Inverse(mat, true));
            cofactor *= det;

            return cofactor;
        }
    };

    class Vector
    {
    public:
        //template <SizeT N>

        // Add
        // Subtract
        // Multiply
        // Divide
        // Negate
        // Project
        // Distance
        // DistanceSq
        // DotProduct
        // HadamardProduct
        // CrossProduct
        // Norm
    };

    // global arithmetic operators
    template <SizeT N>
    VectorN<N> operator +(const VectorN<N>& a, const VectorN<N>& b)
    {
        return a.Added(b);
    }

    template <SizeT N>
    VectorN<N> operator -(const VectorN<N>& a, const VectorN<N>& b)
    {
        return a.Subtracted(b);
    }

    template <SizeT N>
    VectorN<N> operator *(Real real, const VectorN<N>& vector)
    {
        return vector.Multiplied(real);
    }

    template <SizeT N>
    VectorN<N> operator *(const VectorN<N>& vector, Real real)
    {
        return vector.Multiplied(real);
    }

    template <SizeT N>
    VectorN<N> operator /(const VectorN<N>& vector, Real real)
    {
        return vector.Divided(real);
    }

    // global IO operators
    template <SizeT N>
    std::ostream& operator<<(std::ostream& os, const VectorN<N>& rhs)
    {
        os << "{" << rhs[0];
        for (SizeT i = 1; i < N; ++i)
        {
            os << ", " << rhs[i];
        }
        os << "}";

        return os;
    }

    template <SizeT M, SizeT N>
    std::ostream& operator<<(std::ostream& os, const MatrixMxN<M, N>& rhs)
    {
        os << "{ " << rhs[0];
        for (SizeT i = 1; i < M; ++i)
        {
            os << ",\n " << rhs[i];
        }
        os << " }";
        return os;
    }
}
