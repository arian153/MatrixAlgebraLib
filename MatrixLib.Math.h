#pragma once
#include <ostream>

#include "MatrixLib.Math.Constants.h"
#include "MatrixLib.Math.MatrixMxN.h"
#include "MatrixLib.Math.VectorN.h"

// global vector static methods 
namespace MatrixLib
{
    class Vector
    {
    public:
        template <SizeT N>
        const inline static VectorN<N> ZERO = VectorN<N>(Constant::ZERO);

        template <SizeT N>
        const inline static VectorN<N> X_AXIS = VectorN<N>({ Constant::ONE });

        template <SizeT N>
        const inline static VectorN<N> Y_AXIS = VectorN<N>({ Constant::ZERO, Constant::ONE });

        template <SizeT N>
        const inline static VectorN<N> Z_AXIS = VectorN<N>({ Constant::ZERO, Constant::ZERO, Constant::ONE });

        template <SizeT N>
        const inline static VectorN<N> W_AXIS = VectorN<N>({ Constant::ZERO, Constant::ZERO, Constant::ZERO, Constant::ONE });

    public:
        template <SizeT N>
        static Real DotProduct(const VectorN<N>& a, const VectorN<N>& b)
        {
            return VectorN<N>::DotProductImpl(a, b);
        }

        template <SizeT N>
        static VectorN<N> HadamardProduct(const VectorN<N>& a, const VectorN<N>& b)
        {
            return VectorN<N>::HadamardProductImpl(a, b);
        }

        static VectorN<3> CrossProduct(const VectorN<3>& a, const VectorN<3>& b)
        {
            return VectorN<3>::CrossProductImpl(a, b);
        }

        static VectorN<7> CrossProduct(const VectorN<7>& a, const VectorN<7>& b)
        {
            return VectorN<7>::CrossProductImpl(a, b);
        }

        template <SizeT Row, SizeT Col>
        static MatrixMxN<Row, Col> OuterProduct(const VectorN<Row>& vec_a, const VectorN<Col>& vec_b)
        {
            MatrixMxN outer;
            for (SizeT row = 0; row < Row; ++row)
            {
                for (SizeT col = 0; col < Col; ++col)
                {
                    outer[row][col] = vec_a[row] * vec_b[col];
                }
            }

            return outer;
        }

        template <SizeT N>
        static VectorN<N> Project(const VectorN<N>& project_a, const VectorN<N>& onto_b)
        {
            return (Vector::DotProduct(project_a, onto_b) / Vector::DotProduct(onto_b, onto_b)) * onto_b;
        }

        template <SizeT N>
        static Real Distance(const VectorN<N>& a, const VectorN<N>& b)
        {
            return (b - a).Length();
        }

        template <SizeT N>
        static Real DistanceSq(const VectorN<N>& a, const VectorN<N>& b)
        {
            return (b - a).LengthSq();
        }
    };
}

// global matrix static methods 
namespace MatrixLib
{
    class Matrix
    {
    public:
        template <SizeT N>
        const inline static MatrixMxN<N, N> IDENTITY = MatrixMxN<N, N>::Identity();

        template <SizeT Row, SizeT Col>
        const inline static MatrixMxN<Row, Col> ZERO = MatrixMxN<Row, Col>(Constant::ZERO);

    public:
        template <SizeT Row, SizeT Col>
        static MatrixMxN<Row, Col> Add(const MatrixMxN<Row, Col>& mat_a, const MatrixMxN<Row, Col>& mat_b)
        {
            MatrixMxN<Row, Col> added;
            for (SizeT row = 0; row < Row; ++row)
            {
                added[row] = mat_a[row] + mat_b[row];
            }

            return added;
        }

        template <SizeT Row, SizeT Col>
        static MatrixMxN<Row, Col> Subtract(const MatrixMxN<Row, Col>& mat_a, const MatrixMxN<Row, Col>& mat_b)
        {
            MatrixMxN<Row, Col> subtracted;
            for (SizeT row = 0; row < Row; ++row)
            {
                subtracted[row] = mat_a[row] - mat_b[row];
            }

            return subtracted;
        }

        template <SizeT Row, SizeT Col>
        static MatrixMxN<Row, Col> Multiply(const MatrixMxN<Row, Col>& mat, Real r)
        {
            MatrixMxN<Row, Col> multiplied;
            for (SizeT row = 0; row < Row; ++row)
            {
                multiplied[row] = mat[row] * r;
            }

            return multiplied;
        }

        template <SizeT Row, SizeT Col>
        static MatrixMxN<Row, Col> Multiply(Real r, const MatrixMxN<Row, Col>& mat)
        {
            MatrixMxN<Row, Col> multiplied;
            for (SizeT row = 0; row < Row; ++row)
            {
                multiplied[row] = mat[row] * r;
            }

            return multiplied;
        }

        template <SizeT Row, SizeT Col>
        static VectorN<Row> Multiply(const MatrixMxN<Row, Col>& mat, const VectorN<Col>& vec)
        {
            VectorN<Row> multiplied;
            for (SizeT row = 0; row < Row; ++row)
            {
                multiplied[row] = Vector::DotProduct(mat[row], vec);
            }

            return multiplied;
        }

        template <SizeT Row, SizeT Col>
        static VectorN<Col> Multiply(const VectorN<Row>& vec, const MatrixMxN<Row, Col>& mat)
        {
            VectorN<Col> multiplied;
            for (SizeT col = 0; col < Col; ++col)
            {
                multiplied[col] = Vector::DotProduct(vec, mat.Column(col));
            }

            return multiplied;
        }

        template <SizeT RowA, SizeT ColRow, SizeT ColB>
        static MatrixMxN<RowA, ColB> Multiply(MatrixMxN<RowA, ColRow> mat_a, MatrixMxN<ColRow, ColB> mat_b)
        {
            MatrixMxN<RowA, ColB> multiplied;
            for (SizeT row = 0; row < RowA; ++row)
            {
                for (SizeT col = 0; col < ColB; ++col)
                {
                    multiplied(row, col) = Vector::DotProduct(mat_a.Row(row), mat_b.Column(col));
                }
            }

            return multiplied;
        }

        template <SizeT Row, SizeT Col>
        static MatrixMxN<Row, Col> Divide(const MatrixMxN<Row, Col>& mat, Real r)
        {
            MatrixMxN<Row, Col> divided;
            for (SizeT row = 0; row < Row; ++row)
            {
                divided[row] = mat[row] * Tools::Inverse(r);
            }

            return divided;
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
            // Solve A * x = b;
            // x = A^-1 * b;

            MatrixMxN<ColA, Row>  inv_a = Matrix::PseudoInverse(mat_a);
            MatrixMxN<ColA, ColB> x     = Matrix::Multiply(inv_a, mat_b);

            return x;
        }

        template <SizeT Row, SizeT Col>
        static MatrixMxN<Row, Col> HadamardProduct(const MatrixMxN<Row, Col>& mat_a, const MatrixMxN<Row, Col>& mat_b)
        {
            MatrixMxN hadamard;
            for (SizeT row = 0; row < Row; ++row)
            {
                for (SizeT col = 0; col < Col; ++col)
                {
                    hadamard[row][col] = mat_a[row][col] * mat_b[row][col];
                }
            }

            return hadamard;
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
}

// global arithmetic operators
namespace MatrixLib
{
    template <SizeT N>
    VectorN<N> operator +(const VectorN<N>& a, const VectorN<N>& b)
    {
        VectorN<N> added;
        for (SizeT i = 0; i < N; ++i)
        {
            added[i] = a[i] + b[i];
        }

        return added;
    }

    template <SizeT N>
    VectorN<N> operator -(const VectorN<N>& a, const VectorN<N>& b)
    {
        VectorN<N> subtracted;
        for (SizeT i = 0; i < N; ++i)
        {
            subtracted[i] = a[i] - b[i];
        }

        return subtracted;
    }

    template <SizeT N>
    VectorN<N> operator *(Real real, const VectorN<N>& vector)
    {
        VectorN<N> multiplied;
        for (SizeT i = 0; i < N; ++i)
        {
            multiplied[i] = vector[i] * real;
        }

        return multiplied;
    }

    template <SizeT N>
    VectorN<N> operator *(const VectorN<N>& vector, Real real)
    {
        VectorN<N> multiplied;
        for (SizeT i = 0; i < N; ++i)
        {
            multiplied[i] = vector[i] * real;
        }

        return multiplied;
    }

    template <SizeT N>
    VectorN<N> operator /(const VectorN<N>& vector, Real real)
    {
        VectorN<N> multiplied;
        for (SizeT i = 0; i < N; ++i)
        {
            multiplied[i] = vector[i] * Tools::Inverse(real);
        }

        return multiplied;
    }

    template <SizeT Row, SizeT Col>
    MatrixMxN<Row, Col> operator +(const MatrixMxN<Row, Col>& a, const MatrixMxN<Row, Col>& b)
    {
        MatrixMxN<Row, Col> added;
        for (SizeT row = 0; row < Row; ++row)
        {
            added[row] = a[row] + b[row];
        }

        return added;
    }

    template <SizeT Row, SizeT Col>
    MatrixMxN<Row, Col> operator -(const MatrixMxN<Row, Col>& a, const MatrixMxN<Row, Col>& b)
    {
        MatrixMxN<Row, Col> subtracted;
        for (SizeT row = 0; row < Row; ++row)
        {
            subtracted[row] = a[row] - b[row];
        }

        return subtracted;
    }

    template <SizeT Row, SizeT Col>
    static MatrixMxN<Row, Col> operator *(const MatrixMxN<Row, Col>& mat, Real r)
    {
        MatrixMxN<Row, Col> multiplied;
        for (SizeT row = 0; row < Row; ++row)
        {
            multiplied[row] = mat[row] * r;
        }

        return multiplied;
    }

    template <SizeT Row, SizeT Col>
    static MatrixMxN<Row, Col> operator *(Real r, const MatrixMxN<Row, Col>& mat)
    {
        MatrixMxN<Row, Col> multiplied;
        for (SizeT row = 0; row < Row; ++row)
        {
            multiplied[row] = mat[row] * r;
        }

        return multiplied;
    }

    template <SizeT Row, SizeT Col>
    static VectorN<Row> operator *(const MatrixMxN<Row, Col>& mat, const VectorN<Col>& vec)
    {
        VectorN<Row> multiplied;
        for (SizeT row = 0; row < Row; ++row)
        {
            multiplied[row] = Vector::DotProduct(mat[row], vec);
        }

        return multiplied;
    }

    template <SizeT Row, SizeT Col>
    static VectorN<Col> operator *(const VectorN<Row>& vec, const MatrixMxN<Row, Col>& mat)
    {
        VectorN<Col> multiplied;
        for (SizeT col = 0; col < Col; ++col)
        {
            multiplied[col] = Vector::DotProduct(vec, mat.Column(col));
        }

        return multiplied;
    }

    template <SizeT RowA, SizeT ColRow, SizeT ColB>
    static MatrixMxN<RowA, ColB> operator *(MatrixMxN<RowA, ColRow> mat_a, MatrixMxN<ColRow, ColB> mat_b)
    {
        MatrixMxN<RowA, ColB> multiplied;
        for (SizeT row = 0; row < RowA; ++row)
        {
            for (SizeT col = 0; col < ColB; ++col)
            {
                multiplied(row, col) = Vector::DotProduct(mat_a.Row(row), mat_b.Column(col));
            }
        }

        return multiplied;
    }

    template <SizeT Row, SizeT Col>
    static MatrixMxN<Row, Col> operator /(const MatrixMxN<Row, Col>& mat, Real r)
    {
        MatrixMxN<Row, Col> divided;
        for (SizeT row = 0; row < Row; ++row)
        {
            divided[row] = mat[row] * Tools::Inverse(r);
        }

        return divided;
    }
}

// global IO operators
namespace MatrixLib
{
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
            os << ",\n  " << rhs[i];
        }
        os << " }";
        return os;
    }
}
