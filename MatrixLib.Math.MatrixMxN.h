#pragma once
#include <array>
#include <ostream>
#include <ranges>
#include <vector>

#include "MatrixLib.Math.Constants.h"
#include "MatrixLib.Math.Tools.h"
#include "MatrixLib.Math.VectorN.h"

namespace MatrixLib
{
    template <SizeT M, SizeT N>
    class MatrixMxN
    {
    private: // private data
        std::array<VectorN<N>, M> m_elements;

    public: // constructors
        MatrixMxN();
        MatrixMxN(std::initializer_list<Real> values);
        MatrixMxN(std::initializer_list<std::initializer_list<Real>> values);
        MatrixMxN(std::initializer_list<VectorN<N>> values);
        MatrixMxN(const MatrixMxN& rhs);
        ~MatrixMxN() = default;
        explicit MatrixMxN(Real scalar);
        explicit MatrixMxN(const std::vector<Real>& vec);
        explicit MatrixMxN(const std::array<Real, M * N>& arr);
        explicit MatrixMxN(const std::vector<VectorN<N>>& vec);
        explicit MatrixMxN(const std::array<VectorN<N>, M>& arr);

    public: // operators
        // explicit operator
        explicit operator std::vector<Real>() const;

        // access operator
        Real        operator[](SizeT i) const;
        Real&       operator[](SizeT i);
        Real        operator()(SizeT row, SizeT col) const;
        Real&       operator()(SizeT row, SizeT col);
        VectorN<N>  operator()(SizeT row) const;
        VectorN<N>& operator()(SizeT row);

        //// assignment operator
        //MatrixMxN& operator =(const MatrixMxN& rhs);

        //// compare operator
        //bool operator ==(const MatrixMxN& rhs) const;
        //bool operator !=(const MatrixMxN& rhs) const;

        //// arithmetic operator
        //MatrixMxN  operator -() const;
        //MatrixMxN& operator +=(const MatrixMxN& rhs);
        //MatrixMxN& operator -=(const MatrixMxN& rhs);
        //MatrixMxN& operator *=(Real r);
        //MatrixMxN& operator /=(Real r);

    public: // modify methods
        void Set(const MatrixMxN& mat);
        void SetRow(SizeT row, VectorN<N> row_vector);
        void SetColumn(SizeT col, VectorN<M> column_vector);
        void Add(const MatrixMxN& v);
        void Subtract(const MatrixMxN& v);
        void Multiply(Real r);
        void Divide(Real r);
        void Negate();
        void Inverse();
        //void Transpose();
        void ClearDigit(SizeT digit);
        void ClearError(Real epsilon);

    public:
        static bool       LUPDecompose(const MatrixMxN& a, MatrixMxN& lu, std::vector<SizeT>& p, Real tolerance = Constant::EPSILON);
        static VectorN<N> LUPSolve(const MatrixMxN& lu, const std::vector<SizeT>& p, const VectorN<N>& b);
        static MatrixMxN  LUPInverse(const MatrixMxN& lu, const std::vector<SizeT>& p);
        static Real       LUPDeterminant(const MatrixMxN& lu, const std::vector<SizeT>& p);

    public: // range based for loop related methods
        auto begin();
        auto end();
        auto begin() const;
        auto end() const;
    };

    template <SizeT M, SizeT N>
    MatrixMxN<M, N>::MatrixMxN()
    {
        m_elements.fill(VectorN<N>(Constant::ZERO));
    }

    template <SizeT M, SizeT N>
    MatrixMxN<M, N>::MatrixMxN(std::initializer_list<Real> values)
    {
        m_elements.fill(VectorN<N>(Constant::ZERO));
        SizeT count = Tools::Min(M * N, values.size());
        for (SizeT i = 0; i < count; ++i)
        {
            SizeT row            = i / N;
            SizeT col            = i % N;
            m_elements[row][col] = *(values.begin() + i);
        }
    }

    template <SizeT M, SizeT N>
    MatrixMxN<M, N>::MatrixMxN(std::initializer_list<std::initializer_list<Real>> values)
    {
        m_elements.fill(VectorN<N>(Constant::ZERO));

        SizeT row_count = Tools::Min(M, values.size());
        for (SizeT row = 0; row < row_count; ++row)
        {
            auto& elements     = *(values.begin() + row);
            SizeT column_count = Tools::Min(N, elements.size());
            for (SizeT col = 0; col < column_count; ++col)
            {
                m_elements[row][col] = *(elements.begin() + col);
            }
        }
    }

    template <SizeT M, SizeT N>
    MatrixMxN<M, N>::MatrixMxN(std::initializer_list<VectorN<N>> values)
    {
        m_elements.fill(VectorN<N>(Constant::ZERO));
        SizeT row_count = Tools::Min(M, values.size());
        for (SizeT row = 0; row < row_count; ++row)
        {
            auto& elements = *(values.begin() + row);
            for (SizeT col = 0; col < N; ++col)
            {
                m_elements[row][col] = elements[col];
            }
        }
    }

    template <SizeT M, SizeT N>
    MatrixMxN<M, N>::MatrixMxN(const MatrixMxN& rhs)
        : m_elements(rhs.m_elements)
    {
    }

    template <SizeT M, SizeT N>
    MatrixMxN<M, N>::MatrixMxN(Real scalar)
    {
        m_elements.fill(VectorN<N>(scalar));
    }

    template <SizeT M, SizeT N>
    MatrixMxN<M, N>::MatrixMxN(const std::vector<Real>& vec)
    {
        m_elements.fill(VectorN<N>(Constant::ZERO));
        SizeT count = Tools::Min(M * N, vec.size());
        for (SizeT i = 0; i < count; ++i)
        {
            SizeT row            = i / N;
            SizeT col            = i % N;
            m_elements[row][col] = *(vec.begin() + i);
        }
    }

    template <SizeT M, SizeT N>
    MatrixMxN<M, N>::MatrixMxN(const std::array<Real, M * N>& arr)
    {
        m_elements.fill(VectorN<N>(Constant::ZERO));
        for (SizeT i = 0; i < M * N; ++i)
        {
            SizeT row            = i / N;
            SizeT col            = i % N;
            m_elements[row][col] = arr[i];
        }
    }

    template <SizeT M, SizeT N>
    MatrixMxN<M, N>::MatrixMxN(const std::vector<VectorN<N>>& vec)
    {
        m_elements.fill(VectorN<N>(Constant::ZERO));
        SizeT row_count = Tools::Min(M, vec.size());
        for (SizeT row = 0; row < row_count; ++row)
        {
            auto& elements = vec[row];
            for (SizeT col = 0; col < N; ++col)
            {
                m_elements[row][col] = elements[col];
            }
        }
    }

    template <SizeT M, SizeT N>
    MatrixMxN<M, N>::MatrixMxN(const std::array<VectorN<N>, M>& arr)
        : m_elements(arr)
    {
    }

    template <SizeT M, SizeT N>
    MatrixMxN<M, N>::operator std::vector<Real>() const
    {
        std::vector<Real> result;
        result.reserve(M * N);
        for (SizeT i = 0; i < M; ++i)
        {
            result.insert(result.end(), m_elements[i].begin(), m_elements[i].end());
        }

        return result;
    }

    template <SizeT M, SizeT N>
    Real MatrixMxN<M, N>::operator[](SizeT i) const
    {
        SizeT row    = i / N;
        SizeT column = i % N;
        return m_elements[row][column];
    }

    template <SizeT M, SizeT N>
    Real& MatrixMxN<M, N>::operator[](SizeT i)
    {
        SizeT row    = i / N;
        SizeT column = i % N;
        return m_elements[row][column];
    }

    template <SizeT M, SizeT N>
    Real MatrixMxN<M, N>::operator()(SizeT row, SizeT col) const
    {
        return m_elements[row % M][col % N];
    }

    template <SizeT M, SizeT N>
    Real& MatrixMxN<M, N>::operator()(SizeT row, SizeT col)
    {
        return m_elements[row % M][col % N];
    }

    template <SizeT M, SizeT N>
    VectorN<N> MatrixMxN<M, N>::operator()(SizeT row) const
    {
        return m_elements[row % M];
    }

    template <SizeT M, SizeT N>
    VectorN<N>& MatrixMxN<M, N>::operator()(SizeT row)
    {
        return m_elements[row % M];
    }

    template <SizeT M, SizeT N>
    void MatrixMxN<M, N>::Set(const MatrixMxN& mat)
    {
        m_elements = mat.m_elements;
    }

    template <SizeT M, SizeT N>
    void MatrixMxN<M, N>::SetRow(SizeT row, VectorN<N> row_vector)
    {
        m_elements[row] = row_vector;
    }

    template <SizeT M, SizeT N>
    void MatrixMxN<M, N>::SetColumn(SizeT col, VectorN<M> column_vector)
    {
        for (SizeT row = 0; row < M; ++row)
        {
            m_elements[row][col] = column_vector[col];
        }
    }

    template <SizeT M, SizeT N>
    void MatrixMxN<M, N>::Add(const MatrixMxN& v)
    {
        for (SizeT row = 0; row < M; ++row)
        {
            m_elements[row] += v(row);
        }
    }

    template <SizeT M, SizeT N>
    void MatrixMxN<M, N>::Subtract(const MatrixMxN& v)
    {
        for (SizeT row = 0; row < M; ++row)
        {
            m_elements[row] -= v(row);
        }
    }

    template <SizeT M, SizeT N>
    void MatrixMxN<M, N>::Multiply(Real r)
    {
        for (SizeT row = 0; row < M; ++row)
        {
            m_elements[row] *= r;
        }
    }

    template <SizeT M, SizeT N>
    void MatrixMxN<M, N>::Divide(Real r)
    {
        for (SizeT row = 0; row < M; ++row)
        {
            m_elements[row] /= r;
        }
    }

    template <SizeT M, SizeT N>
    void MatrixMxN<M, N>::Negate()
    {
        for (SizeT row = 0; row < M; ++row)
        {
            m_elements[row].Negate();
        }
    }

    template <SizeT M, SizeT N>
    void MatrixMxN<M, N>::Inverse()
    {
        if (M == N)
        {
            MatrixMxN          lu;
            std::vector<SizeT> p;

            if (LUPDecompose(this, lu, p))
            {
                MatrixMxN invert = LUPInverse(lu, p);
                Set(invert);
            }
        }
    }

    template <SizeT M, SizeT N>
    void MatrixMxN<M, N>::ClearDigit(SizeT digit)
    {
        for (SizeT row = 0; row < M; ++row)
        {
            m_elements[row].ClearDigit(digit);
        }
    }

    template <SizeT M, SizeT N>
    void MatrixMxN<M, N>::ClearError(Real epsilon)
    {
        for (SizeT row = 0; row < M; ++row)
        {
            m_elements[row].ClearError(epsilon);
        }
    }

    template <SizeT M, SizeT N>
    bool MatrixMxN<M, N>::LUPDecompose(const MatrixMxN& a, MatrixMxN& lu, std::vector<SizeT>& p, Real tolerance)
    {
        // only decomposable NxN Matrix
        if (M != N)
            return false;

        lu.Set(a);
        p.resize(N + 1);

        // Unit permutation matrix, P[N] initialized with N
        for (SizeT i = 0; i <= N; ++i)
        {
            p[i] = i;
        }

        for (SizeT i = 0; i < N; ++i)
        {
            Real  max_a = 0.0;
            SizeT max_i = i;

            for (SizeT k = i; k < N; k++)
            {
                Real abs_a = Tools::Abs(a[k][i]);
                if (abs_a > max_a)
                {
                    max_a = abs_a;
                    max_i = k;
                }
            }

            // failure, matrix is degenerate
            if (max_a < tolerance)
                return false;

            if (max_i != i)
            {
                //pivoting P
                Tools::Swap(p[i], p[max_i]);

                //pivoting rows of A
                Tools::Swap(lu[i], lu[max_i]);

                //counting pivots starting from N (for determinant)
                ++p[N];
            }

            for (SizeT j = i + 1; j < N; j++)
            {
                lu[j][i] /= lu[i][i];

                for (SizeT k = i + 1; k < N; k++)
                    lu[j][k] -= lu[j][i] * lu[i][k];
            }
        }

        return true;
    }

    template <SizeT M, SizeT N>
    VectorN<N> MatrixMxN<M, N>::LUPSolve(const MatrixMxN& lu, const std::vector<SizeT>& p, const VectorN<N>& b)
    {
        VectorN<N> x;
        for (SizeT i = 0; i < N; ++i)
        {
            x[i] = b[p[i]];

            for (SizeT k = 0; k < i; ++k)
            {
                x[i] -= lu[i][k] * x[k];
            }
        }

        for (int ri = N - 1; ri >= 0; --ri)
        {
            SizeT i = static_cast<SizeT>(ri);
            for (SizeT k = i + 1; k < N; ++k)
            {
                x[i] -= lu[i][k] * x[k];
            }

            x[i] /= lu[i][i];
        }

        return x;
    }

    template <SizeT M, SizeT N>
    MatrixMxN<M, N> MatrixMxN<M, N>::LUPInverse(const MatrixMxN& lu, const std::vector<SizeT>& p)
    {
        MatrixMxN invert;
        for (SizeT j = 0; j < N; ++j)
        {
            for (SizeT i = 0; i < N; ++i)
            {
                invert[i][j] = p[i] == j ? 1.0 : 0.0;

                for (SizeT k = 0; k < i; ++k)
                {
                    invert[i][j] -= lu[i][k] * invert[k][j];
                }
            }

            for (int ri = N - 1; ri >= 0; --ri)
            {
                SizeT i = static_cast<SizeT>(ri);
                for (SizeT k = i + 1; k < N; ++k)
                {
                    invert[i][j] -= lu[i][k] * invert[k][j];
                }

                invert[i][j] /= lu[i][i];
            }
        }

        return invert;
    }

    template <SizeT M, SizeT N>
    Real MatrixMxN<M, N>::LUPDeterminant(const MatrixMxN& lu, const std::vector<SizeT>& p)
    {
        double determinant = lu[0][0];
        for (SizeT i = 1; i < N; ++i)
        {
            determinant *= lu[i][i];
        }

        return (p[N] - N) % 2 == 0 ? determinant : -determinant;
    }

    template <SizeT M, SizeT N>
    auto MatrixMxN<M, N>::begin()
    {
        return m_elements.begin();
    }

    template <SizeT M, SizeT N>
    auto MatrixMxN<M, N>::end()
    {
        return m_elements.end();
    }

    template <SizeT M, SizeT N>
    auto MatrixMxN<M, N>::begin() const
    {
        return m_elements.begin();
    }

    template <SizeT M, SizeT N>
    auto MatrixMxN<M, N>::end() const
    {
        return m_elements.end();
    }

    // IO operator
    template <SizeT M, SizeT N>
    std::ostream& operator<<(std::ostream& os, const MatrixMxN<M, N>& rhs)
    {
        os << "{\n";
        for (auto& row_vec : rhs)
        {
            os << "  " << row_vec << ",\n";
        }
        os << "}";
        return os;
    }
}
