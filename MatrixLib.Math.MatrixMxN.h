#pragma once
#include <array>
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
        Real        operator()(SizeT i) const;
        Real&       operator()(SizeT i);
        Real        operator()(SizeT row, SizeT col) const;
        Real&       operator()(SizeT row, SizeT col);
        VectorN<N>  operator[](SizeT row) const;
        VectorN<N>& operator[](SizeT row);

        // assignment operator
        MatrixMxN& operator =(const MatrixMxN& rhs);
        MatrixMxN& operator =(Real rhs);

        // compare operator
        bool operator ==(const MatrixMxN& rhs) const;
        bool operator !=(const MatrixMxN& rhs) const;

        // arithmetic operator
        MatrixMxN  operator -() const;
        MatrixMxN& operator +=(const MatrixMxN& rhs);
        MatrixMxN& operator -=(const MatrixMxN& rhs);
        MatrixMxN& operator *=(Real r);
        MatrixMxN& operator /=(Real r);

    public: // modify methods
        void Set(const MatrixMxN& mat);
        void SetRow(SizeT row, VectorN<N> row_vector);
        void SetColumn(SizeT col, VectorN<M> column_vector);
        void Negate();
        void ClearDigit(SizeT digit);
        void ClearError(Real epsilon);

    public:
        bool            IsEqual(const MatrixMxN& v) const;
        bool            IsNotEqual(const MatrixMxN& v) const;
        bool            IsZero() const;
        bool            IsIdentity() const;
        bool            IsInvertible() const;
        SizeT           Rank() const;
        Real            Trace() const;
        Real            Determinant() const;
        VectorN<N>      Row(SizeT row) const;
        VectorN<M>      Column(SizeT col) const;
        MatrixMxN<N, M> Invert(bool use_permute = true, Real tolerance = Constant::EPSILON) const;
        MatrixMxN<N, M> Transposed() const;
        MatrixMxN       Negated() const;

    public: // range based for loop related methods
        auto begin();
        auto end();
        auto begin() const;
        auto end() const;

    private: // private implementations 
        static MatrixMxN<N, M> TransposeImpl(const MatrixMxN& input_mat_g);
        static SizeT           RankImpl(const MatrixMxN& input_mat_g, Real tolerance);
        static Real            DeterminantImpl(const MatrixMxN& input_mat_g);
        static MatrixMxN       InverseImpl(const MatrixMxN& input_mat_g, bool use_permute);
        static MatrixMxN<N, M> PseudoInverseImpl(const MatrixMxN& input_mat_g, Real tolerance);

        template <SizeT A, SizeT B, SizeT C>
        static MatrixMxN<B, C> SolveImpl(const MatrixMxN<A, B>& mat_a, const MatrixMxN<A, C>& mat_b);

        template <SizeT A, SizeT B, SizeT C>
        static MatrixMxN<A, C>         MultiplyImpl(MatrixMxN<A, B> mat_a, MatrixMxN<B, C> mat_b);
        static MatrixMxN<M - 1, N - 1> MinorMatrixImpl(const MatrixMxN& input_mat_g, SizeT row, SizeT col);

        static MatrixMxN HadamardProductImpl(const MatrixMxN& mat_a, const MatrixMxN& mat_b);
        static MatrixMxN OuterProductImpl(const VectorN<M>& vec_a, const VectorN<N>& vec_b);

    private: // Free-Size Matrix for some calculations;
        class MatrixFxF
        {
        public:
            MatrixFxF()  = default;
            ~MatrixFxF() = default;

            MatrixFxF(SizeT row, SizeT col)
            {
                data = std::vector(row, std::vector(col, Constant::ZERO));
            }

            template <SizeT Row, SizeT Col>
            explicit MatrixFxF(const MatrixMxN<Row, Col>& mat_mxn)
            {
                data = std::vector(Row, std::vector(Col, Constant::ZERO));
                for (SizeT row = 0; row < Row; ++row)
                {
                    for (SizeT col = 0; col < Col; ++col)
                    {
                        data[row][col] = mat_mxn[row][col];
                    }
                }
            }

        public:
            static MatrixFxF Multiply(const MatrixFxF& a, const MatrixFxF& b)
            {
                SizeT a_row_count = a.data.size();
                SizeT a_col_count = a.data[0].size();
                SizeT b_row_count = b.data.size();
                SizeT b_col_count = b.data[0].size();

                if (a_col_count != b_row_count)
                {
                    return MatrixFxF();
                }

                MatrixFxF multiplied(a_row_count, b_col_count);
                for (SizeT i = 0; i < a_row_count; ++i)
                {
                    for (SizeT j = 0; j < b_col_count; ++j)
                    {
                        for (SizeT k = 0; k < a_col_count; ++k)
                        {
                            multiplied.data[i][j] += a.data[i][k] * b.data[k][j];
                        }
                    }
                }

                return multiplied;
            }

            static MatrixFxF Transpose(const MatrixFxF& a)
            {
                SizeT row_count = a.data.size();
                SizeT col_count = a.data[0].size();

                MatrixFxF transposed(col_count, row_count);

                for (SizeT i = 0; i < row_count; ++i)
                {
                    for (SizeT j = 0; j < col_count; ++j)
                    {
                        transposed.data[j][i] = a.data[i][j];
                    }
                }
                return transposed;
            }

            template <SizeT Row, SizeT Col>
            MatrixMxN<Row, Col> ToMatrixMxN()
            {
                MatrixMxN<Row, Col> result;
                for (SizeT row = 0; row < Row; ++row)
                {
                    for (SizeT col = 0; col < Col; ++col)
                    {
                        result[row][col] = data[row][col];
                    }
                }

                return result;
            }

        public:
            std::vector<std::vector<Real>> data;
        };

    public: // friend
        friend class Matrix;
        friend class Vector;
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
    Real MatrixMxN<M, N>::operator()(SizeT i) const
    {
        SizeT row    = i / N;
        SizeT column = i % N;
        return m_elements[row][column];
    }

    template <SizeT M, SizeT N>
    Real& MatrixMxN<M, N>::operator()(SizeT i)
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
    VectorN<N> MatrixMxN<M, N>::operator[](SizeT row) const
    {
        return m_elements[row % M];
    }

    template <SizeT M, SizeT N>
    VectorN<N>& MatrixMxN<M, N>::operator[](SizeT row)
    {
        return m_elements[row % M];
    }

    template <SizeT M, SizeT N>
    MatrixMxN<M, N>& MatrixMxN<M, N>::operator=(const MatrixMxN& rhs)
    {
        if (this != &rhs)
            Set(rhs);

        return *this;
    }

    template <SizeT M, SizeT N>
    MatrixMxN<M, N>& MatrixMxN<M, N>::operator=(Real rhs)
    {
        m_elements.fill(VectorN<N>(rhs));

        return *this;
    }

    template <SizeT M, SizeT N>
    bool MatrixMxN<M, N>::operator==(const MatrixMxN& rhs) const
    {
        return IsEqual(rhs);
    }

    template <SizeT M, SizeT N>
    bool MatrixMxN<M, N>::operator!=(const MatrixMxN& rhs) const
    {
        return IsNotEqual(rhs);
    }

    template <SizeT M, SizeT N>
    MatrixMxN<M, N> MatrixMxN<M, N>::operator-() const
    {
        return Negated();
    }

    template <SizeT M, SizeT N>
    MatrixMxN<M, N>& MatrixMxN<M, N>::operator+=(const MatrixMxN& rhs)
    {
        for (SizeT row = 0; row < M; ++row)
        {
            m_elements[row] += rhs[row];
        }

        return *this;
    }

    template <SizeT M, SizeT N>
    MatrixMxN<M, N>& MatrixMxN<M, N>::operator-=(const MatrixMxN& rhs)
    {
        for (SizeT row = 0; row < M; ++row)
        {
            m_elements[row] -= rhs[row];
        }

        return *this;
    }

    template <SizeT M, SizeT N>
    MatrixMxN<M, N>& MatrixMxN<M, N>::operator*=(Real r)
    {
        for (SizeT row = 0; row < M; ++row)
        {
            m_elements[row] *= r;
        }

        return *this;
    }

    template <SizeT M, SizeT N>
    MatrixMxN<M, N>& MatrixMxN<M, N>::operator/=(Real r)
    {
        for (SizeT row = 0; row < M; ++row)
        {
            m_elements[row] /= r;
        }

        return *this;
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
            m_elements[row][col] = column_vector[row];
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
    bool MatrixMxN<M, N>::IsEqual(const MatrixMxN& v) const
    {
        for (SizeT i = 0; i < M; ++i)
        {
            if (!Tools::IsEqual(m_elements[i], v.m_elements[i]))
                return false;
        }

        return true;
    }

    template <SizeT M, SizeT N>
    bool MatrixMxN<M, N>::IsNotEqual(const MatrixMxN& v) const
    {
        for (SizeT i = 0; i < M; ++i)
        {
            if (!Tools::IsEqual(m_elements[i], v.m_elements[i]))
                return true;
        }

        return false;
    }

    template <SizeT M, SizeT N>
    bool MatrixMxN<M, N>::IsZero() const
    {
        for (SizeT i = 0; i < M; ++i)
        {
            if (!Tools::IsZero(m_elements[i]))
                return false;
        }

        return true;
    }

    template <SizeT M, SizeT N>
    bool MatrixMxN<M, N>::IsIdentity() const
    {
        for (SizeT row = 0; row < M; ++row)
        {
            for (SizeT col = 0; col < N; ++col)
            {
                if (row == col && !Tools::IsEqual(m_elements[row][col], Constant::ONE))
                {
                    return false;
                }

                if (!Tools::IsZero(m_elements[row][col]))
                    return false;
            }
        }

        return true;
    }

    template <SizeT M, SizeT N>
    bool MatrixMxN<M, N>::IsInvertible() const
    {
        return !Tools::IsZero(Determinant());
    }

    template <SizeT M, SizeT N>
    SizeT MatrixMxN<M, N>::Rank() const
    {
        return RankImpl(*this);
    }

    template <SizeT M, SizeT N>
    Real MatrixMxN<M, N>::Trace() const
    {
        // https://en.wikipedia.org/wiki/Trace_(linear_algebra)

        Real trace = Constant::ZERO;
        if constexpr (M == N)
        {
            for (SizeT i = 0; i < M; ++i)
            {
                trace += m_elements[i][i];
            }
        }

        return trace;
    }

    template <SizeT M, SizeT N>
    Real MatrixMxN<M, N>::Determinant() const
    {
        return DeterminantImpl(*this);
    }

    template <SizeT M, SizeT N>
    VectorN<N> MatrixMxN<M, N>::Row(SizeT row) const
    {
        return m_elements[row];
    }

    template <SizeT M, SizeT N>
    VectorN<M> MatrixMxN<M, N>::Column(SizeT col) const
    {
        VectorN<M> column;
        for (SizeT row = 0; row < M; ++row)
        {
            column[row] = m_elements[row][col];
        }

        return column;
    }

    template <SizeT M, SizeT N>
    MatrixMxN<N, M> MatrixMxN<M, N>::Invert(bool use_permute, Real tolerance) const
    {
        if constexpr (M == N)
            return InverseImpl(*this, use_permute);
        else
            return PseudoInverseImpl(*this, tolerance);
    }

    template <SizeT M, SizeT N>
    MatrixMxN<N, M> MatrixMxN<M, N>::Transposed() const
    {
        return TransposeImpl(*this);
    }

    template <SizeT M, SizeT N>
    MatrixMxN<M, N> MatrixMxN<M, N>::Negated() const
    {
        MatrixMxN negated;
        for (SizeT i = 0; i < N; ++i)
        {
            negated.m_elements[i] = -m_elements[i];
        }

        return negated;
    }

    template <SizeT M, SizeT N>
    MatrixMxN<N, M> MatrixMxN<M, N>::TransposeImpl(const MatrixMxN& input_mat_g)
    {
        MatrixMxN<N, M> transposed;
        for (SizeT i = 0; i < M; ++i)
        {
            transposed.SetColumn(i, input_mat_g.Row(i));
        }

        return transposed;
    }

    template <SizeT M, SizeT N>
    SizeT MatrixMxN<M, N>::RankImpl(const MatrixMxN& input_mat_g, Real tolerance)
    {
        // https://en.wikipedia.org/wiki/Cholesky_decomposition
        // Cholesky decomposition	
        constexpr SizeT       size = M < N ? M : N;
        MatrixMxN<size, size> mat_a;

        if (M < N)
        {
            // A = G * G'
            for (SizeT i = 0; i < size; ++i)
            {
                for (SizeT j = 0; j < size; ++j)
                {
                    for (SizeT k = 0; k < N; ++k)
                    {
                        mat_a[i][j] += input_mat_g[i][k] * input_mat_g[j][k];
                    }
                }
            }
        }
        else
        {
            // A = G' * G
            for (SizeT i = 0; i < size; ++i)
            {
                for (SizeT j = 0; j < size; ++j)
                {
                    for (SizeT k = 0; k < N; ++k)
                    {
                        mat_a[i][j] += input_mat_g[k][i] * input_mat_g[k][j];
                    }
                }
            }
        }

        // Full Rank Cholesky decomposition of A
        Real tol = Tools::Abs(mat_a[0][0]);

        for (SizeT i = 0; i < size; ++i)
        {
            if (mat_a[i][i] > 0)
            {
                tol = Tools::Min(Tools::Abs(mat_a[i][i]), tol);
            }
        }
        tol *= tolerance;

        MatrixMxN<size, size> mat_l;

        SizeT rank_a = 0;
        for (SizeT k = 0; k < size; ++k)
        {
            for (SizeT i = k; i < size; ++i)
            {
                mat_l[i][rank_a] = mat_a[i][k];
                for (SizeT j = 0; j < rank_a; ++j)
                {
                    mat_l[i][rank_a] -= mat_l[i][j] * mat_l[k][j];
                }
            }

            if (mat_l[k][rank_a] > tol)
            {
                mat_l[k][rank_a] = std::sqrt(mat_l[k][rank_a]);
                if (k < size)
                {
                    for (SizeT j = k + 1; j < size; ++j)
                    {
                        mat_l[j][rank_a] /= mat_l[k][rank_a];
                    }
                }

                ++rank_a;
            }
        }

        return rank_a;
    }

    template <SizeT M, SizeT N>
    Real MatrixMxN<M, N>::DeterminantImpl(const MatrixMxN& input_mat_g)
    {
        // https://en.wikipedia.org/wiki/LU_decomposition

        // 0. Check Square Matrix
        if constexpr (M != N)
            return Constant::ZERO;

        // 1. Row Permutation
        MatrixMxN mat_lu;
        bool      change_sign = false;
        {
            std::vector<SizeT> permute_lu;
            for (SizeT i = 0; i < N; ++i)
            {
                permute_lu.push_back(i);
            }

            for (SizeT j = 0; j < N; ++j)
            {
                Real max_value = Constant::ZERO;
                for (SizeT i = j; i < N; ++i)
                {
                    Real current_v = Tools::Abs(input_mat_g[permute_lu[i]][j]);
                    if (current_v > max_value)
                    {
                        max_value = current_v;
                        if (permute_lu[i] != permute_lu[j])
                        {
                            change_sign = !change_sign;

                            // swap rows
                            SizeT temp    = permute_lu[j];
                            permute_lu[j] = permute_lu[i];
                            permute_lu[i] = temp;
                        }
                    }
                }
            }

            for (SizeT i = 0; i < N; ++i)
            {
                mat_lu.SetRow(i, input_mat_g[permute_lu[i]]);
            }
        }

        // 2. LU decomposition
        {
            if (Tools::IsZero(mat_lu[0][0]))
            {
                // Singular matrix
                return Constant::ZERO;
            }

            for (SizeT i = 1; i < N; ++i)
            {
                mat_lu[i][0] /= mat_lu[0][0];
            }

            for (SizeT i = 1; i < N; ++i)
            {
                for (SizeT j = i; j < N; ++j)
                {
                    for (SizeT k = 0; k < i; ++k)
                    {
                        // Calculate U matrix
                        mat_lu[i][j] -= mat_lu[i][k] * mat_lu[k][j];
                    }
                }

                if (Tools::IsZero(mat_lu[i][i]))
                {
                    // Singular matrix
                    return Constant::ZERO;
                }

                for (SizeT k = i + 1; k < N; ++k)
                {
                    for (SizeT j = 0; j < i; ++j)
                    {
                        // Calculate L matrix
                        mat_lu[k][i] -= mat_lu[k][j] * mat_lu[j][i];
                    }
                    mat_lu[k][i] /= mat_lu[i][i];
                }
            }
        }

        Real determinant = change_sign ? -Constant::ONE : Constant::ONE;
        for (SizeT i = 0; i < N; ++i)
        {
            determinant *= mat_lu[i][i];
        }

        return determinant;
    }

    template <SizeT M, SizeT N>
    MatrixMxN<M, N> MatrixMxN<M, N>::InverseImpl(const MatrixMxN& input_mat_g, bool use_permute)
    {
        // https://en.wikipedia.org/wiki/LU_decomposition
        // 0. Check Square Matrix
        if constexpr (M != N)
            return MatrixMxN();

        // 1. Row Permutation   
        MatrixMxN          mat_lu;
        std::vector<SizeT> permute_lu;
        {
            for (SizeT i = 0; i < N; ++i)
            {
                // Set up row index
                permute_lu.push_back(i);
            }

            if (use_permute)
            {
                // Sort rows by pivot element
                for (SizeT j = 0; j < N; ++j)
                {
                    Real max_v = Constant::ZERO;
                    for (SizeT i = j; i < N; ++i)
                    {
                        Real current_v = Tools::Abs(input_mat_g[permute_lu[i]][j]);
                        if (current_v > max_v)
                        {
                            max_v = current_v;

                            // Swap rows
                            SizeT temp    = permute_lu[j];
                            permute_lu[j] = permute_lu[i];
                            permute_lu[i] = temp;
                        }
                    }
                }

                for (SizeT i = 0; i < N; ++i)
                {
                    // Make a permuted matrix with new row order
                    mat_lu.SetRow(i, input_mat_g[permute_lu[i]]);
                }
            }
            else
            {
                // duplicate matrix
                mat_lu = input_mat_g;
            }
        }

        // 2. LU decomposition
        {
            if (Tools::IsZero(mat_lu[0][0]))
            {
                // Singular matrix
                return MatrixMxN();
            }

            // Initialize first column of L matrix
            for (SizeT i = 1; i < N; ++i)
            {
                mat_lu[i][0] /= mat_lu[0][0];
            }

            for (SizeT i = 1; i < N; ++i)
            {
                for (SizeT j = i; j < N; ++j)
                {
                    for (SizeT k = 0; k < i; ++k)
                    {
                        // Calculate U matrix
                        mat_lu[i][j] -= mat_lu[i][k] * mat_lu[k][j];
                    }
                }

                if (Tools::IsZero(mat_lu[i][i]))
                {
                    // Singular matrix
                    return MatrixMxN();
                }

                for (SizeT k = i + 1; k < N; ++k)
                {
                    for (SizeT j = 0; j < i; ++j)
                    {
                        // Calculate L matrix
                        mat_lu[k][i] -= mat_lu[k][j] * mat_lu[j][i];
                    }
                    mat_lu[k][i] /= mat_lu[i][i];
                }
            }
        }

        // 3. L & U inversion
        MatrixMxN mat_inv_lu;
        {
            // mat L inverse & mat U inverse
            for (SizeT i = 0; i < N; ++i)
            {
                // L matrix inverse, omit diagonal ones
                mat_inv_lu[i][i] = Constant::ONE;
                for (SizeT k = i + 1; k < N; ++k)
                {
                    for (SizeT j = i; j <= k - 1; ++j)
                    {
                        mat_inv_lu[k][i] -= mat_lu[k][j] * mat_inv_lu[j][i];
                    }
                }
                // U matrix inverse
                mat_inv_lu[i][i] = Tools::Inverse(mat_lu[i][i]);
                for (SizeT k = i; k > 0; --k)
                {
                    for (SizeT j = k; j <= i; ++j)
                    {
                        mat_inv_lu[k - 1][i] -= mat_lu[k - 1][j] * mat_inv_lu[j][i];
                    }
                    mat_inv_lu[k - 1][i] /= mat_lu[k - 1][k - 1];
                }
            }
        }

        // 4. Calculate G^-1 = U^-1 * L^-1
        {
            // 4.1. Lower part product
            for (SizeT i = 1; i < N; ++i)
            {
                for (SizeT j = 0; j < i; ++j)
                {
                    // Permute column back
                    SizeT j_p      = permute_lu[j];
                    mat_lu[i][j_p] = Constant::ZERO;
                    for (SizeT k = i; k < N; ++k)
                    {
                        mat_lu[i][j_p] += mat_inv_lu[i][k] * mat_inv_lu[k][j];
                    }
                }
            }

            // 4.2. Upper part product
            for (SizeT i = 0; i < N; ++i)
            {
                for (SizeT j = i; j < N; ++j)
                {
                    // Permute column back
                    SizeT j_p      = permute_lu[j];
                    mat_lu[i][j_p] = mat_inv_lu[i][j];
                    for (SizeT k = j + 1; k < N; ++k)
                    {
                        mat_lu[i][j_p] += mat_inv_lu[i][k] * mat_inv_lu[k][j];
                    }
                }
            }
        }

        // inverse of input matrix
        return mat_lu;
    }

    template <SizeT M, SizeT N>
    MatrixMxN<N, M> MatrixMxN<M, N>::PseudoInverseImpl(const MatrixMxN& input_mat_g, Real tolerance)
    {
        // https://en.wikipedia.org/wiki/Moore-Penrose_inverse
        constexpr SizeT size = M < N ? M : N;

        MatrixFxF mat_a;
        MatrixFxF mat_g(input_mat_g);
        MatrixFxF mat_g_transposed = MatrixFxF::Transpose(mat_g);

        if constexpr (M < N)
        {
            // A = G * Gt
            mat_a = MatrixFxF::Multiply(mat_g, mat_g_transposed);
        }
        else
        {
            // A = Gt * G
            mat_a = MatrixFxF::Multiply(mat_g_transposed, mat_g);
        }

        // Full rank Cholesky decomposition of A
        MatrixFxF mat_l(size, size);
        SizeT     rank_a = 0;
        {
            Real tol = Tools::Abs(mat_a.data[0][0]);
            for (SizeT i = 0; i < size; ++i)
            {
                if (mat_a.data[i][i] > Constant::ZERO)
                {
                    tol = Tools::Min(mat_a.data[i][i], tol);
                }
            }

            tol *= tolerance;

            for (SizeT k = 0; k < size; ++k)
            {
                for (SizeT i = k; i < size; ++i)
                {
                    mat_l.data[i][rank_a] = mat_a.data[i][k];
                    for (SizeT j = 0; j < rank_a; ++j)
                    {
                        mat_l.data[i][rank_a] -= mat_l.data[i][j] * mat_l.data[k][j];
                    }
                }

                if (mat_l.data[k][rank_a] > tol)
                {
                    mat_l.data[k][rank_a] = Tools::Sqrt(mat_l.data[k][rank_a]);

                    if (k < size)
                    {
                        for (SizeT j = k + 1; j < size; ++j)
                        {
                            mat_l.data[j][rank_a] /= mat_l.data[k][rank_a];
                        }
                    }

                    ++rank_a;
                }
            }
        }

        // Zero Rank No Inverse
        if (rank_a == 0)
        {
            // All-Zero transposed
            return MatrixMxN<N, M>();
        }

        // Slice L = L(:, 0:r);
        for (SizeT i = 0; i < size; ++i)
        {
            for (SizeT k = 0; k < size - rank_a; ++k)
            {
                mat_l.data[i].pop_back();
            }
        }

        // Generalized inverse
        MatrixFxF             mat_l_transposed = MatrixFxF::Transpose(mat_l);
        MatrixFxF             mat_lt_l         = MatrixFxF::Multiply(mat_l_transposed, mat_l);
        MatrixMxN<size, size> mat_lt_l_fixed   = mat_lt_l.template ToMatrixMxN<size, size>();
        MatrixMxN<size, size> mat_inv_lt_l     = mat_lt_l_fixed.Invert(false);

        // M = inv(Lt * L)
        MatrixFxF mat_m(mat_inv_lt_l);

        // inv = L * M * M * Lt
        MatrixFxF pseudo_invert = MatrixFxF::Multiply(MatrixFxF::Multiply(MatrixFxF::Multiply(mat_l, mat_m), mat_m), mat_l_transposed);

        if constexpr (M < N)
        {
            // pseudo_inverse(G) = Gt * (L * M * M * Lt)
            return MatrixFxF::Multiply(mat_g_transposed, pseudo_invert).template ToMatrixMxN<N, M>();
        }
        else
        {
            // pseudo_inverse(G) = (L * M * M * Lt) * Gt
            return MatrixFxF::Multiply(pseudo_invert, mat_g_transposed).template ToMatrixMxN<N, M>();
        }
    }

    template <SizeT M, SizeT N>
    template <SizeT A, SizeT B, SizeT C>
    MatrixMxN<B, C> MatrixMxN<M, N>::SolveImpl(const MatrixMxN<A, B>& mat_a, const MatrixMxN<A, C>& mat_b)
    {
        // Solve A * x = b;
        // x = A^-1 * b;
        MatrixMxN<B, A> inv_a = MatrixMxN::PseudoInverseImpl(mat_a, Constant::EPSILON);
        MatrixMxN<B, C> x     = MatrixMxN::MultiplyImpl(inv_a, mat_b);
        return x;
    }

    template <SizeT M, SizeT N>
    template <SizeT A, SizeT B, SizeT C>
    MatrixMxN<A, C> MatrixMxN<M, N>::MultiplyImpl(MatrixMxN<A, B> mat_a, MatrixMxN<B, C> mat_b)
    {
        MatrixMxN<A, C> multiplied;
        for (SizeT row = 0; row < A; ++row)
        {
            for (SizeT col = 0; col < C; ++col)
            {
                multiplied(row, col) = Constant::ZERO;
                VectorN<B> row_vec   = mat_a.Row(row);
                VectorN<B> col_vec   = mat_b.Column(col);

                for (SizeT i = 0; i < B; ++i)
                {
                    multiplied(row, col) += row_vec[i] * col_vec[i];
                }
            }
        }

        return multiplied;
    }

    template <SizeT M, SizeT N>
    MatrixMxN<M - 1, N - 1> MatrixMxN<M, N>::MinorMatrixImpl(const MatrixMxN& input_mat_g, SizeT row, SizeT col)
    {
        // https://en.wikipedia.org/wiki/Minor_(linear_algebra)
        MatrixMxN<M - 1, N - 1> minor;
        for (SizeT i = 0; i < M; ++i)
        {
            SizeT mi = i - row;
            if (mi == 0)
            {
                // skip this row
                continue;
            }

            SizeT minor_i = 0 < mi && mi < M ? i - 1 : i;

            for (SizeT j = 0; j < N; ++j)
            {
                SizeT mj = j - col;
                if (mj == 0)
                {
                    // skip this col
                    continue;
                }

                SizeT minor_j = 0 < mj && mj < N ? j - 1 : j;

                minor[minor_i][minor_j] = input_mat_g[i][j];
            }
        }

        return minor;
    }

    template <SizeT M, SizeT N>
    MatrixMxN<M, N> MatrixMxN<M, N>::HadamardProductImpl(const MatrixMxN& mat_a, const MatrixMxN& mat_b)
    {
        MatrixMxN hadamard;
        for (SizeT row = 0; row < M; ++row)
        {
            for (SizeT col = 0; col < N; ++col)
            {
                hadamard[row][col] = mat_a[row][col] * mat_b[row][col];
            }
        }

        return hadamard;
    }

    template <SizeT M, SizeT N>
    MatrixMxN<M, N> MatrixMxN<M, N>::OuterProductImpl(const VectorN<M>& vec_a, const VectorN<N>& vec_b)
    {
        MatrixMxN outer;
        for (SizeT row = 0; row < M; ++row)
        {
            for (SizeT col = 0; col < N; ++col)
            {
                outer[row][col] = vec_a[row] * vec_b[col];
            }
        }

        return outer;
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
}
