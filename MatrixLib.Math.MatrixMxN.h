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
        Real        operator()(SizeT i) const;
        Real&       operator()(SizeT i);
        Real        operator()(SizeT row, SizeT col) const;
        Real&       operator()(SizeT row, SizeT col);
        VectorN<N>  operator[](SizeT row) const;
        VectorN<N>& operator[](SizeT row);

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

        //TODO
        // Transpose
        // Trace
        // Determinant
        // Inverse
        // Adjoint
        // Hadamard Product
        // Multiply MxN * NxM

        VectorN<N> Row(SizeT row) const;
        VectorN<M> Column(SizeT col) const;

    public:
        MatrixMxN Invert() const;

    public:
        static MatrixMxN<N, M> Transpose(const MatrixMxN& mat);
        static SizeT           Rank(const MatrixMxN& mat_g, Real tolerance = Constant::EPSILON);
        static Real            Determinant(const MatrixMxN& mat_g);
        static MatrixMxN       Inverse(const MatrixMxN& mat_g, bool use_permute = true);
        static MatrixMxN<N, M> PseudoInverse(const MatrixMxN& input, Real tolerance = Constant::EPSILON);
        static MatrixMxN       SolveAxEqualB(const MatrixMxN& mat_a, const MatrixMxN& mat_b);

        template <SizeT O>
        static MatrixMxN<M, O> Multiply(MatrixMxN a, MatrixMxN<N, O> b);

    public: // range based for loop related methods
        auto begin();
        auto end();
        auto begin() const;
        auto end() const;

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
    MatrixMxN<M, N> MatrixMxN<M, N>::Invert() const
    {
        MatrixMxN invert;

        return invert;
    }

    template <SizeT M, SizeT N>
    MatrixMxN<N, M> MatrixMxN<M, N>::Transpose(const MatrixMxN& mat)
    {
        MatrixMxN<N, M> transposed;
        for (SizeT i = 0; i < M; ++i)
        {
            transposed.SetColumn(i, mat.Row(i));
        }

        return transposed;
    }

    template <SizeT M, SizeT N>
    SizeT MatrixMxN<M, N>::Rank(const MatrixMxN& mat_g, Real tolerance)
    {
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
                        mat_a[i][j] += mat_g[i][k] * mat_g[j][k];
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
                        mat_a[i][j] += mat_g[k][i] * mat_g[k][j];
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
    Real MatrixMxN<M, N>::Determinant(const MatrixMxN& mat_g)
    {
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
                    Real current_v = Tools::Abs(mat_g[permute_lu[i]][j]);
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
                mat_lu.SetRow(i, mat_g[permute_lu[i]]);
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
    MatrixMxN<M, N> MatrixMxN<M, N>::Inverse(const MatrixMxN& mat_g, bool use_permute)
    {
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
                        Real current_v = Tools::Abs(mat_g[permute_lu[i]][j]);
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
                    mat_lu.SetRow(i, mat_g[permute_lu[i]]);
                }
            }
            else
            {
                // duplicate matrix
                mat_lu = mat_g;
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
    MatrixMxN<N, M> MatrixMxN<M, N>::PseudoInverse(const MatrixMxN& input, Real tolerance)
    {
        constexpr SizeT size = M < N ? M : N;

        MatrixFxF mat_a;
        MatrixFxF mat_g(input);
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
        MatrixMxN<size, size> mat_inv_lt_l     = MatrixMxN<size, size>::Inverse(mat_lt_l.template ToMatrixMxN<size, size>(), false);

        // M = inv(L' * L)
        MatrixFxF mat_m(mat_inv_lt_l);

        // A = L * M * M * L'
        MatrixFxF pseudo_invert = MatrixFxF::Multiply(MatrixFxF::Multiply(MatrixFxF::Multiply(mat_l, mat_m), mat_m), mat_l_transposed);

        if constexpr (M < N)
        {
            // pseudo_inverse(G) = G' * (L * M * M * L')
            return MatrixFxF::Multiply(mat_g_transposed, pseudo_invert).template ToMatrixMxN<N, M>();
        }
        else
        {
            // pseudo_inverse(G) = (L * M * M * L') * G'
            return MatrixFxF::Multiply(pseudo_invert, mat_g_transposed).template ToMatrixMxN<N, M>();
        }
    }

    template <SizeT M, SizeT N>
    MatrixMxN<M, N> MatrixMxN<M, N>::SolveAxEqualB(const MatrixMxN& mat_a, const MatrixMxN& mat_b)
    {
        // Solve A * x = b;
        // x = A^-1 * b;
        MatrixMxN x = Multiply(PseudoInverse(mat_a), mat_b);
        return x;
    }

    template <SizeT M, SizeT N>
    template <SizeT O>
    MatrixMxN<M, O> MatrixMxN<M, N>::Multiply(MatrixMxN a, MatrixMxN<N, O> b)
    {
        MatrixMxN<M, O> multiplied;
        for (SizeT row = 0; row < M; ++row)
        {
            for (SizeT col = 0; col < O; ++col)
            {
                multiplied(row, col) = VectorN<N>::DotProduct(a.Row(row), b.Column(col));
            }
        }

        return multiplied;
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

    // ari

    template <SizeT M, SizeT N, SizeT O>
    MatrixMxN<M, O> Multiply(MatrixMxN<M, N> a, MatrixMxN<N, O> b)
    {
        MatrixMxN<M, O> multiplied;
        for (SizeT row = 0; row < M; ++row)
        {
            for (SizeT col = 0; col < O; ++col)
            {
                multiplied(row, col) = VectorN<N>::DotProduct(a.Row(row), b.Column(col));
            }
        }

        return multiplied;
    }

    // IO operator
    template <SizeT M, SizeT N>
    std::ostream& operator<<(std::ostream& os, const MatrixMxN<M, N>& rhs)
    {
        os << "{\n " << rhs[0];
        for (SizeT i = 1; i < M; ++i)
        {
            os << ",\n " << rhs[i];
        }
        os << "\n}";
        return os;
    }
}
