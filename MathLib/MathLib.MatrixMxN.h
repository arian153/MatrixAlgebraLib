#pragma once
#include <ranges>
#include <vector>

#include "MathLib.Constants.h"
#include "MathLib.Tools.h"
#include "MathLib.VectorN.h"

namespace MathLib
{
    class MatrixMxN
    {
    private: // private data
        std::vector<VectorN> m_elements;

    public: // constructors
        MatrixMxN();
        MatrixMxN(std::initializer_list<Real> values);
        MatrixMxN(std::initializer_list<std::initializer_list<Real>> values);
        MatrixMxN(std::initializer_list<VectorN> values);
        MatrixMxN(const MatrixMxN& rhs);
        ~MatrixMxN() = default;
        explicit MatrixMxN(Real scalar);
        explicit MatrixMxN(const std::vector<Real>& vec);
        explicit MatrixMxN(const std::vector<VectorN>& vec);

    public: // operators
        // explicit operator
        explicit operator std::vector<Real>() const;

        // access operator
        Real     operator()(SizeT i) const;
        Real&    operator()(SizeT i);
        Real     operator()(SizeT row, SizeT col) const;
        Real&    operator()(SizeT row, SizeT col);
        VectorN  operator[](SizeT row) const;
        VectorN& operator[](SizeT row);

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
        void SetRow(SizeT row, const VectorN& row_vector);
        void SetColumn(SizeT col, const VectorN& column_vector);
        void Negate();
        void ClearDigit(SizeT digit);
        void ClearError(Real epsilon);

    public:
        bool      IsEqual(const MatrixMxN& v) const;
        bool      IsNotEqual(const MatrixMxN& v) const;
        bool      IsZero() const;
        bool      IsIdentity() const;
        bool      IsInvertible() const;
        SizeT     Rank() const;
        SizeT     SmallestIdx() const;
        SizeT     LargestIdx() const;
        Real      Trace() const;
        Real      Determinant() const;
        VectorN   Row(SizeT row) const;
        VectorN   Column(SizeT col) const;
        MatrixMxN Invert(bool use_permute = true, Real tolerance = Constant::EPSILON) const;
        MatrixMxN Transposed() const;
        MatrixMxN Negated() const;

        // const
        static constexpr MatrixMxN Identity();

    public: // range based for loop related methods
        auto begin();
        auto end();
        auto begin() const;
        auto end() const;

    private: // private implementations 
        static MatrixMxN TransposeImpl(const MatrixMxN& input_mat_g);
        static SizeT     RankImpl(const MatrixMxN& input_mat_g, Real tolerance);
        static Real      DeterminantImpl(const MatrixMxN& input_mat_g);
        static MatrixMxN InverseImpl(const MatrixMxN& input_mat_g, bool use_permute);
        static MatrixMxN PseudoInverseImpl(const MatrixMxN& input_mat_g, Real tolerance);
        static MatrixMxN MinorMatrixImpl(const MatrixMxN& input_mat_g, SizeT row, SizeT col);
    };

   
}
