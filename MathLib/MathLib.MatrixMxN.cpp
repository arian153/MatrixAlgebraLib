#include "MathLib.MatrixMxN.h"

namespace MathLib
{
    MatrixMxN::MatrixMxN()
    {
    }

    MatrixMxN::MatrixMxN(std::initializer_list<Real> values)
    {
    }

    MatrixMxN::MatrixMxN(std::initializer_list<std::initializer_list<Real>> values)
    {
    }

    MatrixMxN::MatrixMxN(std::initializer_list<VectorN> values)
    {
    }

    MatrixMxN::MatrixMxN(const MatrixMxN& rhs)
    {
    }

    MatrixMxN::MatrixMxN(Real scalar)
    {
    }

    MatrixMxN::MatrixMxN(const std::vector<Real>& vec)
    {
    }

    MatrixMxN::MatrixMxN(const std::vector<VectorN>& vec)
    {
    }

    MatrixMxN::operator std::vector<float>() const
    {
    }

    Real MatrixMxN::operator()(SizeT i) const
    {
    }

    Real& MatrixMxN::operator()(SizeT i)
    {
    }

    Real MatrixMxN::operator()(SizeT row, SizeT col) const
    {
    }

    Real& MatrixMxN::operator()(SizeT row, SizeT col)
    {
    }

    VectorN MatrixMxN::operator[](SizeT row) const
    {
    }

    VectorN& MatrixMxN::operator[](SizeT row)
    {
    }

    MatrixMxN& MatrixMxN::operator=(const MatrixMxN& rhs)
    {
    }

    MatrixMxN& MatrixMxN::operator=(Real rhs)
    {
    }

    bool MatrixMxN::operator==(const MatrixMxN& rhs) const
    {
    }

    bool MatrixMxN::operator!=(const MatrixMxN& rhs) const
    {
    }

    MatrixMxN MatrixMxN::operator-() const
    {
    }

    MatrixMxN& MatrixMxN::operator+=(const MatrixMxN& rhs)
    {
    }

    MatrixMxN& MatrixMxN::operator-=(const MatrixMxN& rhs)
    {
    }

    MatrixMxN& MatrixMxN::operator*=(Real r)
    {
    }

    MatrixMxN& MatrixMxN::operator/=(Real r)
    {
    }

    void MatrixMxN::Set(const MatrixMxN& mat)
    {
    }

    void MatrixMxN::SetRow(SizeT row, const VectorN& row_vector)
    {
    }

    void MatrixMxN::SetColumn(SizeT col, const VectorN& column_vector)
    {
    }

    void MatrixMxN::Negate()
    {
    }

    void MatrixMxN::ClearDigit(SizeT digit)
    {
    }

    void MatrixMxN::ClearError(Real epsilon)
    {
    }

    bool MatrixMxN::IsEqual(const MatrixMxN& v) const
    {
    }

    bool MatrixMxN::IsNotEqual(const MatrixMxN& v) const
    {
    }

    bool MatrixMxN::IsZero() const
    {
    }

    bool MatrixMxN::IsIdentity() const
    {
    }

    bool MatrixMxN::IsInvertible() const
    {
    }

    SizeT MatrixMxN::Rank() const
    {
    }

    SizeT MatrixMxN::SmallestIdx() const
    {
    }

    SizeT MatrixMxN::LargestIdx() const
    {
    }

    Real MatrixMxN::Trace() const
    {
    }

    Real MatrixMxN::Determinant() const
    {
    }

    VectorN MatrixMxN::Row(SizeT row) const
    {
    }

    VectorN MatrixMxN::Column(SizeT col) const
    {
    }

    MatrixMxN MatrixMxN::Invert(bool use_permute, Real tolerance) const
    {
    }

    MatrixMxN MatrixMxN::Transposed() const
    {
    }

    MatrixMxN MatrixMxN::Negated() const
    {
    }

    constexpr MatrixMxN MatrixMxN::Identity()
    {
    }

    auto MatrixMxN::begin()
    {
    }

    auto MatrixMxN::end()
    {
    }

    auto MatrixMxN::begin() const
    {
    }

    auto MatrixMxN::end() const
    {
    }

    MatrixMxN MatrixMxN::TransposeImpl(const MatrixMxN& input_mat_g)
    {
    }

    SizeT MatrixMxN::RankImpl(const MatrixMxN& input_mat_g, Real tolerance)
    {
    }

    Real MatrixMxN::DeterminantImpl(const MatrixMxN& input_mat_g)
    {
    }

    MatrixMxN MatrixMxN::InverseImpl(const MatrixMxN& input_mat_g, bool use_permute)
    {
    }

    MatrixMxN MatrixMxN::PseudoInverseImpl(const MatrixMxN& input_mat_g, Real tolerance)
    {
    }

    MatrixMxN MatrixMxN::MinorMatrixImpl(const MatrixMxN& input_mat_g, SizeT row, SizeT col)
    {
    }
}
