#include "MathLib.VectorN.h"

namespace MathLib
{
    VectorN::VectorN()
    {
    }

    VectorN::VectorN(std::initializer_list<Real> values)
    {
    }

    VectorN::VectorN(const VectorN& rhs)
    {
    }

    VectorN::VectorN(Real scalar)
    {
    }

    VectorN::VectorN(const std::vector<Real>& vec)
    {
    }

    VectorN::operator std::vector<float>() const
    {
    }

    Real VectorN::operator[](SizeT i) const
    {
    }

    Real& VectorN::operator[](SizeT i)
    {
    }

    Real VectorN::operator()(SizeT i) const
    {
    }

    Real& VectorN::operator()(SizeT i)
    {
    }

    VectorN& VectorN::operator=(const VectorN& rhs)
    {
    }

    VectorN& VectorN::operator=(Real rhs)
    {
    }

    bool VectorN::operator==(const VectorN& rhs) const
    {
    }

    bool VectorN::operator!=(const VectorN& rhs) const
    {
    }

    VectorN VectorN::operator-() const
    {
    }

    VectorN& VectorN::operator+=(const VectorN& rhs)
    {
    }

    VectorN& VectorN::operator-=(const VectorN& rhs)
    {
    }

    VectorN& VectorN::operator*=(Real r)
    {
    }

    VectorN& VectorN::operator/=(Real r)
    {
    }

    void VectorN::Set(const VectorN& v)
    {
    }

    void VectorN::SetScalar(Real value)
    {
    }

    void VectorN::SetZero()
    {
    }

    void VectorN::Negate()
    {
    }

    void VectorN::Inverse()
    {
    }

    void VectorN::Normalize()
    {
    }

    void VectorN::ClearDigit(SizeT digit)
    {
    }

    void VectorN::ClearError(Real epsilon)
    {
    }

    bool VectorN::IsEqual(const VectorN& v) const
    {
    }

    bool VectorN::IsNotEqual(const VectorN& v) const
    {
    }

    bool VectorN::IsZero() const
    {
    }

    bool VectorN::IsNormalized() const
    {
    }

    VectorN VectorN::Negated() const
    {
    }

    VectorN VectorN::Invert() const
    {
    }

    VectorN VectorN::Normalized() const
    {
    }

    VectorN VectorN::Scaled(Real r) const
    {
    }

    Real VectorN::Length() const
    {
    }

    Real VectorN::LengthSq() const
    {
    }

    SizeT VectorN::SmallestIdx() const
    {
    }

    SizeT VectorN::LargestIdx() const
    {
    }

    VectorN VectorN::Half() const
    {
    }

    VectorN VectorN::Absolute() const
    {
    }

    Real VectorN::Norm(Real p) const
    {
    }

    VectorN VectorN::Swizzle(std::initializer_list<SizeT> indices_list) const
    {
    }

    Real VectorN::DotProduct(const VectorN& a, const VectorN& b)
    {
    }

    VectorN VectorN::HadamardProduct(const VectorN& a, const VectorN& b)
    {
    }

    VectorN VectorN::CrossProduct3D(const VectorN& a, const VectorN& b)
    {
    }

    VectorN VectorN::CrossProduct7D(const VectorN& a, const VectorN& b)
    {
    }

    auto VectorN::begin()
    {
    }

    auto VectorN::end()
    {
    }

    auto VectorN::begin() const
    {
    }

    auto VectorN::end() const
    {
    }
}
