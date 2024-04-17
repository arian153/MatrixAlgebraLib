#include "MatrixLib.Math.VectorN.h"
#include "MatrixLib.Math.Tools.h"

namespace MatrixLib
{
    VectorN::VectorN(const VectorN& rhs)
        : m_data(rhs.m_data), m_fixed_size(rhs.m_fixed_size)
    {
    }

    VectorN::VectorN(const std::vector<Real>& v)
        : m_data(v), m_fixed_size(v.size())
    {
    }

    VectorN::VectorN(SizeT size_of_vector, Real value)
        : m_fixed_size(size_of_vector)
    {
        m_data.resize(size_of_vector);
        SetScalar(value);
    }

    VectorN::operator std::vector<Real>() const
    {
        return { m_data };
    }

    Real VectorN::operator[](size_t i) const
    {
        return m_data[i];
    }

    Real& VectorN::operator[](size_t i)
    {
        return m_data[i];
    }

    Real VectorN::operator()(size_t i) const
    {
        return m_data[i];
    }

    Real& VectorN::operator()(size_t i)
    {
        return m_data[i];
    }

    VectorN& VectorN::operator=(const VectorN& rhs)
    {
        if (this != &rhs)
            Set(rhs);

        return *this;
    }

    VectorN& VectorN::operator=(Real rhs)
    {
        SetScalar(rhs);
        return *this;
    }

    bool VectorN::operator==(const VectorN& rhs) const
    {
        return IsEqual(rhs);
    }

    bool VectorN::operator!=(const VectorN& rhs) const
    {
        return IsNotEqual(rhs);
    }

    VectorN VectorN::operator-() const
    {
        return Negated();
    }

    VectorN& VectorN::operator+=(const VectorN& rhs)
    {
        Add(rhs);
        return *this;
    }

    VectorN& VectorN::operator-=(const VectorN& rhs)
    {
        Subtract(rhs);
        return *this;
    }

    VectorN& VectorN::operator+=(Real r)
    {
        Add(r);
        return *this;
    }

    VectorN& VectorN::operator-=(Real r)
    {
        Subtract(r);
        return *this;
    }

    VectorN& VectorN::operator*=(Real r)
    {
        Multiply(r);
        return *this;
    }

    VectorN& VectorN::operator/=(Real r)
    {
        Divide(r);
        return *this;
    }

    void VectorN::Set(const VectorN& v)
    {
        m_data       = v.m_data;
        m_fixed_size = v.m_fixed_size;
    }

    void VectorN::SetScalar(Real value)
    {
        for (SizeT i = 0; i < m_fixed_size; ++i)
        {
            m_data[i] = value;
        }
    }

    void VectorN::SetZero()
    {
        for (SizeT i = 0; i < m_fixed_size; ++i)
        {
            m_data[i] = Constant::ZERO;
        }
    }

    void VectorN::Negate()
    {
        for (SizeT i = 0; i < m_fixed_size; ++i)
        {
            m_data[i] = -m_data[i];
        }
    }

    void VectorN::Inverse()
    {
        for (SizeT i = 0; i < m_fixed_size; ++i)
        {
            m_data[i] = Tools::Inverse(m_data[i]);
        }
    }

    void VectorN::Add(const VectorN& v)
    {
        if (m_fixed_size == v.m_fixed_size)
        {
            for (SizeT i = 0; i < m_fixed_size; ++i)
            {
                m_data[i] += v.m_data[i];
            }
        }
    }

    void VectorN::Subtract(const VectorN& v)
    {
        if (m_fixed_size == v.m_fixed_size)
        {
            for (SizeT i = 0; i < m_fixed_size; ++i)
            {
                m_data[i] -= v.m_data[i];
            }
        }
    }

    void VectorN::Add(Real r)
    {
        for (SizeT i = 0; i < m_fixed_size; ++i)
        {
            m_data[i] += r;
        }
    }

    void VectorN::Subtract(Real r)
    {
        for (SizeT i = 0; i < m_fixed_size; ++i)
        {
            m_data[i] -= r;
        }
    }

    void VectorN::Multiply(Real r)
    {
        for (SizeT i = 0; i < m_fixed_size; ++i)
        {
            m_data[i] *= r;
        }
    }

    void VectorN::Divide(Real r)
    {
        for (SizeT i = 0; i < m_fixed_size; ++i)
        {
            m_data[i] /= r;
        }
    }

    bool VectorN::IsEqual(const VectorN& v) const
    {
        if (this->m_fixed_size != v.m_fixed_size)
            return false;

        for (SizeT i = 0; i < m_fixed_size; ++i)
        {
            if (!Tools::IsEqual(m_data[i], v.m_data[i]))
                return false;
        }

        return true;
    }

    bool VectorN::IsNotEqual(const VectorN& v) const
    {
        if (this->m_fixed_size != v.m_fixed_size)
            return true;

        for (SizeT i = 0; i < m_fixed_size; ++i)
        {
            if (!Tools::IsEqual(m_data[i], v.m_data[i]))
                return true;
        }

        return false;
    }

    VectorN VectorN::Negated() const
    {
        VectorN negated(m_fixed_size);
        for (SizeT i = 0; i < m_fixed_size; ++i)
        {
            negated.m_data[i] = -m_data[i];
        }

        return negated;
    }

    VectorN VectorN::Invert() const
    {
        VectorN invert(m_fixed_size);
        for (SizeT i = 0; i < m_fixed_size; ++i)
        {
            invert.m_data[i] = Tools::Inverse(m_data[i]);
        }

        return invert;
    }

    VectorN VectorN::Added(const VectorN& v) const
    {
        VectorN added(m_fixed_size);
        if (m_fixed_size == v.m_fixed_size)
        {
            for (SizeT i = 0; i < m_fixed_size; ++i)
            {
                added.m_data[i] = m_data[i] + v.m_data[i];
            }
        }

        return added;
    }

    VectorN VectorN::Subtracted(const VectorN& v) const
    {
        VectorN subtracted(m_fixed_size);
        if (m_fixed_size == v.m_fixed_size)
        {
            for (SizeT i = 0; i < m_fixed_size; ++i)
            {
                subtracted.m_data[i] = m_data[i] - v.m_data[i];
            }
        }

        return subtracted;
    }

    VectorN VectorN::Multiplied(Real r) const
    {
        VectorN added(m_fixed_size);
        for (SizeT i = 0; i < m_fixed_size; ++i)
        {
            added.m_data[i] = m_data[i] * r;
        }

        return added;
    }

    VectorN VectorN::Divided(Real r) const
    {
        VectorN added(m_fixed_size);
        for (SizeT i = 0; i < m_fixed_size; ++i)
        {
            added.m_data[i] = m_data[i] / r;
        }

        return added;
    }

    SizeT VectorN::N() const
    {
        return m_fixed_size;
    }

    Real VectorN::PNorm(Real p) const
    {
        // compute p-norm
        Real sum = Constant::ZERO;
        for (SizeT i = 0; i < m_fixed_size; ++i)
        {
            sum += Tools::Pow(Tools::Abs(m_data[i]), p);
        }

        return Tools::Pow(sum, Constant::ONE / p);
    }

    Real VectorN::EuclideanNorm() const
    {
        Real sum = Constant::ZERO;
        for (SizeT i = 0; i < m_fixed_size; ++i)
        {
            sum += m_data[i] * m_data[i];
        }

        return Tools::Sqrt(sum);
    }

    VectorN VectorN::ComputeConvexCombination2(const VectorN& a, const VectorN& b, Real u)
    {
        VectorN computed(a.m_fixed_size);

        if (a.m_fixed_size == b.m_fixed_size)
        {
            Real v = Constant::ONE - u;
            for (SizeT i = 0; i < computed.m_fixed_size; ++i)
            {
                computed.m_data[i] = v * a.m_data[i] + u * b.m_data[i];
            }
        }

        return computed;
    }

    VectorN VectorN::ComputeConvexCombination3(const VectorN& a, const VectorN& b, const VectorN& c, Real u, Real v)
    {
        VectorN computed(a.m_fixed_size);

        if (a.m_fixed_size == b.m_fixed_size && b.m_fixed_size == c.m_fixed_size)
        {
            Real w = Constant::ONE - u - v;
            for (SizeT i = 0; i < computed.m_fixed_size; ++i)
            {
                computed.m_data[i] = u * a.m_data[i] + v * b.m_data[i] + w * c.m_data[i];
            }
        }

        return computed;
    }

    Real VectorN::DotProduct(const VectorN& a, const VectorN& b)
    {
        if (a.m_fixed_size != b.m_fixed_size)
            return Constant::RNAN;

        Real sum = Constant::ZERO;
        for (SizeT i = 0; i < a.m_fixed_size; ++i)
        {
            sum += a.m_data[i] * b.m_data[i];
        }

        return sum;
    }

    VectorN VectorN::Project(const VectorN& a, const VectorN& b)
    {
        return (DotProduct(a, b) / DotProduct(b, b)) * b;
    }

    VectorN operator+(const VectorN& a, const VectorN& b)
    {
        return a.Added(b);
    }

    VectorN operator-(const VectorN& a, const VectorN& b)
    {
        return a.Subtracted(b);
    }

    VectorN operator*(Real real, const VectorN& vector)
    {
        return vector.Multiplied(real);
    }

    VectorN operator*(const VectorN& vector, Real real)
    {
        return vector.Multiplied(real);
    }

    VectorN operator/(const VectorN& vector, Real real)
    {
        return vector.Divided(real);
    }

    std::ostream& operator<<(std::ostream& os, const VectorN& rhs)
    {
        SizeT size = rhs.N();

        os << "[";
        for (SizeT i = 0; i < size - 1; ++i)
        {
            os << rhs[i] << ", ";
        }
        os << rhs[size - 1] << "]";

        return os;
    }
}
