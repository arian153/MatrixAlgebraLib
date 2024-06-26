#pragma once
#include <array>
#include <ostream>
#include <vector>

#include "MatrixLib.Math.Constants.h"
#include "MatrixLib.Math.Tools.h"

namespace MatrixLib
{
    template <SizeT N>
    class VectorN
    {
    private: // private data
        std::array<Real, N> m_elements;

    public: // constructors
        VectorN();
        VectorN(std::initializer_list<Real> values);
        VectorN(const VectorN& rhs);
        ~VectorN() = default;
        explicit VectorN(Real scalar);
        explicit VectorN(const std::vector<Real>& vec);
        explicit VectorN(const std::array<Real, N>& arr);
        explicit VectorN(Real arr[N]);

    public: // operators
        // explicit operator
        explicit operator std::vector<Real>() const;

        // access operator
        Real  operator[](SizeT i) const;
        Real& operator[](SizeT i);
        Real  operator()(SizeT i) const;
        Real& operator()(SizeT i);

        // assignment operator
        VectorN& operator =(const VectorN& rhs);
        VectorN& operator =(Real rhs);

        // compare operator
        bool operator ==(const VectorN& rhs) const;
        bool operator !=(const VectorN& rhs) const;

        // arithmetic operator
        VectorN  operator -() const;
        VectorN& operator +=(const VectorN& rhs);
        VectorN& operator -=(const VectorN& rhs);
        VectorN& operator *=(Real r);
        VectorN& operator /=(Real r);

    public: // modify methods
        void Set(const VectorN& v);
        void SetScalar(Real value);
        void SetZero();
        void Negate();
        void Inverse();
        void Normalize();
        void ClearDigit(SizeT digit);
        void ClearError(Real epsilon);

    public: // info methods
        bool IsEqual(const VectorN& v) const;
        bool IsNotEqual(const VectorN& v) const;
        bool IsZero() const;
        bool IsNormalized() const;

        VectorN Negated() const;
        VectorN Invert() const;
        VectorN Normalized() const;
        VectorN Scaled(Real r) const;
        Real    Length() const;
        Real    LengthSq() const;
        SizeT   SmallestIdx() const;
        SizeT   LargestIdx() const;
        VectorN Half() const;
        VectorN Absolute() const;
        Real    Norm(Real p) const;

    public: // template methods
        template <SizeT M>
        VectorN<M> Swizzle(std::initializer_list<SizeT> indices_list) const;

    private: // static methods
        static Real       DotProductImpl(const VectorN& a, const VectorN& b);
        static VectorN    HadamardProductImpl(const VectorN& a, const VectorN& b);
        static VectorN<3> CrossProductImpl(const VectorN<3>& a, const VectorN<3>& b);
        static VectorN<7> CrossProductImpl(const VectorN<7>& a, const VectorN<7>& b);

    public: // range based for loop related methods
        auto begin();
        auto end();
        auto begin() const;
        auto end() const;

    public: // friend    
        friend class Vector;
        friend class Matrix;
    };

    template <SizeT N>
    VectorN<N>::VectorN()
    {
        m_elements.fill(Constant::ZERO);
    }

    template <SizeT N>
    VectorN<N>::VectorN(std::initializer_list<Real> values)
    {
        m_elements.fill(Constant::ZERO);

        SizeT valid_compo_count = Tools::Min(N, values.size());
        for (SizeT i = 0; i < valid_compo_count; ++i)
        {
            m_elements[i] = *(values.begin() + i);
        }
    }

    template <SizeT N>
    VectorN<N>::VectorN(const VectorN& rhs)
        : m_elements(rhs.m_elements)
    {
    }

    template <SizeT N>
    VectorN<N>::VectorN(Real scalar)
    {
        m_elements.fill(scalar);
    }

    template <SizeT N>
    VectorN<N>::VectorN(const std::vector<Real>& vec)
    {
        m_elements.fill(Constant::ZERO);

        SizeT valid_compo_count = Tools::Min(N, vec.size());
        for (SizeT i = 0; i < valid_compo_count; ++i)
        {
            m_elements[i] = vec[i];
        }
    }

    template <SizeT N>
    VectorN<N>::VectorN(const std::array<Real, N>& arr)
        : m_elements(arr)
    {
    }

    template <SizeT N>
    VectorN<N>::VectorN(Real arr[N])
    {
        for (SizeT i = 0; i < N; ++i)
        {
            m_elements[i] = arr[i];
        }
    }

    template <SizeT N>
    VectorN<N>::operator std::vector<Real>() const
    {
        return std::vector<Real>(m_elements);
    }

    template <SizeT N>
    Real VectorN<N>::operator[](SizeT i) const
    {
        return m_elements[i];
    }

    template <SizeT N>
    Real& VectorN<N>::operator[](SizeT i)
    {
        return m_elements[i];
    }

    template <SizeT N>
    Real VectorN<N>::operator()(SizeT i) const
    {
        return m_elements[i];
    }

    template <SizeT N>
    Real& VectorN<N>::operator()(SizeT i)
    {
        return m_elements[i];
    }

    template <SizeT N>
    VectorN<N>& VectorN<N>::operator=(const VectorN& rhs)
    {
        if (this != &rhs)
            Set(rhs);

        return *this;
    }

    template <SizeT N>
    VectorN<N>& VectorN<N>::operator=(Real rhs)
    {
        SetScalar(rhs);
        return *this;
    }

    template <SizeT N>
    bool VectorN<N>::operator==(const VectorN& rhs) const
    {
        return IsEqual(rhs);
    }

    template <SizeT N>
    bool VectorN<N>::operator!=(const VectorN& rhs) const
    {
        return IsNotEqual(rhs);
    }

    template <SizeT N>
    VectorN<N> VectorN<N>::operator-() const
    {
        return Negated();
    }

    template <SizeT N>
    VectorN<N>& VectorN<N>::operator+=(const VectorN& rhs)
    {
        for (SizeT i = 0; i < N; ++i)
        {
            m_elements[i] += rhs.m_elements[i];
        }

        return *this;
    }

    template <SizeT N>
    VectorN<N>& VectorN<N>::operator-=(const VectorN& rhs)
    {
        for (SizeT i = 0; i < N; ++i)
        {
            m_elements[i] -= rhs.m_elements[i];
        }

        return *this;
    }

    template <SizeT N>
    VectorN<N>& VectorN<N>::operator*=(Real r)
    {
        for (SizeT i = 0; i < N; ++i)
        {
            m_elements[i] *= r;
        }

        return *this;
    }

    template <SizeT N>
    VectorN<N>& VectorN<N>::operator/=(Real r)
    {
        for (SizeT i = 0; i < N; ++i)
        {
            m_elements[i] /= r;
        }

        return *this;
    }

    template <SizeT N>
    void VectorN<N>::Set(const VectorN& v)
    {
        m_elements = v.m_elements;
    }

    template <SizeT N>
    void VectorN<N>::SetScalar(Real value)
    {
        m_elements.fill(value);
    }

    template <SizeT N>
    void VectorN<N>::SetZero()
    {
        m_elements.fill(Constant::ZERO);
    }

    template <SizeT N>
    void VectorN<N>::Negate()
    {
        for (SizeT i = 0; i < N; ++i)
        {
            m_elements[i] = -m_elements[i];
        }
    }

    template <SizeT N>
    void VectorN<N>::Inverse()
    {
        for (SizeT i = 0; i < N; ++i)
        {
            m_elements[i] = Tools::Inverse(m_elements[i]);
        }
    }

    template <SizeT N>
    void VectorN<N>::Normalize()
    {
        Real inv_length = Tools::Inverse(Length());
        for (SizeT i = 0; i < N; ++i)
        {
            m_elements[i] *= inv_length;
        }
    }

    template <SizeT N>
    void VectorN<N>::ClearDigit(SizeT digit)
    {
        for (SizeT i = 0; i < N; ++i)
        {
            m_elements[i] = Tools::ClearDigit(m_elements[i], digit);
            if (Tools::IsZero(m_elements[i]))
                m_elements[i] = Constant::ZERO;
        }
    }

    template <SizeT N>
    void VectorN<N>::ClearError(Real epsilon)
    {
        for (SizeT i = 0; i < N; ++i)
        {
            m_elements[i] = Tools::ClearError(m_elements[i], epsilon);
            if (Tools::IsZero(m_elements[i]))
                m_elements[i] = Constant::ZERO;
        }
    }

    template <SizeT N>
    bool VectorN<N>::IsEqual(const VectorN& v) const
    {
        for (SizeT i = 0; i < N; ++i)
        {
            if (!Tools::IsEqual(m_elements[i], v.m_elements[i]))
                return false;
        }

        return true;
    }

    template <SizeT N>
    bool VectorN<N>::IsNotEqual(const VectorN& v) const
    {
        for (SizeT i = 0; i < N; ++i)
        {
            if (!Tools::IsEqual(m_elements[i], v.m_elements[i]))
                return true;
        }

        return false;
    }

    template <SizeT N>
    bool VectorN<N>::IsZero() const
    {
        for (SizeT i = 0; i < N; ++i)
        {
            if (!Tools::IsZero(m_elements[i]))
                return false;
        }

        return true;
    }

    template <SizeT N>
    bool VectorN<N>::IsNormalized() const
    {
        return Tools::IsEqual(Length(), Constant::ONE);
    }

    template <SizeT N>
    VectorN<N> VectorN<N>::Negated() const
    {
        VectorN negated;
        for (SizeT i = 0; i < N; ++i)
        {
            negated.m_elements[i] = -m_elements[i];
        }

        return negated;
    }

    template <SizeT N>
    VectorN<N> VectorN<N>::Invert() const
    {
        VectorN invert;
        for (SizeT i = 0; i < N; ++i)
        {
            invert.m_elements[i] = Tools::Inverse(m_elements[i]);
        }

        return invert;
    }

    template <SizeT N>
    VectorN<N> VectorN<N>::Normalized() const
    {
        VectorN normalized;
        Real    inv_length = Tools::Inverse(Length());
        for (SizeT i = 0; i < N; ++i)
        {
            normalized.m_elements[i] = m_elements[i] * inv_length;
        }

        return normalized;
    }

    template <SizeT N>
    VectorN<N> VectorN<N>::Scaled(Real r) const
    {
        VectorN multiplied;
        for (SizeT i = 0; i < N; ++i)
        {
            multiplied.m_elements[i] = m_elements[i] * r;
        }

        return multiplied;
    }

    template <SizeT N>
    Real VectorN<N>::Length() const
    {
        return Tools::Sqrt(LengthSq());
    }

    template <SizeT N>
    Real VectorN<N>::LengthSq() const
    {
        Real sum = Constant::ZERO;
        for (SizeT i = 0; i < N; ++i)
        {
            sum += m_elements[i] * m_elements[i];
        }

        return sum;
    }

    template <SizeT N>
    SizeT VectorN<N>::SmallestIdx() const
    {
        SizeT smallest_idx = 0;
        for (SizeT i = 1; i < N; ++i)
        {
            if (m_elements[i] < m_elements[smallest_idx])
                smallest_idx = i;
        }

        return smallest_idx;
    }

    template <SizeT N>
    SizeT VectorN<N>::LargestIdx() const
    {
        SizeT largest_idx = 0;
        for (SizeT i = 1; i < N; ++i)
        {
            if (m_elements[i] > m_elements[largest_idx])
                largest_idx = i;
        }

        return largest_idx;
    }

    template <SizeT N>
    VectorN<N> VectorN<N>::Half() const
    {
        VectorN half;
        for (SizeT i = 0; i < N; ++i)
        {
            half.m_elements[i] = static_cast<Real>(0.5) * m_elements[i];
        }

        return half;
    }

    template <SizeT N>
    VectorN<N> VectorN<N>::Absolute() const
    {
        VectorN absolute;
        for (SizeT i = 0; i < N; ++i)
        {
            absolute.m_elements[i] = Tools::Abs(m_elements[i]);
        }

        return absolute;
    }

    template <SizeT N>
    Real VectorN<N>::Norm(Real p) const
    {
        Real sum = Constant::ZERO;
        for (SizeT i = 0; i < N; ++i)
        {
            sum += Tools::Pow(Tools::Abs(m_elements[i]), p);
        }

        return Tools::Pow(sum, Tools::Inverse(p));
    }

    template <SizeT N>
    template <SizeT M>
    VectorN<M> VectorN<N>::Swizzle(std::initializer_list<SizeT> indices_list) const
    {
        VectorN<M> selected;
        SizeT      valid_compo_count = Tools::Min(M, indices_list.size());
        for (SizeT i = 0; i < valid_compo_count; ++i)
        {
            selected[i] = m_elements[*(indices_list.begin() + i)];
        }

        return selected;
    }

    template <SizeT N>
    Real VectorN<N>::DotProductImpl(const VectorN& a, const VectorN& b)
    {
        Real dot = Constant::ZERO;
        for (SizeT i = 0; i < N; ++i)
        {
            dot += a.m_elements[i] * b.m_elements[i];
        }

        return dot;
    }

    template <SizeT N>
    VectorN<N> VectorN<N>::HadamardProductImpl(const VectorN& a, const VectorN& b)
    {
        VectorN hadamard;
        for (SizeT i = 0; i < N; ++i)
        {
            hadamard = a.m_elements[i] * b.m_elements[i];
        }

        return hadamard;
    }

    inline VectorN<3> VectorN<3>::CrossProductImpl(const VectorN& a, const VectorN& b)
    {
        VectorN cross;
        cross[0] = a[1] * b[2] - a[2] * b[1];
        cross[1] = a[2] * b[0] - a[0] * b[2];
        cross[2] = a[0] * b[1] - a[1] * b[0];

        return cross;
    }

    inline VectorN<7> VectorN<7>::CrossProductImpl(const VectorN& a, const VectorN& b)
    {
        VectorN cross;
        for (SizeT i = 0; i < 7; ++i)
        {
            SizeT e1 = i;
            SizeT e2 = (i + 1) % 7;
            SizeT e3 = (i + 2) % 7;
            SizeT e4 = (i + 3) % 7;
            SizeT e5 = (i + 4) % 7;
            SizeT e6 = (i + 5) % 7;
            SizeT e7 = (i + 6) % 7;

            cross[e1] = a[e2] * b[e4] - a[e4] * b[e2]
                    + a[e3] * b[e7] - a[e7] * b[e3]
                    + a[e5] * b[e6] - a[e6] * b[e5];
        }

        return cross;
    }

    // range based for
    template <SizeT N>
    auto VectorN<N>::begin()
    {
        return m_elements.begin();
    }

    template <SizeT N>
    auto VectorN<N>::end()
    {
        return m_elements.end();
    }

    template <SizeT N>
    auto VectorN<N>::begin() const
    {
        return m_elements.begin();
    }

    template <SizeT N>
    auto VectorN<N>::end() const
    {
        return m_elements.end();
    }
}
