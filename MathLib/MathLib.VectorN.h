#pragma once
#include <ostream>
#include <vector>

#include "MathLib.Constants.h"
#include "MathLib.Tools.h"

namespace MathLib
{
    class VectorN
    {
    private: // private data
        std::vector<Real> m_elements;

    public: // constructors
        VectorN();
        VectorN(std::initializer_list<Real> values);
        VectorN(const VectorN& rhs);
        ~VectorN() = default;
        explicit VectorN(Real scalar);
        explicit VectorN(const std::vector<Real>& vec);

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
        VectorN Swizzle(std::initializer_list<SizeT> indices_list) const;

    public: // static methods
        static Real    DotProduct(const VectorN& a, const VectorN& b);
        static VectorN HadamardProduct(const VectorN& a, const VectorN& b);
        static VectorN CrossProduct3D(const VectorN& a, const VectorN& b);
        static VectorN CrossProduct7D(const VectorN& a, const VectorN& b);

    public: // range based for loop related methods
        auto begin();
        auto end();
        auto begin() const;
        auto end() const;
    };
}
