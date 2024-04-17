#pragma once
#include <ostream>
#include <vector>

#include "MatrixLib.Math.Constants.h"

namespace MatrixLib
{
    class VectorN
    {
    public: // defaults    
        VectorN(const VectorN& rhs);
        ~VectorN() = default;
        explicit VectorN(const std::vector<Real>& v);
        explicit VectorN(SizeT size_of_vector, Real value = Constant::ZERO);

    public: // operators
        // explicit operator
        explicit operator std::vector<Real>() const;

        // access operator
        Real  operator[](size_t i) const;
        Real& operator[](size_t i);
        Real  operator()(size_t i) const;
        Real& operator()(size_t i);

        // assignment operator
        VectorN& operator =(const VectorN& rhs);
        VectorN& operator =(Real rhs);

        // compare operator
        bool operator ==(const VectorN& rhs) const;
        bool operator !=(const VectorN& rhs) const;

        // arithmetic operator
        VectorN operator -() const;

        VectorN& operator +=(const VectorN& rhs);
        VectorN& operator -=(const VectorN& rhs);
        VectorN& operator +=(Real r);
        VectorN& operator -=(Real r);
        VectorN& operator *=(Real r);
        VectorN& operator /=(Real r);

    public: // modify methods     
        void Set(const VectorN& v);
        void SetScalar(Real value);
        void SetZero();

        void Negate();
        void Inverse();
        void Add(const VectorN& v);
        void Subtract(const VectorN& v);
        void Add(Real r);
        void Subtract(Real r);
        void Multiply(Real r);
        void Divide(Real r);

    public: // info methods

        bool IsEqual(const VectorN& v) const;
        bool IsNotEqual(const VectorN& v) const;

        VectorN Negated() const;
        VectorN Invert() const;
        VectorN Added(const VectorN& v) const;
        VectorN Subtracted(const VectorN& v) const;
        VectorN Multiplied(Real r) const;
        VectorN Divided(Real r) const;

        SizeT N() const;
        Real  PNorm(Real p) const;
        Real  EuclideanNorm() const;

        // TODO
        // Length 2Norm
        // LengthSq
        // Smallest Idx
        // Largest Idx
        // Half
        // Absolute
        // Scale
        // IsZero
        // IsNormalized
        // Normalize
        // Normalized

    public: // static methods
        static VectorN ComputeConvexCombination2(const VectorN& a, const VectorN& b, Real u);
        static VectorN ComputeConvexCombination3(const VectorN& a, const VectorN& b, const VectorN& c, Real u, Real v);
        static Real    DotProduct(const VectorN& a, const VectorN& b);
        static VectorN Project(const VectorN& a, const VectorN& b);

        // TODO
        //static Real    Distance(const VectorN& a, const VectorN& b);
        //static Real    DistanceSq(const VectorN& a, const VectorN& b);
        //static VectorN CrossProduct3D(const VectorN& a, const VectorN& b);
        //static VectorN CrossProduct7D(const VectorN& a, const VectorN& b);
        //static VectorN HadamardProduct(const VectorN& a, const VectorN& b);
        //static MatrixNM OuterProduct(const VectorN& a, const VectorN& b);

    private: // private data
        std::vector<Real> m_data;
        SizeT             m_fixed_size = 0;

    public: // std/stl related methods
        auto begin()
        {
            return m_data.begin();
        }

        auto end()
        {
            return m_data.end();
        }

        auto begin() const
        {
            return m_data.begin();
        }

        auto end() const
        {
            return m_data.end();
        }
    };

    VectorN operator +(const VectorN& a, const VectorN& b);
    VectorN operator -(const VectorN& a, const VectorN& b);
    VectorN operator *(Real real, const VectorN& vector);
    VectorN operator *(const VectorN& vector, Real real);
    VectorN operator /(const VectorN& vector, Real real);
    std::ostream& operator<<(std::ostream& os, const VectorN& rhs);
}
