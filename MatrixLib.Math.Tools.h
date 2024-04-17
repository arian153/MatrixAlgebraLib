#pragma once
#include <cmath>
#include "MatrixLib.Math.Constants.h"

namespace MatrixLib::Tools
{
    inline Real Min(Real a, Real b)
    {
        return a < b ? a : b;
    }

    inline Real Max(Real a, Real b)
    {
        return a > b ? a : b;
    }

    inline Real Clamp(Real x, Real low, Real high)
    {
        return x < low ? low : (x > high ? high : x);
    }

    inline Real Abs(Real rhs)
    {
        return std::abs(rhs);
    }

    inline Real Cos(Real rhs)
    {
        return std::cos(rhs);
    }

    inline Real Sin(Real rhs)
    {
        return std::sin(rhs);
    }

    inline Real Tan(Real rhs)
    {
        return std::tan(rhs);
    }

    inline Real Acos(Real rhs)
    {
        return std::acos(rhs);
    }

    inline Real Asin(Real rhs)
    {
        return std::asin(rhs);
    }

    inline Real Atan(Real x, Real y)
    {
        return std::atan2(y, x);
    }

    inline Real Pow(Real a, Real b)
    {
        return std::pow(a, b);
    }

    // Return t 0.0 ~ 1.0
    inline Real Normalize(Real t, Real min, Real max)
    {
        return (t - min) / (max - min);
    }

    inline bool IsValid(Real rhs)
    {
        return !isnan(rhs) && isfinite(rhs);
    }

    inline bool IsInfinite(Real rhs)
    {
        return !isfinite(rhs);
    }

    inline bool IsZero(Real rhs)
    {
        return (std::abs(rhs) < Constant::EPSILON);
    }

    inline bool IsEqual(Real a, Real b)
    {
        return (std::abs(a - b) < Constant::EPSILON);
    }

    inline bool IsNotEqual(Real a, Real b)
    {
        return std::abs(a - b) >= Constant::EPSILON;
    }

    inline bool IsLess(Real a, Real b)
    {
        return ((a + Constant::EPSILON) < (b - Constant::EPSILON));
    }

    inline bool IsGreator(Real a, Real b)
    {
        return ((a - Constant::EPSILON) > (b + Constant::EPSILON));
    }

    inline bool IsGreatorThanZero(Real a)
    {
        return (a >= Constant::EPSILON);
    }

    inline Real Inverse(Real rhs)
    {
        return IsZero(rhs) ? Constant::ZERO : Constant::ONE / rhs;
    }

    inline Real InvSqrt(Real rhs)
    {
        return IsZero(rhs) ? Constant::ZERO : Constant::ONE / std::sqrt(rhs);
    }

    inline Real Sqrt(Real rhs)
    {
        return std::sqrt(rhs);
    }

    inline Real Cbrt(Real rhs)
    {
        return std::cbrt(rhs);
    }

    inline Real Signum(Real param)
    {
        if (IsZero(param))
        {
            return 0.0;
        }

        return param >= Constant::EPSILON ? Constant::ONE : -Constant::ONE;
    }

    inline int RoundToInt(Real param)
    {
        return static_cast<int>(std::floor(param + static_cast<Real>(0.5)));
    }

    // Convert Degree to Radian
    inline Real DegreeToRadian(Real degrees)
    {
        return (degrees * Constant::RADIAN);
    }

    // Convert Radian to Degree
    inline Real RadianToDegree(Real radians)
    {
        return (radians / Constant::RADIAN);
    }

    // Get angle from x, y coordinates
    inline Real ToAngle(Real x, Real y)
    {
        Real theta = std::atan2(y, x);
        if (theta > Constant::PI)
            theta = theta + Constant::TWO_PI;
        if (theta > Constant::TWO_PI)
            theta -= Constant::TWO_PI;
        if (theta < Constant::ZERO)
            theta += Constant::TWO_PI;
        return theta;
    }

    // Round to nearest value.  
    inline Real ClearError(Real value, Real digit = Constant::EPSILON)
    {
        return std::round(value / digit) * digit;
    }

    inline bool SolveQuadraticEquation(Real a, Real b, Real c, Real& root0, Real& root1)
    {
        if (IsZero(a))
        {
            if (IsZero(b))
            {
                root0 = -c;
                root1 = -c;
            }
            else
            {
                root0 = -c / b;
                root1 = -c / b;
            }

            return false;
        }

        Real discriminant = (b * b) - (static_cast<Real>(4.0) * a * c);

        if (discriminant > Constant::ZERO)
        {
            root0 = (-b + std::sqrt(discriminant)) / (static_cast<Real>(2.0) * a);
            root1 = (-b - std::sqrt(discriminant)) / (static_cast<Real>(2.0) * a);
        }
        else if (IsZero(discriminant) == true)
        {
            Real result = -b / (static_cast<Real>(2.0) * a);
            root0       = result;
            root1       = result;
        }
        else
        {
            Real real      = -b / (static_cast<Real>(2.0) * a);
            Real imaginary = std::sqrt(-discriminant) / (static_cast<Real>(2.0) * a);
            root0          = real;
            root1          = imaginary;
            return false;
        }

        return true;
    }

    inline int SolveCubicEquation(Real a, Real b, Real c, Real& root0, Real& root1, Real& root2)
    {
        Real a2 = a * a;
        Real q  = (a2 - static_cast<Real>(3.0) * b) / static_cast<Real>(9.0);
        Real r  = (a * (static_cast<Real>(2.0) * a2 - static_cast<Real>(9.0) * b) + static_cast<Real>(27.0) * c) / static_cast<Real>(54.0);

        // equation y^3 - 3q*y + r/2 = 0 where x = y-a/3
        if (std::abs(q) < Constant::EPSILON)
        {
            if (std::abs(r) < Constant::EPSILON)
            {	// three identical roots
                root0 = root1 = root2 = -a / static_cast<Real>(3.0);

                return 3;
            }

            // y^3 =-r/2
            root0 = std::cbrt(-r / static_cast<Real>(2.0));
            root1 = root0 * static_cast<Real>(0.5);
            root2 = root0 * std::sqrt(static_cast<Real>(3.0)) / static_cast<Real>(2.0);

            return 1;
        }

        Real r2 = r * r;
        Real q3 = q * q * q;
        if (r2 <= q3 + Constant::EPSILON)
        {
            Real t = r / std::sqrt(q3);
            t      = Clamp(t, static_cast<Real>(-1.0), static_cast<Real>(1.0));

            t = std::acos(t);
            a /= static_cast<Real>(3.0);
            q = static_cast<Real>(-2.0) * std::sqrt(q);

            root0 = q * std::cos(t / static_cast<Real>(3.0)) - a;
            root1 = q * std::cos((t + Constant::TWO_PI) / static_cast<Real>(3.0)) - a;
            root2 = q * std::cos((t - Constant::TWO_PI) / static_cast<Real>(3.0)) - a;

            return 3;
        }

        Real a3 = -std::cbrt(std::abs(r) + std::sqrt(r2 - q3));

        if (r < Constant::ZERO)
            a3 = -a3;

        Real b3 = a3 == Constant::ZERO ? Constant::ZERO : b3 = q / a3;

        a /= static_cast<Real>(3.0);

        root0 = a3 + b3 - a;
        root1 = -static_cast<Real>(0.5) * (a3 + b3) - a;
        root2 = static_cast<Real>(0.5) * std::sqrt(static_cast<Real>(3.0)) * (a3 - b3);

        if (std::abs(root2) < Constant::EPSILON)
        {
            root2 = root1;

            return 2;
        }

        return 1;
    }
}
