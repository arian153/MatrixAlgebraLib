#pragma once
#include <limits>

namespace MathLib
{
    using Real32 = float;
    using Real64 = double;
    using Int8 = __int8;
    using Int16 = __int16;
    using Int32 = __int32;
    using Int64 = __int64;
    using UInt8 = unsigned __int8;
    using UInt16 = unsigned __int16;
    using UInt32 = unsigned __int32;
    using UInt64 = unsigned __int64;

    // 32bit
    using Real = Real32;
    using UInt = UInt32;

    // size of stl data
    using SizeT = size_t;

    namespace Constant
    {
        constexpr Real PI         = static_cast<Real>(3.14159265358979323846264);
        constexpr Real HALF_PI    = PI * static_cast<Real>(0.5);
        constexpr Real QUARTER_PI = PI * static_cast<Real>(0.25);
        constexpr Real TWO_PI     = PI * static_cast<Real>(2.0);
        constexpr Real PI_DIV_2   = HALF_PI;
        constexpr Real PI_DIV_3   = PI / static_cast<Real>(3.0);
        constexpr Real PI_DIV_4   = QUARTER_PI;
        constexpr Real PI_DIV_6   = PI / static_cast<Real>(6.0);
        constexpr Real RADIAN     = PI / static_cast<Real>(180.0);
        constexpr Real ZERO       = static_cast<Real>(0.0);
        constexpr Real ONE        = static_cast<Real>(1.0);

        constexpr Real EPSILON         = static_cast<Real>(0.00001);
        constexpr Real EPSILON_SQUARED = EPSILON * EPSILON;
        constexpr Real EPSILON_BIAS    = static_cast<Real>(1.00001);

        constexpr Real REAL_MAX = std::numeric_limits<Real>::max();
        constexpr Real REAL_MIN = std::numeric_limits<Real>::min();

        constexpr Real REAL_POS_MAX = REAL_MAX; //go to +max
        constexpr Real REAL_POS_MIN = REAL_MIN; //near to 0
        constexpr Real REAL_NEG_MAX = -REAL_MAX; //go to -max
        constexpr Real REAL_NEG_MIN = -REAL_MIN; //near to 0

        constexpr Real ROOT_TWO        = static_cast<Real>(1.41421356237309504880168);
        constexpr Real ROOT_THREE      = static_cast<Real>(1.73205080756887729352744);
        constexpr Real ROOT_FIVE       = static_cast<Real>(2.23606797749978969640917);
        constexpr Real ROOT_TEN        = static_cast<Real>(3.16227766016837933199889);
        constexpr Real CUBE_ROOT_TWO   = static_cast<Real>(1.25992104989487316476721);
        constexpr Real CUBE_ROOT_THREE = static_cast<Real>(1.25992104989487316476721);
        constexpr Real FORTH_ROOT_TWO  = static_cast<Real>(1.18920711500272106671749);

        constexpr Real E_MATHEMATICAL = static_cast<Real>(2.71828182845904523536);
        constexpr Real LN_TWO         = static_cast<Real>(0.69314718055994530941723);
        constexpr Real LN_THREE       = static_cast<Real>(1.09861228866810969139524);
        constexpr Real LN_TEN         = static_cast<Real>(2.30258509299404568401799);

        const Real RNAN = std::numeric_limits<Real>::quiet_NaN();
        const Real INF  = std::numeric_limits<Real>::infinity();
    }}
