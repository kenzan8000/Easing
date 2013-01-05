#include "Easing.h"
#include <math.h>


/* ***** Linear interpolation (no easing) ***** */

/// Modeled after the line: y = x
float kLinearInterpolation(float p) { return p; }


/* ***** Quadratic easing p^2 ***** */

/// Modeled after the parabola: y = x^2
float kQuadraticEaseIn(float p) { return p * p; }

/// Modeled after the parabola: y = -x^2 + 2x
float kQuadraticEaseOut(float p) { return -(p * (p - 2)); }

/// Modeled after the piecewise quadratic: y = (1/2)((2x)^2) [0, 0.5] | y = -(1/2)((2x-1)*(2x-3) - 1) [0.5, 1]
float kQuadraticEaseInOut(float p) { return (p < 0.5) ? 2 * p * p : (-2 * p * p) + (4 * p) - 1; }


/* ***** Cubic easing p^3 ***** */

/// Modeled after the cubic: y = x^3
float kCubicEaseIn(float p) { return p * p * p; }

/// Modeled after the cubic: y = (x - 1)^3 + 1
float kCubicEaseOut(float p) { return (p - 1) * (p - 1) * (p - 1) + 1; }

/// Modeled after the piecewise cubic: y = (1/2)((2x)^3) [0, 0.5) | y = (1/2)((2x-2)^3 + 2) [0.5, 1]
float kCubicEaseInOut(float p) { return (p < 0.5) ? 4 * p * p * p : 0.5 * ((2 * p) - 2) * ((2 * p) - 2) * ((2 * p) - 2) + 1; }


/* ***** Quartic easing p^4 ***** */

/// Modeled after the quartic: y = x^4
float kQuarticEaseIn(float p) { return p * p * p * p; }

/// Modeled after the quartic: y = 1 - (x - 1)^4
float kQuarticEaseOut(float p) { return (p - 1) * (p - 1) * (p - 1) * (1 - p) + 1; }

/// Modeled after the piecewise quartic: y = (1/2)((2x)^4) [0, 0.5] | y = -(1/2)((2x-2)^4 - 2) [0.5, 1]
float kQuarticEaseInOut(float p) { return (p < 0.5) ? 8 * p * p * p * p : -8 * (p - 1) * (p - 1) * (p - 1) * (p - 1) + 1; }


/* ***** Quintic easing p^5 ***** */

/// Modeled after the quintic: y = x^5
float kQuinticEaseIn(float p) { return p * p * p * p * p; }

/// Modeled after the quintic: y = (x - 1)^5 + 1
float kQuinticEaseOut(float p) { return (p - 1) * (p - 1) * (p - 1) * (p - 1) * (p - 1) + 1; }

/// Modeled after the piecewise quintic: y = (1/2)((2x)^5) [0, 0.5] | y = (1/2)((2x-2)^5 + 2) [0.5, 1]
float kQuinticEaseInOut(float p) { return (p < 0.5) ? 16 * p * p * p * p * p : 0.5 * ((2 * p) - 2) * ((2 * p) - 2) * ((2 * p) - 2) * ((2 * p) - 2) * ((2 * p) - 2) + 1; }


/* ***** Sine wave easing sin(p * PI/2) ***** */

/// Modeled after quarter-cycle of sine wave
float kSineEaseIn(float p) { return sin((p - 1) * M_PI_2) + 1; }

/// Modeled after quarter-cycle of sine wave (different phase)
float kSineEaseOut(float p) { return sin(p * M_PI_2); }

/// Modeled after half sine wave
float kSineEaseInOut(float p) { return 0.5 * (1 - cos(p * M_PI)); }


/* ***** Circular easing sqrt(1 - p^2) ***** */

/// Modeled after shifted quadrant IV of unit circle
float kCircularEaseIn(float p) { return 1 - sqrt(1 - (p * p)); }

/// Modeled after shifted quadrant II of unit circle
float kCircularEaseOut(float p) { return sqrt((2 - p) * p); }

/// Modeled after the piecewise circular function: y = (1/2)(1 - sqrt(1 - 4x^2)) [0, 0.5] | y = (1/2)(sqrt(-(2x - 3)*(2x - 1)) + 1) [0.5, 1]
float kCircularEaseInOut(float p) { return (p < 0.5) ? 0.5 * (1 - sqrt(1 - 4 * (p * p))) : 0.5 * (sqrt(-((2 * p) - 3) * ((2 * p) - 1)) + 1); }


/* ***** Exponential easing, base 2 ***** */

/// Modeled after the exponential function: y = 2^(10(x - 1))
float kExponentialEaseIn(float p) { return (p == 0.0) ? p : pow(2, 10 * (p - 1)); }

/// Modeled after the exponential function: y = -2^(-10x) + 1
float kExponentialEaseOut(float p) { return (p == 1.0) ? p : 1 - pow(2, -10 * p); }

/// Modeled after the piecewise exponential: y = (1/2)2^(10(2x - 1)) [0,0.5] | y = -(1/2)*2^(-10(2x - 1))) + 1 [0.5,1]
float kExponentialEaseInOut(float p) { return (p == 0.0 || p == 1.0) ? p : (p < 0.5) ? 0.5 * pow(2, (20 * p) - 10) : -0.5 * pow(2, (-20 * p) + 10) + 1; }


/* ***** Exponentially-damped sine wave easing ***** */

/// Modeled after the damped sine wave: y = sin(13pi/2*x)*pow(2, 10 * (x - 1))
float kElasticEaseIn(float p) { return sin(13 * M_PI_2 * p) * pow(2, 10 * (p - 1)); }

/// Modeled after the damped sine wave: y = sin(-13pi/2*(x + 1))*pow(2, -10x) + 1
float kElasticEaseOut(float p) { return sin(-13 * M_PI_2 * (p + 1)) * pow(2, -10 * p) + 1; }

/// Modeled after the piecewise exponentially-damped sine wave: y = (1/2)*sin(13pi/2*(2*x))*pow(2, 10 * ((2*x) - 1)) [0,0.5] | y = (1/2)*(sin(-13pi/2*((2x-1)+1))*pow(2,-10(2*x-1)) + 2) [0.5, 1]
float kElasticEaseInOut(float p) { return (p < 0.5) ? 0.5 * sin(13 * M_PI_2 * (2 * p)) * pow(2, 10 * ((2 * p) - 1)) : 0.5 * (sin(-13 * M_PI_2 * ((2 * p - 1) + 1)) * pow(2, -10 * (2 * p - 1)) + 2); }


/* ***** Overshooting cubic easing ***** */
/// Modeled after the overshooting cubic: y = x^3-x*sin(x*pi)
float kBackEaseIn(float p) { return p * p * p - p * sin(p * M_PI); }

/// Modeled after overshooting cubic: y = 1-((1-x)^3-(1-x)*sin((1-x)*pi))
float kBackEaseOut(float p) { return 1 - ((1 - p) * (1 - p) * (1 - p) - (1 - p) * sin((1 - p) * M_PI)); }

/// Modeled after the piecewise overshooting cubic function: y = (1/2)*((2x)^3-(2x)*sin(2*x*pi)) [0, 0.5] | y = (1/2)*(1-((1-x)^3-(1-x)*sin((1-x)*pi))+1) [0.5, 1]
float kBackEaseInOut(float p) { return (p < 0.5) ? 0.5 * (8 * p * p * p - (2 * p) * sin((2 * p) * M_PI)) : 0.5 * (1 - ((1 - (2*p - 1)) * (1 - (2*p - 1)) * (1 - (2*p - 1)) - (1 - (2*p - 1)) * sin((1 - (2*p - 1)) * M_PI))) + 0.5; }


/* ***** Exponentially-decaying bounce easing ***** */
float kBounceEaseOut(float p);

///
float kBounceEaseIn(float p) { return 1 - kBounceEaseOut(1 - p); }

///
float kBounceEaseOut(float p) { return (p < 4/11.0) ? (121 * p * p)/16.0 : (p < 8/11.0) ? (363/40.0 * p * p) - (99/10.0 * p) + 17/5.0 : (p < 9/10.0) ? (4356/361.0 * p * p) - (35442/1805.0 * p) + 16061/1805.0 : (54/5.0 * p * p) - (513/25.0 * p) + 268/25.0; }

///
float kBounceEaseInOut(float p) { return (p < 0.5) ? 0.5 * kBounceEaseIn(p*2) : 0.5 * kBounceEaseOut(p * 2 - 1) + 0.5; }
