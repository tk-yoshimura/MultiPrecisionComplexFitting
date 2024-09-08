using MultiPrecision;
using MultiPrecisionComplex;
using MultiPrecisionComplexAlgebra;
using System;
using System.Collections.Generic;
using System.Linq;

namespace MultiPrecisionComplexFitting {
    public static class ComplexCurveFittingUtils {
        public static MultiPrecision<N> RelativeError<N>(Complex<N> expected, Complex<N> actual) where N : struct, IConstant {
            if (Complex<N>.IsZero(expected)) {
                return Complex<N>.IsZero(actual) ? MultiPrecision<N>.Zero : MultiPrecision<N>.PositiveInfinity;
            }
            else {
                Complex<N> diff = expected - actual;

                MultiPrecision<N> error = MultiPrecision<N>.Max(MultiPrecision<N>.Abs(diff.R), MultiPrecision<N>.Abs(diff.I)) / expected.Magnitude;

                return error;
            }
        }

        public static MultiPrecision<N> AbsoluteError<N>(Complex<N> expected, Complex<N> actual) where N : struct, IConstant {
            Complex<N> diff = expected - actual;

            MultiPrecision<N> error = MultiPrecision<N>.Max(MultiPrecision<N>.Abs(diff.R), MultiPrecision<N>.Abs(diff.I));

            return error;
        }

        public static MultiPrecision<N> MaxRelativeError<N>(ComplexVector<N> expected, ComplexVector<N> actual) where N : struct, IConstant {
            if (expected.Dim != actual.Dim) {
                throw new ArgumentException("mismatch size", $"{nameof(expected)},{nameof(actual)}");
            }

            MultiPrecision<N> max_error = MultiPrecision<N>.Zero;

            for (int i = 0, n = expected.Dim; i < n; i++) {
                max_error = MultiPrecision<N>.Max(max_error, RelativeError(expected[i], actual[i]));
            }

            return max_error;
        }

        public static MultiPrecision<N> MaxAbsoluteError<N>(ComplexVector<N> expected, ComplexVector<N> actual) where N : struct, IConstant {
            if (expected.Dim != actual.Dim) {
                throw new ArgumentException("mismatch size", $"{nameof(expected)},{nameof(actual)}");
            }

            MultiPrecision<N> max_error = MultiPrecision<N>.Zero;

            for (int i = 0, n = expected.Dim; i < n; i++) {
                max_error = MultiPrecision<N>.Max(max_error, AbsoluteError(expected[i], actual[i]));
            }

            return max_error;
        }

        public static IEnumerable<(Complex<N> numer, Complex<N> denom)> EnumeratePadeCoef<N>(ComplexVector<N> param, int m, int n) where N : struct, IConstant {
            if (param.Dim != checked(m + n)) {
                throw new ArgumentException("invalid param dims", nameof(param));
            }

            if (m >= n) {
                for (int i = 0; i < n; i++) {
                    yield return (param[..m][i], param[m..][i]);
                }
                for (int i = 0; i < m - n; i++) {
                    yield return (param[..m][n + i], 0);
                }
            }
            else {
                for (int i = 0; i < m; i++) {
                    yield return (param[..m][i], param[m..][i]);
                }
                for (int i = 0; i < n - m; i++) {
                    yield return (0, param[m..][m + i]);
                }
            }
        }

        public static (long exp_scale, ComplexVector<N> v_standardized) StandardizeExponent<N>(ComplexVector<N> v) where N : struct, IConstant {
            if (v.Dim <= 0 || v.All(v => Complex<N>.IsZero(v.val))) {
                throw new ArgumentException("invalid vector because it zero vector", nameof(v));
            }

            if (v.Any(v => !Complex<N>.IsFinite(v.val))) {
                throw new ArgumentException("invalid vector because it contains infinite values", nameof(v));
            }

            long exp_scale = v.Max(v => long.Max(v.val.R.Exponent, v.val.I.Exponent));
            ComplexVector<N> v_standardized = (x => (MultiPrecision<N>.Ldexp(x.R, -exp_scale), MultiPrecision<N>.Ldexp(x.I, -exp_scale)), v);

            return (exp_scale, v_standardized);
        }
    }
}
