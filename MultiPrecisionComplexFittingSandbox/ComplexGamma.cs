using MultiPrecision;
using MultiPrecisionComplex;

namespace MultiPrecisionComplexFittingSandbox {
    public static class ComplexGammaN16 {
        private static readonly int stirling_convergence_norm = 16384;
        private static readonly List<MultiPrecision<Pow2.N16>> coefs = [];

        static ComplexGammaN16() {
            for (int k = 1; k <= 64; k++) {
                MultiPrecision<Pow2.N16> coef = MultiPrecision<Pow2.N16>.BernoulliSequence(k)
                    / checked((2 * k) * (2 * k - 1));

                coefs.Add(coef);
            }
        }

        public static Complex<Pow2.N16> StirlingTerm(Complex<Pow2.N16> z) {
            Complex<Pow2.N16> s = 0, z2 = z * z, w = z;

            for (int i = 0; i < coefs.Count; i++) {
                Complex<Pow2.N16> ds = coefs[i] / w;

                if (s.Magnitude.Exponent > ds.Magnitude.Exponent + MultiPrecision<Pow2.N16>.Bits) {
                    return s;
                }

                w *= z2;
                s += ds;
            }

            return Complex<Pow2.N16>.NaN;
        }

        public static Complex<Pow2.N16> Gamma(Complex<Pow2.N16> z) {
            if (z.R < MultiPrecision<Pow2.N16>.Point5) {
                Complex<Pow2.N16> y = MultiPrecision<Pow2.N16>.PI / (Complex<Pow2.N16>.SinPI(z) * Gamma(1 - z));

                return y;
            }

            static Complex<Pow2.N16> stirling(Complex<Pow2.N16> z) {
                Complex<Pow2.N16> s = StirlingTerm(z);

                Complex<Pow2.N16> y =
                    Complex<Pow2.N16>.Sqrt(2 * MultiPrecision<Pow2.N16>.PI / z) *
                    Complex<Pow2.N16>.Pow(z / MultiPrecision<Pow2.N16>.E, z) *
                    Complex<Pow2.N16>.Exp(s);

                return y;
            }

            if (z.Norm > stirling_convergence_norm) {
                Complex<Pow2.N16> s = StirlingTerm(z);
                Complex<Pow2.N16> y =
                    Complex<Pow2.N16>.Sqrt(2 * MultiPrecision<Pow2.N16>.PI / z) *
                    Complex<Pow2.N16>.Pow(z / MultiPrecision<Pow2.N16>.E, z) *
                    Complex<Pow2.N16>.Exp(s);

                return y;
            }
            else {
                int rn = (int)double.Floor((double)z.R);
                MultiPrecision<Pow2.N16> rf = z.R - rn;

                double zid = (double)z.I, zrd = double.Sqrt(stirling_convergence_norm - zid * zid);
                int rk = (int)double.Ceiling(zrd);

                if (double.IsNaN(zrd) || rn >= rk) {
                    return stirling(z);
                }

                Complex<Pow2.N16> c = stirling((rf + rk, z.I));
                Complex<Pow2.N16> m = z;
                for (int k = rn + 1; k < rk; k++) {
                    m *= (rf + k, z.I);
                }

                Complex<Pow2.N16> y = c / m;

                return y;
            }
        }
    }
}
