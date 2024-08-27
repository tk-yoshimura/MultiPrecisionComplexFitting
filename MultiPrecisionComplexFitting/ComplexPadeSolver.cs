using MultiPrecision;
using MultiPrecisionComplexAlgebra;

namespace MultiPrecisionComplexFitting {
    public static class ComplexPadeSolver<N> where N : struct, IConstant {
        public static (ComplexVector<N> ms, ComplexVector<N> ns) Solve(ComplexVector<N> cs, int m, int n) {
            if (m < 0) {
                throw new ArgumentOutOfRangeException(nameof(m));
            }
            if (n < 0) {
                throw new ArgumentOutOfRangeException(nameof(n));
            }
            if (cs.Dim != checked(m + n + 1)) {
                throw new ArgumentException("invalid size", nameof(cs));
            }

            int k = m + n;

            ComplexMatrix<N> a = ComplexMatrix<N>.Zero(k, k);
            ComplexVector<N> c = cs[1..];

            for (int i = 0; i < m; i++) {
                a[i, i] = 1;
            }

            for (int i = m; i < k; i++) {
                for (int j = i - m, r = 0; j < k; j++, r++) {
                    a[j, i] = -cs[r];
                }
            }

            ComplexVector<N> v = ComplexMatrix<N>.Solve(a, c);
            ComplexVector<N> ms = ComplexVector<N>.Concat(cs[0], v[..m]), ns = ComplexVector<N>.Concat(1, v[m..]);

            return (ms, ns);
        }
    }
}
