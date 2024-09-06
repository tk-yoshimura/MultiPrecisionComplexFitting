using MultiPrecision;
using MultiPrecisionAlgebra;
using MultiPrecisionComplex;
using MultiPrecisionComplexAlgebra;

namespace MultiPrecisionComplexFitting {
    public class ComplexPadeFitter<N> : ComplexFitter<N> where N : struct, IConstant {

        private readonly ComplexSumTable<N> sum_table;
        private readonly Complex<N>? intercept;

        public int Numer { get; private set; }

        public int Denom { get; private set; }

        public ComplexPadeFitter(ComplexVector<N> xs, ComplexVector<N> ys, int numer, int denom, Complex<N>? intercept = null)
            : base(xs, ys,
                  parameters:
                  (numer >= 2 && denom >= 2)
                      ? (numer + denom)
                      : throw new ArgumentOutOfRangeException($"{nameof(numer)},{nameof(denom)}")) {

            this.sum_table = new(X, Y);
            this.intercept = intercept;
            this.Numer = numer;
            this.Denom = denom;
        }

        public ComplexPadeFitter(ComplexSumTable<N> sum_table, int numer, int denom, Complex<N>? intercept = null)
            : base(sum_table.X, sum_table.Y,
                  parameters:
                  (numer >= 2 && denom >= 2)
                      ? (numer + denom)
                      : throw new ArgumentOutOfRangeException($"{nameof(numer)},{nameof(denom)}")) {

            this.sum_table = sum_table;
            this.intercept = intercept;
            this.Numer = numer;
            this.Denom = denom;
        }

        public override Complex<N> Regress(Complex<N> x, ComplexVector<N> parameters) {
            if (parameters.Dim != Parameters) {
                throw new ArgumentException("invalid size", nameof(parameters));
            }

            (Complex<N> numer, Complex<N> denom) = Fraction(x, parameters);

            Complex<N> y = numer / denom;

            return y;
        }

        public (Complex<N> numer, Complex<N> denom) Fraction(Complex<N> x, ComplexVector<N> parameters) {
            Complex<N> n = ComplexVector<N>.Polynomial(x, parameters[..Numer]);
            Complex<N> d = ComplexVector<N>.Polynomial(x, parameters[Numer..]);

            return (n, d);
        }

        public ComplexVector<N> Fit(Vector<N>? weights = null, MultiPrecision<N>? norm_cost = null) {
            sum_table.W = weights;
            (ComplexMatrix<N> m, ComplexVector<N> v) = GenerateTable(sum_table, Numer, Denom);

            if (norm_cost is not null) {
                Complex<N> c = norm_cost * sum_table[0, 0, 0, 0];

                for (int i = 0; i < m.Rows; i++) {
                    m[i, i] += c;
                }
            }

            if (intercept is null) {
                ComplexVector<N> x = ComplexMatrix<N>.SolvePositiveSymmetric(m, v);

                ComplexVector<N> parameters = ComplexVector<N>.Concat(x[..Numer], 1, x[Numer..]);

                return parameters;
            }
            else {
                v = v[1..] - intercept * m[0, 1..].Conj;
                m = m[1.., 1..];

                ComplexVector<N> x = ComplexMatrix<N>.SolvePositiveSymmetric(m, v);

                ComplexVector<N> parameters = ComplexVector<N>.Concat(intercept, x[..(Numer - 1)], 1, x[(Numer - 1)..]);

                return parameters;
            }
        }

        internal static (ComplexMatrix<N> m, ComplexVector<N> v) GenerateTable(ComplexSumTable<N> sum_table, int numer, int denom) {
            int dim = numer + denom - 1;

            Complex<N>[,] m = new Complex<N>[dim, dim];
            for (int i = 0, n = numer; i < n; i++) {
                for (int j = i; j < n; j++) {
                    m[i, j] = sum_table[j, i, 0, 0];

                    if (i != j) {
                        m[j, i] = Complex<N>.Conjugate(m[i, j]);
                    }
                    else {
                        m[j, i] = m[j, i].R;
                    }
                }
            }
            for (int i = numer, n = dim; i < n; i++) {
                for (int j = 0; j < numer; j++) {
                    m[i, j] = -sum_table[j + 1, i - numer, 0, 1];

                    if (i != j) {
                        m[j, i] = Complex<N>.Conjugate(m[i, j]);
                    }
                }
            }
            for (int i = numer, n = dim; i < n; i++) {
                for (int j = i; j < n; j++) {
                    m[i, j] = sum_table[j - numer + 1, i - numer + 1, 1, 1];

                    if (i != j) {
                        m[j, i] = Complex<N>.Conjugate(m[i, j]);
                    }
                    else {
                        m[j, i] = m[j, i].R;
                    }
                }
            }

            Complex<N>[] v = new Complex<N>[numer + denom - 1];
            for (int i = 0; i < numer; i++) {
                v[i] = sum_table[0, i, 1, 0];
            }
            for (int i = numer; i < dim; i++) {
                v[i] = -sum_table[1, i - numer, 1, 1];
            }

            return (m, v);
        }
    }
}
