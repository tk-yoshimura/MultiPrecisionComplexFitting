using MultiPrecision;
using MultiPrecisionAlgebra;
using MultiPrecisionComplex;
using MultiPrecisionComplexAlgebra;
using System;
using System.Linq;

namespace MultiPrecisionComplexFitting {

    public class ComplexPolynomialFitter<N> : ComplexFitter<N> where N : struct, IConstant {

        private readonly ComplexSumTable<N> sum_table;
        private readonly Complex<N>? intercept;

        public int Degree {
            get; private set;
        }

        public ComplexPolynomialFitter(ComplexVector<N> xs, ComplexVector<N> ys, int degree, Complex<N>? intercept = null)
            : base(xs, ys, parameters: checked(degree + 1)) {

            this.sum_table = new(X, (intercept is null) ? ys : ys.Select(y => y.val - intercept).ToArray());
            this.intercept = intercept;
            this.Degree = degree;
        }

        public override Complex<N> Regress(Complex<N> x, ComplexVector<N> parameters) {
            if (parameters.Dim != Parameters) {
                throw new ArgumentException("invalid size", nameof(parameters));
            }

            Complex<N> y = ComplexVector<N>.Polynomial(x, parameters);

            return y;
        }

        public ComplexVector<N> Fit(Vector<N>? weights = null) {
            sum_table.W = weights;
            (ComplexMatrix<N> m, ComplexVector<N> v) = GenerateTable(sum_table, Degree, enable_intercept: intercept is null);

            if (intercept is null) {

                ComplexVector<N> parameters = ComplexMatrix<N>.SolvePositiveSymmetric(m, v,
#if DEBUG
                    enable_check_hermitian: true
#else
                    enable_check_hermitian: false
#endif
                );

                return parameters;
            }
            else {

                ComplexVector<N> parameters = ComplexVector<N>.Concat(
                    intercept,
                    ComplexMatrix<N>.SolvePositiveSymmetric(m, v,
#if DEBUG
                    enable_check_hermitian: true
#else
                    enable_check_hermitian: false
#endif
                    )
                );

                return parameters;
            }
        }

        internal static (ComplexMatrix<N> m, ComplexVector<N> v) GenerateTable(ComplexSumTable<N> sum_table, int degree, bool enable_intercept) {
            int dim = degree + (enable_intercept ? 1 : 0);

            Complex<N>[,] m = new Complex<N>[dim, dim];
            Complex<N>[] v = new Complex<N>[dim];

            if (enable_intercept) {
                for (int i = 0; i < dim; i++) {
                    for (int j = i; j < dim; j++) {
                        m[i, j] = sum_table[j, i, 0, 0];

                        if (i != j) {
                            m[j, i] = Complex<N>.Conjugate(m[i, j]);
                        }
                        else {
                            m[j, i] = m[j, i].R;
                        }
                    }
                }

                for (int i = 0; i < dim; i++) {
                    v[i] = sum_table[0, i, 1, 0];
                }
            }
            else {
                for (int i = 0; i < dim; i++) {
                    for (int j = i; j < dim; j++) {
                        m[i, j] = sum_table[j + 1, i + 1, 0, 0];

                        if (i != j) {
                            m[j, i] = Complex<N>.Conjugate(m[i, j]);
                        }
                        else {
                            m[j, i] = m[j, i].R;
                        }
                    }
                }

                for (int i = 0; i < dim; i++) {
                    v[i] = sum_table[0, i + 1, 1, 0];
                }
            }

            return (m, v);
        }
    }
}
