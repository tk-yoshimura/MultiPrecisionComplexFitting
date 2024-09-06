using MultiPrecision;
using MultiPrecisionComplex;
using MultiPrecisionComplexAlgebra;
using System;

namespace MultiPrecisionComplexFitting {
    public abstract class ComplexFitter<N> where N : struct, IConstant {
        public ComplexVector<N> X { get; private set; }

        public ComplexVector<N> Y { get; private set; }

        public int Points { get; private set; }

        public int Parameters { get; private set; }

        public ComplexFitter(ComplexVector<N> xs, ComplexVector<N> ys, int parameters) {
            if (xs.Dim < parameters || xs.Dim != ys.Dim) {
                throw new ArgumentException("mismatch size", $"{nameof(xs)},{nameof(ys)}");
            }
            ArgumentOutOfRangeException.ThrowIfLessThan(parameters, 1, nameof(parameters));

            this.X = xs.Copy();
            this.Y = ys.Copy();
            this.Points = xs.Dim;
            this.Parameters = parameters;
        }

        public MultiPrecision<N> Cost(ComplexVector<N> parameters) {
            if (parameters.Dim != Parameters) {
                throw new ArgumentException("invalid size", nameof(parameters));
            }

            ComplexVector<N> errors = Error(parameters);
            MultiPrecision<N> cost = errors.SquareNorm;

            return cost;
        }

        public ComplexVector<N> Error(ComplexVector<N> parameters) {
            if (parameters.Dim != Parameters) {
                throw new ArgumentException("invalid size", nameof(parameters));
            }

            ComplexVector<N> errors = Regress(X, parameters) - Y;

            return errors;
        }

        public abstract Complex<N> Regress(Complex<N> x, ComplexVector<N> parameters);

        public ComplexVector<N> Regress(ComplexVector<N> xs, ComplexVector<N> parameters) {
            Complex<N>[] ys = new Complex<N>[xs.Dim];

            for (int i = 0; i < xs.Dim; i++) {
                ys[i] = Regress(xs[i], parameters);
            }

            return ys;
        }
    }
}
