using MultiPrecision;
using MultiPrecisionAlgebra;
using MultiPrecisionComplex;
using MultiPrecisionComplexAlgebra;

namespace MultiPrecisionComplexFitting {
    public class ComplexSumTable<N> where N : struct, IConstant {
        private readonly List<ComplexVector<N>> xs = new(), ys = new(), xs_conj = new(), ys_conj = new();
        private Dictionary<(int xn, int xn_conj, int yn, int yn_conj), Complex<N>> table;

        private Vector<N>? w = null;

        internal ComplexVector<N> X { get; }
        internal ComplexVector<N> Y { get; }

        public ComplexSumTable(ComplexVector<N> x, ComplexVector<N> y) {
            if (x.Dim != y.Dim) {
                throw new ArgumentException("invalid size", $"{nameof(x)},{nameof(y)}");
            }

            this.xs.Add(x);
            this.ys.Add(y);
            this.xs_conj.Add(x.Conj);
            this.ys_conj.Add(y.Conj);
            this.table = new() {
                { (0, 0, 0, 0), x.Dim },
            };

            this.X = x;
            this.Y = y;
        }

        public Complex<N> this[int xn, int xn_conj, int yn, int yn_conj] {
            get {
                if (xn < 0 || yn < 0) {
                    throw new ArgumentOutOfRangeException($"{nameof(xn)},{nameof(yn)}");
                }
                if (xn_conj < 0 || yn_conj < 0) {
                    throw new ArgumentOutOfRangeException($"{nameof(xn_conj)},{nameof(yn_conj)}");
                }

                for (int i = xs.Count; i < xn; i++) {
                    int xn0 = (i + 1) / 2 - 1, xn1 = i - xn0 - 1;

                    xs.Add(xs[xn0] * xs[xn1]);
                }

                for (int i = xs_conj.Count; i < xn_conj; i++) {
                    int xn0 = (i + 1) / 2 - 1, xn1 = i - xn0 - 1;

                    xs_conj.Add(xs_conj[xn0] * xs_conj[xn1]);
                }

                for (int i = ys.Count; i < yn; i++) {
                    int yn0 = (i + 1) / 2 - 1, yn1 = i - yn0 - 1;

                    ys.Add(ys[yn0] * ys[yn1]);
                }

                for (int i = ys_conj.Count; i < yn_conj; i++) {
                    int yn0 = (i + 1) / 2 - 1, yn1 = i - yn0 - 1;

                    ys_conj.Add(ys_conj[yn0] * ys_conj[yn1]);
                }

                if (!table.ContainsKey((xn, xn_conj, yn, yn_conj))) {
                    Complex<N> s;

                    ComplexVector<N> v = ComplexVector<N>.Fill(X.Dim, 1);

                    if (xn > 0) {
                        v *= xs[xn - 1];
                    }
                    if (xn_conj > 0) {
                        v *= xs_conj[xn_conj - 1];
                    }
                    if (yn > 0) {
                        v *= ys[yn - 1];
                    }
                    if (yn_conj > 0) {
                        v *= ys_conj[yn_conj - 1];
                    }

                    s = w is null ? v.Sum : (v * w).Sum;

                    table.Add((xn, xn_conj, yn, yn_conj), s);
                }

                return table[(xn, xn_conj, yn, yn_conj)];
            }
        }

        public Vector<N>? W {
            get => w;
            set {
                if (value is not null && xs[0].Dim != value.Dim) {
                    throw new ArgumentException("invalid size", nameof(W));
                }

                this.w = value;
                this.table = new() {
                    { (0, 0, 0, 0), w is null ? xs[0].Dim : w.Sum },
                };
            }
        }
    }
}
