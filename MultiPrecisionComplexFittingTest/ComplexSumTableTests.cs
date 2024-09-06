using MultiPrecision;
using MultiPrecisionComplex;
using MultiPrecisionComplexAlgebra;
using MultiPrecisionComplexFitting;

namespace MultiPrecisionComplexFittingTests {
    [TestClass()]
    public class ComplexSumTableTests {
        [TestMethod()]
        public void IndexerTest() {
            (Complex<Pow2.N4> x, Complex<Pow2.N4> y)[] vs = [
                ((2, 1), (11, -1)), ((3, -2), (13, 3)), ((5, 4), (17, -3)), ((7, -5), (19, -6))
            ];
            Complex<Pow2.N4> s(int xn, int xn_conj, int yn, int yn_conj) {
                return (new ComplexVector<Pow2.N4>(vs.Select(
                    v =>
                    Complex<Pow2.N4>.Pow(v.x, xn) *
                    Complex<Pow2.N4>.Pow(Complex<Pow2.N4>.Conjugate(v.x), xn_conj) *
                    Complex<Pow2.N4>.Pow(v.y, yn) *
                    Complex<Pow2.N4>.Pow(Complex<Pow2.N4>.Conjugate(v.y), yn_conj)
                ))).Sum;
            };

            ComplexSumTable<Pow2.N4> table = new(vs.Select(v => v.x).ToArray(), vs.Select(v => v.y).ToArray());

            for (int i = 0; i <= 12; i++) {
                for (int j = 0; j <= i; j++) {
                    for (int k = 0; k <= j; k++) {
                        for (int m = 0; m <= k; m++) {
                            Assert.IsTrue((s(i, j, k, m) - table[i, j, k, m]).Norm < 1e-25 * s(i, j, k, m).Norm);
                            Assert.IsTrue((s(j, i, m, k) - table[j, i, m, k]).Norm < 1e-25 * s(j, i, m, k).Norm);
                        }
                    }
                }
            }
        }

        [TestMethod()]
        public void ReverseIndexerTest() {
            (Complex<Pow2.N4> x, Complex<Pow2.N4> y)[] vs = [
                ((2, 1), (11, -1)), ((3, -2), (13, 3)), ((5, 4), (17, -3)), ((7, -5), (19, -6))
            ];
            Complex<Pow2.N4> s(int xn, int xn_conj, int yn, int yn_conj) {
                return (new ComplexVector<Pow2.N4>(vs.Select(
                    v =>
                    Complex<Pow2.N4>.Pow(v.x, xn) *
                    Complex<Pow2.N4>.Pow(Complex<Pow2.N4>.Conjugate(v.x), xn_conj) *
                    Complex<Pow2.N4>.Pow(v.y, yn) *
                    Complex<Pow2.N4>.Pow(Complex<Pow2.N4>.Conjugate(v.y), yn_conj)
                ))).Sum;
            };

            ComplexSumTable<Pow2.N4> table = new(vs.Select(v => v.x).ToArray(), vs.Select(v => v.y).ToArray());

            for (int i = 12; i >= 0; i--) {
                for (int j = i; j >= 0; j--) {
                    for (int k = j; k >= 0; k--) {
                        for (int m = k; m >= 0; m--) {
                            Assert.IsTrue((s(i, j, k, m) - table[i, j, k, m]).Norm < 1e-25 * s(i, j, k, m).Norm);
                            Assert.IsTrue((s(j, i, m, k) - table[j, i, m, k]).Norm < 1e-25 * s(j, i, m, k).Norm);
                        }
                    }
                }
            }
        }
    }
}