﻿using MultiPrecision;
using MultiPrecisionAlgebra;
using MultiPrecisionComplex;
using MultiPrecisionComplexAlgebra;
using MultiPrecisionComplexFitting;

namespace MultiPrecisionComplexFittingTests {
    [TestClass()]
    public class ComplexPadeFitterTests {
        [TestMethod()]
        public void FitWithInterceptTest() {
            Complex<Pow2.N8>[] xs = ComplexMatrix<Pow2.N8>.Flatten(ComplexMatrix<Pow2.N8>.Grid((-7, 9), (-6, 10)) / 16);
            Complex<Pow2.N8>[] ys = ComplexVector<Pow2.N8>.Func(x => Complex<Pow2.N8>.Cos(x * (0.5, 0.25)) - 0.25 + Complex<Pow2.N8>.ImaginaryOne, xs);

            ComplexPadeFitter<Pow2.N8> fitter = new(xs, ys, intercept: 0.75 + Complex<Pow2.N8>.ImaginaryOne, numer: 8, denom: 6);

            ComplexVector<Pow2.N8> parameters = fitter.Fit();

            Console.WriteLine($"Numer : {parameters[..fitter.Numer]}");
            Console.WriteLine($"Denom : {parameters[fitter.Numer..]}");

            Assert.AreEqual(0.75 + Complex<Pow2.N8>.ImaginaryOne, fitter.Regress(0, parameters));

            for (int i = 0; i < xs.Length; i++) {
                Assert.IsTrue(MultiPrecision<Pow2.N8>.Abs(ys[i].R - fitter.Regress(xs[i], parameters).R) < 1e-5,
                    $"\nexpected : {ys[i].R}\n actual  : {fitter.Regress(xs[i], parameters).R}"
                );
                Assert.IsTrue(MultiPrecision<Pow2.N8>.Abs(ys[i].I - fitter.Regress(xs[i], parameters).I) < 1e-5,
                    $"\nexpected : {ys[i].I}\n actual  : {fitter.Regress(xs[i], parameters).I}"
                );
            }

            Assert.IsTrue(fitter.Error(parameters).Norm < 1e-4);
        }

        [TestMethod()]
        public void FitWithoutInterceptTest() {
            Complex<Pow2.N8>[] xs = ComplexMatrix<Pow2.N8>.Flatten(ComplexMatrix<Pow2.N8>.Grid((-7, 9), (-6, 10)) / 16);
            Complex<Pow2.N8>[] ys = ComplexVector<Pow2.N8>.Func(x => Complex<Pow2.N8>.Cos(x * (0.5, 0.25)) - 0.25 + Complex<Pow2.N8>.ImaginaryOne, xs);

            ComplexPadeFitter<Pow2.N8> fitter = new(xs, ys, numer: 8, denom: 6);

            ComplexVector<Pow2.N8> parameters = fitter.Fit();

            Console.WriteLine($"Numer : {parameters[..fitter.Numer]}");
            Console.WriteLine($"Denom : {parameters[fitter.Numer..]}");

            for (int i = 0; i < xs.Length; i++) {
                Assert.IsTrue(MultiPrecision<Pow2.N8>.Abs(ys[i].R - fitter.Regress(xs[i], parameters).R) < 1e-5,
                    $"\nexpected : {ys[i].R}\n actual  : {fitter.Regress(xs[i], parameters).R}"
                );
                Assert.IsTrue(MultiPrecision<Pow2.N8>.Abs(ys[i].I - fitter.Regress(xs[i], parameters).I) < 1e-5,
                    $"\nexpected : {ys[i].I}\n actual  : {fitter.Regress(xs[i], parameters).I}"
                );
            }

            Assert.IsTrue(fitter.Error(parameters).Norm < 1e-4);
        }

        [TestMethod()]
        public void ExecuteWeightedFittingWithInterceptTest() {
            Complex<Pow2.N8>[] xs = ComplexMatrix<Pow2.N8>.Flatten(ComplexMatrix<Pow2.N8>.Grid((-7, 9), (-6, 10)) / 16);
            Complex<Pow2.N8>[] ys = ComplexVector<Pow2.N8>.Func(x => Complex<Pow2.N8>.Cos(x * (0.5, 0.25)) - 0.25 + Complex<Pow2.N8>.ImaginaryOne, xs);
            MultiPrecision<Pow2.N8>[] ws = Vector<Pow2.N8>.Fill(xs.Length, value: 0.5);

            ys[256] = 1e+8;
            ws[256] = 0;

            ComplexPadeFitter<Pow2.N8> fitter = new(xs, ys, intercept: 0.75 + Complex<Pow2.N8>.ImaginaryOne, numer: 8, denom: 6);

            ComplexVector<Pow2.N8> parameters = fitter.Fit(ws);

            Console.WriteLine($"Numer : {parameters[..fitter.Numer]}");
            Console.WriteLine($"Denom : {parameters[fitter.Numer..]}");

            Assert.AreEqual(0.75 + Complex<Pow2.N8>.ImaginaryOne, fitter.Regress(0, parameters));

            for (int i = 0; i < xs.Length; i++) {
                if (i == 256) {
                    continue;
                }

                Assert.IsTrue(MultiPrecision<Pow2.N8>.Abs(ys[i].R - fitter.Regress(xs[i], parameters).R) < 1e-5,
                    $"\nexpected : {ys[i].R}\n actual  : {fitter.Regress(xs[i], parameters).R}"
                );
                Assert.IsTrue(MultiPrecision<Pow2.N8>.Abs(ys[i].I - fitter.Regress(xs[i], parameters).I) < 1e-5,
                    $"\nexpected : {ys[i].I}\n actual  : {fitter.Regress(xs[i], parameters).I}"
                );
            }
        }

        [TestMethod()]
        public void ExecuteWeightedFittingWithoutInterceptTest() {
            Complex<Pow2.N8>[] xs = ComplexMatrix<Pow2.N8>.Flatten(ComplexMatrix<Pow2.N8>.Grid((-7, 9), (-6, 10)) / 16);
            Complex<Pow2.N8>[] ys = ComplexVector<Pow2.N8>.Func(x => Complex<Pow2.N8>.Cos(x * (0.5, 0.25)) - 0.25 + Complex<Pow2.N8>.ImaginaryOne, xs);
            MultiPrecision<Pow2.N8>[] ws = Vector<Pow2.N8>.Fill(xs.Length, value: 0.5);

            ys[256] = 1e+8;
            ws[256] = 0;

            ComplexPadeFitter<Pow2.N8> fitter = new(xs, ys, numer: 8, denom: 6);

            ComplexVector<Pow2.N8> parameters = fitter.Fit(ws);

            Console.WriteLine($"Numer : {parameters[..fitter.Numer]}");
            Console.WriteLine($"Denom : {parameters[fitter.Numer..]}");

            for (int i = 0; i < xs.Length; i++) {
                if (i == 256) {
                    continue;
                }

                Assert.IsTrue(MultiPrecision<Pow2.N8>.Abs(ys[i].R - fitter.Regress(xs[i], parameters).R) < 1e-5,
                    $"\nexpected : {ys[i].R}\n actual  : {fitter.Regress(xs[i], parameters).R}"
                );
                Assert.IsTrue(MultiPrecision<Pow2.N8>.Abs(ys[i].I - fitter.Regress(xs[i], parameters).I) < 1e-5,
                    $"\nexpected : {ys[i].I}\n actual  : {fitter.Regress(xs[i], parameters).I}"
                );
            }
        }

        [TestMethod()]
        public void FitWithInterceptRealTest() {
            Complex<Pow2.N8>[] xs = (ComplexVector<Pow2.N8>)Vector<Pow2.N8>.Arange(0, 256) / 256;
            Complex<Pow2.N8>[] ys = ComplexVector<Pow2.N8>.Func(x => Complex<Pow2.N8>.Cos(x) - 0.25, xs);

            ComplexPadeFitter<Pow2.N8> fitter = new(xs, ys, intercept: 0.75, numer: 8, denom: 6);

            ComplexVector<Pow2.N8> parameters = fitter.Fit();

            Console.WriteLine($"Numer : {parameters[..fitter.Numer]}");
            Console.WriteLine($"Denom : {parameters[fitter.Numer..]}");

            Assert.AreEqual(0.75, fitter.Regress(0, parameters));

            for (int i = 0; i < xs.Length; i++) {
                Assert.IsTrue(MultiPrecision<Pow2.N8>.Abs(ys[i].R - fitter.Regress(xs[i], parameters).R) < 1e-5,
                    $"\nexpected : {ys[i].R}\n actual  : {fitter.Regress(xs[i], parameters).R}"
                );
                Assert.IsTrue(MultiPrecision<Pow2.N8>.Abs(ys[i].I - fitter.Regress(xs[i], parameters).I) < 1e-5,
                    $"\nexpected : {ys[i].I}\n actual  : {fitter.Regress(xs[i], parameters).I}"
                );
            }

            Assert.IsTrue(fitter.Error(parameters).Norm < 1e-4);
        }

        [TestMethod()]
        public void FitWithoutInterceptRealTest() {
            Complex<Pow2.N8>[] xs = (ComplexVector<Pow2.N8>)Vector<Pow2.N8>.Arange(0, 256) / 256;
            Complex<Pow2.N8>[] ys = ComplexVector<Pow2.N8>.Func(x => Complex<Pow2.N8>.Cos(x) - 0.25, xs);

            ComplexPadeFitter<Pow2.N8> fitter = new(xs, ys, numer: 8, denom: 6);

            ComplexVector<Pow2.N8> parameters = fitter.Fit();

            Console.WriteLine($"Numer : {parameters[..fitter.Numer]}");
            Console.WriteLine($"Denom : {parameters[fitter.Numer..]}");

            for (int i = 0; i < xs.Length; i++) {
                Assert.IsTrue(MultiPrecision<Pow2.N8>.Abs(ys[i].R - fitter.Regress(xs[i], parameters).R) < 1e-5,
                    $"\nexpected : {ys[i].R}\n actual  : {fitter.Regress(xs[i], parameters).R}"
                );
                Assert.IsTrue(MultiPrecision<Pow2.N8>.Abs(ys[i].I - fitter.Regress(xs[i], parameters).I) < 1e-5,
                    $"\nexpected : {ys[i].I}\n actual  : {fitter.Regress(xs[i], parameters).I}"
                );
            }

            Assert.IsTrue(fitter.Error(parameters).Norm < 1e-4);
        }

        [TestMethod()]
        public void ExecuteWeightedFittingWithInterceptRealTest() {
            Complex<Pow2.N8>[] xs = (ComplexVector<Pow2.N8>)Vector<Pow2.N8>.Arange(0, 256) / 256;
            Complex<Pow2.N8>[] ys = ComplexVector<Pow2.N8>.Func(x => Complex<Pow2.N8>.Cos(x) - 0.25, xs);
            MultiPrecision<Pow2.N8>[] ws = Vector<Pow2.N8>.Fill(xs.Length, value: 0.5);

            ys[64] = 1e+8;
            ws[64] = 0;

            ComplexPadeFitter<Pow2.N8> fitter = new(xs, ys, intercept: 0.75, numer: 8, denom: 6);

            ComplexVector<Pow2.N8> parameters = fitter.Fit(ws);

            Console.WriteLine($"Numer : {parameters[..fitter.Numer]}");
            Console.WriteLine($"Denom : {parameters[fitter.Numer..]}");

            Assert.AreEqual(0.75, fitter.Regress(0, parameters));

            for (int i = 0; i < xs.Length; i++) {
                if (i == 64) {
                    continue;
                }

                Assert.IsTrue(MultiPrecision<Pow2.N8>.Abs(ys[i].R - fitter.Regress(xs[i], parameters).R) < 1e-5,
                    $"\nexpected : {ys[i].R}\n actual  : {fitter.Regress(xs[i], parameters).R}"
                );
                Assert.IsTrue(MultiPrecision<Pow2.N8>.Abs(ys[i].I - fitter.Regress(xs[i], parameters).I) < 1e-5,
                    $"\nexpected : {ys[i].I}\n actual  : {fitter.Regress(xs[i], parameters).I}"
                );
            }
        }

        [TestMethod()]
        public void ExecuteWeightedFittingWithoutInterceptRealTest() {
            Complex<Pow2.N8>[] xs = (ComplexVector<Pow2.N8>)Vector<Pow2.N8>.Arange(0, 256) / 256;
            Complex<Pow2.N8>[] ys = ComplexVector<Pow2.N8>.Func(x => Complex<Pow2.N8>.Cos(x) - 0.25, xs);
            MultiPrecision<Pow2.N8>[] ws = Vector<Pow2.N8>.Fill(xs.Length, value: 0.5);

            ys[64] = 1e+8;
            ws[64] = 0;

            ComplexPadeFitter<Pow2.N8> fitter = new(xs, ys, numer: 8, denom: 6);

            ComplexVector<Pow2.N8> parameters = fitter.Fit(ws);

            Console.WriteLine($"Numer : {parameters[..fitter.Numer]}");
            Console.WriteLine($"Denom : {parameters[fitter.Numer..]}");

            for (int i = 0; i < xs.Length; i++) {
                if (i == 64) {
                    continue;
                }

                Assert.IsTrue(MultiPrecision<Pow2.N8>.Abs(ys[i].R - fitter.Regress(xs[i], parameters).R) < 1e-5,
                    $"\nexpected : {ys[i].R}\n actual  : {fitter.Regress(xs[i], parameters).R}"
                );
                Assert.IsTrue(MultiPrecision<Pow2.N8>.Abs(ys[i].I - fitter.Regress(xs[i], parameters).I) < 1e-5,
                    $"\nexpected : {ys[i].I}\n actual  : {fitter.Regress(xs[i], parameters).I}"
                );
            }
        }

        [TestMethod()]
        public void FitWithInterceptImagTest() {
            Complex<Pow2.N8>[] xs = (ComplexVector<Pow2.N8>)Vector<Pow2.N8>.Arange(0, 256) / 256;
            Complex<Pow2.N8>[] ys = ComplexVector<Pow2.N8>.Func(x => (Complex<Pow2.N8>.Cos(x) - 0.25) * (0, 1), xs);

            ComplexPadeFitter<Pow2.N8> fitter = new(xs, ys, intercept: (0, 0.75), numer: 8, denom: 6);

            ComplexVector<Pow2.N8> parameters = fitter.Fit();

            Console.WriteLine($"Numer : {parameters[..fitter.Numer]}");
            Console.WriteLine($"Denom : {parameters[fitter.Numer..]}");

            Assert.AreEqual((0, 0.75), fitter.Regress(0, parameters));

            for (int i = 0; i < xs.Length; i++) {
                Assert.IsTrue(MultiPrecision<Pow2.N8>.Abs(ys[i].R - fitter.Regress(xs[i], parameters).R) < 1e-5,
                    $"\nexpected : {ys[i].R}\n actual  : {fitter.Regress(xs[i], parameters).R}"
                );
                Assert.IsTrue(MultiPrecision<Pow2.N8>.Abs(ys[i].I - fitter.Regress(xs[i], parameters).I) < 1e-5,
                    $"\nexpected : {ys[i].I}\n actual  : {fitter.Regress(xs[i], parameters).I}"
                );
            }

            Assert.IsTrue(fitter.Error(parameters).Norm < 1e-4);
        }

        [TestMethod()]
        public void FitWithoutInterceptImagTest() {
            Complex<Pow2.N8>[] xs = (ComplexVector<Pow2.N8>)Vector<Pow2.N8>.Arange(0, 256) / 256;
            Complex<Pow2.N8>[] ys = ComplexVector<Pow2.N8>.Func(x => (Complex<Pow2.N8>.Cos(x) - 0.25) * (0, 1), xs);

            ComplexPadeFitter<Pow2.N8> fitter = new(xs, ys, numer: 8, denom: 6);

            ComplexVector<Pow2.N8> parameters = fitter.Fit();

            Console.WriteLine($"Numer : {parameters[..fitter.Numer]}");
            Console.WriteLine($"Denom : {parameters[fitter.Numer..]}");

            for (int i = 0; i < xs.Length; i++) {
                Assert.IsTrue(MultiPrecision<Pow2.N8>.Abs(ys[i].R - fitter.Regress(xs[i], parameters).R) < 1e-5,
                    $"\nexpected : {ys[i].R}\n actual  : {fitter.Regress(xs[i], parameters).R}"
                );
                Assert.IsTrue(MultiPrecision<Pow2.N8>.Abs(ys[i].I - fitter.Regress(xs[i], parameters).I) < 1e-5,
                    $"\nexpected : {ys[i].I}\n actual  : {fitter.Regress(xs[i], parameters).I}"
                );
            }

            Assert.IsTrue(fitter.Error(parameters).Norm < 1e-4);
        }

        [TestMethod()]
        public void ExecuteWeightedFittingWithInterceptImagTest() {
            Complex<Pow2.N8>[] xs = (ComplexVector<Pow2.N8>)Vector<Pow2.N8>.Arange(0, 256) / 256;
            Complex<Pow2.N8>[] ys = ComplexVector<Pow2.N8>.Func(x => (Complex<Pow2.N8>.Cos(x) - 0.25) * (0, 1), xs);
            MultiPrecision<Pow2.N8>[] ws = Vector<Pow2.N8>.Fill(xs.Length, value: 0.5);

            ys[64] = 1e+8;
            ws[64] = 0;

            ComplexPadeFitter<Pow2.N8> fitter = new(xs, ys, intercept: (0, 0.75), numer: 8, denom: 6);

            ComplexVector<Pow2.N8> parameters = fitter.Fit(ws);

            Console.WriteLine($"Numer : {parameters[..fitter.Numer]}");
            Console.WriteLine($"Denom : {parameters[fitter.Numer..]}");

            Assert.AreEqual((0, 0.75), fitter.Regress(0, parameters));

            for (int i = 0; i < xs.Length; i++) {
                if (i == 64) {
                    continue;
                }

                Assert.IsTrue(MultiPrecision<Pow2.N8>.Abs(ys[i].R - fitter.Regress(xs[i], parameters).R) < 1e-5,
                    $"\nexpected : {ys[i].R}\n actual  : {fitter.Regress(xs[i], parameters).R}"
                );
                Assert.IsTrue(MultiPrecision<Pow2.N8>.Abs(ys[i].I - fitter.Regress(xs[i], parameters).I) < 1e-5,
                    $"\nexpected : {ys[i].I}\n actual  : {fitter.Regress(xs[i], parameters).I}"
                );
            }
        }

        [TestMethod()]
        public void ExecuteWeightedFittingWithoutInterceptImagTest() {
            Complex<Pow2.N8>[] xs = (ComplexVector<Pow2.N8>)Vector<Pow2.N8>.Arange(0, 256) / 256;
            Complex<Pow2.N8>[] ys = ComplexVector<Pow2.N8>.Func(x => (Complex<Pow2.N8>.Cos(x) - 0.25) * (0, 1), xs);
            MultiPrecision<Pow2.N8>[] ws = Vector<Pow2.N8>.Fill(xs.Length, value: 0.5);

            ys[64] = 1e+8;
            ws[64] = 0;

            ComplexPadeFitter<Pow2.N8> fitter = new(xs, ys, numer: 8, denom: 6);

            ComplexVector<Pow2.N8> parameters = fitter.Fit(ws);

            Console.WriteLine($"Numer : {parameters[..fitter.Numer]}");
            Console.WriteLine($"Denom : {parameters[fitter.Numer..]}");

            for (int i = 0; i < xs.Length; i++) {
                if (i == 64) {
                    continue;
                }

                Assert.IsTrue(MultiPrecision<Pow2.N8>.Abs(ys[i].R - fitter.Regress(xs[i], parameters).R) < 1e-5,
                    $"\nexpected : {ys[i].R}\n actual  : {fitter.Regress(xs[i], parameters).R}"
                );
                Assert.IsTrue(MultiPrecision<Pow2.N8>.Abs(ys[i].I - fitter.Regress(xs[i], parameters).I) < 1e-5,
                    $"\nexpected : {ys[i].I}\n actual  : {fitter.Regress(xs[i], parameters).I}"
                );
            }
        }
    }
}