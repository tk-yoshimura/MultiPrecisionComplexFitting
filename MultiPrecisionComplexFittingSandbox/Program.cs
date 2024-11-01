using MultiPrecision;
using MultiPrecisionComplex;
using MultiPrecisionComplexAlgebra;
using MultiPrecisionComplexFitting;
using MultiPrecisionCurveFitting;

namespace MultiPrecisionComplexFittingSandbox {
    internal class Program {
        static void Main() {
            ComplexVector<Pow2.N16> xs = ComplexMatrix<Pow2.N16>.Flatten(ComplexMatrix<Pow2.N16>.Grid((-32, 32), (-32, 32)) / 256);
            ComplexVector<Pow2.N16> ys = (x => x != 0 ? Complex<Pow2.N16>.Log(1 + x) / x : 1, xs);

            ComplexSumTable<Pow2.N64> table = new(xs.Convert<Pow2.N64>(), ys.Convert<Pow2.N64>());

            bool finished = false; 
            for (int coefs = 5; coefs <= 256 && !finished; coefs++) {
                foreach ((int m, int n) in CurveFittingUtils.EnumeratePadeDegree(coefs, 2)) {
                    ComplexPadeFitter<Pow2.N64> pade = new(table, m, n, intercept: 1);

                    ComplexVector<Pow2.N64> param = pade.Fit();
                    ComplexVector<Pow2.N64> errs = pade.Error(param);

                    ComplexVector<Pow2.N16> vs = pade.Regress(xs.Convert<Pow2.N64>(), param).Convert<Pow2.N16>();

                    MultiPrecision<Pow2.N16> max_rateerr = ComplexCurveFittingUtils.MaxRelativeError(ys, vs);

                    Console.WriteLine($"m={m},n={n}");
                    Console.WriteLine($"{max_rateerr:e20}");

                    if (max_rateerr > "1e-25") {
                        coefs += 4;
                        break;
                    }

                    if (max_rateerr < "1e-40") {
                        break;
                    }

                    if (max_rateerr < "1e-30") {
                        finished = true;

                        using (StreamWriter sw = new("../../../../results_disused/complex_log_precision16.csv")) {
                            sw.WriteLine($"x=[{xs[0]},{xs[^1]}]");
                            sw.WriteLine($"samples={xs.Dim}");
                            sw.WriteLine($"m={m},n={n}");
                            sw.WriteLine("numer");
                            foreach (var (_, val) in param[..m]) {
                                sw.WriteLine($"{val:e38}");
                            }
                            sw.WriteLine("denom");
                            foreach (var (_, val) in param[m..]) {
                                sw.WriteLine($"{val:e38}");
                            }

                            sw.WriteLine("relative err");
                            sw.WriteLine($"{max_rateerr:e20}");
                            sw.Flush();
                        }

                        break;
                    }
                }
            }

            Console.WriteLine("END");
            Console.Read();
        }
    }
}
