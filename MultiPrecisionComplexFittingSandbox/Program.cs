using MultiPrecision;
using MultiPrecisionComplex;

namespace MultiPrecisionComplexFittingSandbox {
    internal class Program {
        static void Main() {
            for (int x = 256; x >= 8; x--) {
                Complex<Pow2.N16> y = ComplexGammaN16.Gamma((x + 0.5, 1));

                Console.WriteLine($"{(x + 0.5, 1)},{y}");

                if (Complex<Pow2.N16>.IsNaN(y)) {
                    break;
                }
            }

            Console.WriteLine("END");
            Console.Read();
        }
    }
}
