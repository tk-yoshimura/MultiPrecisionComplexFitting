# MultiPrecisionComplexFitting
 MultiPrecision Complex Fitting Implements 

## Requirement
.NET 8.0  
AVX2 suppoted CPU. (Intel:Haswell(2013)-, AMD:Excavator(2015)-)  
[MultiPrecision](https://github.com/tk-yoshimura/MultiPrecision)  
[MultiPrecisionComplex](https://github.com/tk-yoshimura/MultiPrecisionComplex)  
[MultiPrecisionAlgebra](https://github.com/tk-yoshimura/MultiPrecisionAlgebra)  
[MultiPrecisionComplexAlgebra](https://github.com/tk-yoshimura/MultiPrecisionComplexAlgebra)  

## Install

[Download DLL](https://github.com/tk-yoshimura/MultiPrecisionComplexFitting/releases)  
[Download Nuget](https://www.nuget.org/packages/tyoshimura.multiprecision.complexfitting/)

## Usage

```csharp
Complex<Pow2.N8>[] xs = ComplexMatrix<Pow2.N8>.Flatten(ComplexMatrix<Pow2.N8>.Grid((-7, 9), (-6, 10)) / 16);
Complex<Pow2.N8>[] ys = ComplexVector<Pow2.N8>.Func(x => Complex<Pow2.N8>.Cos(x * (0.5, 0.25)) - 0.25 + Complex<Pow2.N8>.ImaginaryOne, xs);

ComplexPadeFitter<Pow2.N8> fitter = new(xs, ys, numer: 8, denom: 6);

ComplexVector<Pow2.N8> parameters = fitter.Fit();

Console.WriteLine($"Numer : {parameters[..fitter.Numer]}");
Console.WriteLine($"Denom : {parameters[fitter.Numer..]}");

ComplexVector<Pow2.N8> vs = fitter.Regress(xs, parameters);
```

## Licence
[MIT](https://github.com/tk-yoshimura/MultiPrecisionComplexFitting/blob/master/LICENSE)

## Author

[T.Yoshimura](https://github.com/tk-yoshimura)
