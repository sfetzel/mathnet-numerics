// <copyright file="AmvwSolverTest.cs" company="Math.NET">
// Math.NET Numerics, part of the Math.NET Project
// http://numerics.mathdotnet.com
// http://github.com/mathnet/mathnet-numerics
//
// Copyright (c) 2009-2016 Math.NET
//
// Permission is hereby granted, free of charge, to any person
// obtaining a copy of this software and associated documentation
// files (the "Software"), to deal in the Software without
// restriction, including without limitation the rights to use,
// copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the
// Software is furnished to do so, subject to the following
// conditions:
//
// The above copyright notice and this permission notice shall be
// included in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
// EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
// OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
// NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
// HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
// WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
// OTHER DEALINGS IN THE SOFTWARE.
// </copyright>

using System;
using System.Collections.Generic;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Complex;
using MathNet.Numerics.LinearAlgebra.Complex.Solvers;
using MathNet.Numerics.Providers.LinearAlgebra;
using MathNet.Numerics;
using NUnit.Framework;
using System.Data;
using System.Collections.Immutable;

namespace MathNet.Numerics.Tests.LinearAlgebraTests.Complex
{
    using static MathNet.Numerics.LinearAlgebra.Complex.Solvers.AmvwSolver;
    using Complex = System.Numerics.Complex;

    /// <summary>
    /// Amvw solver tests.
    /// </summary>
    [TestFixture, Category("LA")]
    public class AmvwSolverTest
    {
        // f(x) = (1+2j) + (3+4j)*x + (5+6j)*x^2 + x^3
        // f(x) = (x - (-4.3533-5.9811i))*(x-(-0.37265+0.43338i))*(x-(-0.27401-0.45232i))
        // n = 3.
        public Complex[] GetTestCoefficients() => new Complex[]
        {
            new Complex(1, 2),
            new Complex(3, 4),
            new Complex(5, 6)
        };

        public Complex[] GetTestSolution() => new Complex[]
        {
            new Complex(-4.353335455010091, -5.981059347778331),
            new Complex(-0.372651320121572, 0.433379276936973),
            new Complex(-0.274013224868337, -0.452319929158646)
        };

        [Test]
        public void ApplyGivensRotation()
        {
            var rotation = new GivensRotation(new Complex(0.5919, -0.0143), 0.8059);
            var vector = (new Complex(3, 4), new Complex(7, 8));
            (var a, var b) = rotation * vector;

            Assert.That(a.Real, Is.EqualTo(-3.8083).Within(1e-3));
            Assert.That(a.Imaginary, Is.EqualTo(-4.1221).Within(1e-3));
            Assert.That(b.Real, Is.EqualTo(6.4470).Within(1e-3));
            Assert.That(b.Imaginary, Is.EqualTo(8.0587).Within(1e-3));
        }

        [Test]
        public void GetMatrix()
        {
            var cosine = new Complex(1, 2);
            var sine = 3;
            var rotation = new GivensRotation(cosine, sine);
            var matrix = rotation.GetMatrix(4, 1);
            Assert.Multiple(() =>
            {
                Assert.That(matrix[0, 0], Is.EqualTo(Complex.One));
                Assert.That(matrix[3, 3], Is.EqualTo(Complex.One));

                Assert.That(matrix[1, 1], Is.EqualTo(cosine));
                Assert.That(matrix[2, 2], Is.EqualTo(cosine.Conjugate()));
                Assert.That(matrix[1, 2], Is.EqualTo(new Complex(-sine, 0)));
                Assert.That(matrix[2, 1], Is.EqualTo(new Complex(sine, 0)));
            });
        }

        /// <summary>
        /// Checks if the method GivensRotation.Create calculates a rotation,
        /// which zeros the second component of the vector (1+2j, 3+4j).
        /// </summary>
        [TestCase(1, 2, 3, 4)]
        [TestCase(1, 2, 0, 0)]
        [TestCase(3, 4, 7, 8)]
        [TestCase(-2, 1, -3, 9)]
        [TestCase(-2, -1, -3, -9)]
        public void CalculateGivensRotationTest(double ar, double ai, double br, double bi)
        {
            (var a, var b) = (new Complex(ar, ai), new Complex(br, bi));
            (var rotation, var norm) = AmvwSolver.GivensRotation.Create(ref a, ref b);
            Assert.That(norm, Is.EqualTo(Math.Sqrt(ar * ar + ai * ai + br * br + bi * bi)).Within(1e-12));
            (var c, var d) = rotation * (new Complex(ar, ai), new Complex(br, bi));

            Assert.That(a.Magnitude, Is.EqualTo(norm).Within(1e-12));
            Assert.That(b.Magnitude, Is.EqualTo(0).Within(1e-12));

            (a, b) = (new Complex(ar, ai), new Complex(br, bi));
            Assert.Multiple(() =>
            {
                Assert.That(rotation.Cosine.MagnitudeSquared() + rotation.Sine * rotation.Sine, Is.EqualTo(1).Within(1e-12));

                Assert.That(c.Magnitude, Is.EqualTo(norm).Within(1e-12));
                Assert.That(c.Real, Is.EqualTo((rotation.Cosine * a - rotation.Sine * b).Real).Within(1e-12));
                Assert.That(c.Imaginary, Is.EqualTo((rotation.Cosine * a - rotation.Sine * b).Imaginary).Within(1e-12));

                Assert.That(d.Magnitude, Is.EqualTo(0).Within(1e-12));
                Assert.That(d.Real, Is.EqualTo((rotation.Cosine.Conjugate() * b + rotation.Sine * a).Real).Within(1e-12));
                Assert.That(d.Imaginary, Is.EqualTo((rotation.Cosine.Conjugate() * b + rotation.Sine * a).Imaginary).Within(1e-12));
            });
        }

        [TestCase(1, 2, 3)]
        [TestCase(-0.40166, -0.073029, 0)]
        [TestCase(-0.40166, 0.073029, 0.91287)]
        [TestCase(0.0123, 0.5, 0)]
        [TestCase(-0.0123, 0.5, 0)]
        [TestCase(0.0123, -0.5, 0)]
        [TestCase(-0.0123, -0.5, 0)]
        [TestCase(2, 3, 1)]
        [TestCase(2, 3, -1)]
        [TestCase(-2, 3, -1)]
        [TestCase(2, -3, 1)]
        [TestCase(-2, 3, 1)]
        [TestCase(3, 1, 2)]
        [TestCase(-0.18257418583505544, -0.3651483716701108, -0.9128709291752769)]
        public void CalculateGivensRotationShouldBePositive(double ar, double ai, double br)
        {
            (var a, var b) = (new Complex(ar, ai), new Complex(br, 0));
            (var rotation, var norm) = AmvwSolver.GivensRotation.CreatePositive(a, b);
            (var c, var d) = rotation * (new Complex(ar, ai), new Complex(br, 0));

            Assert.That(c.Real, Is.GreaterThan(0));
            Assert.That(c.Imaginary, Is.EqualTo(0).Within(1e-12));
        }

        [Test]
        public void CalculateGivensRotationTest2()
        {
            var a = new Complex(1, 2);
            var b = new Complex(1, 0);
            var norm = Math.Sqrt(1 + 1 + 4);
            (var rotation, var alpha) = GivensRotation.Create(a, b);
            Assert.That(rotation.Cosine.Real, Is.EqualTo(0.40824829046386307).Within(1e-12));
            Assert.That(rotation.Cosine.Imaginary, Is.EqualTo(-0.81649658092772615).Within(1e-12));
            Assert.That(rotation.Sine, Is.EqualTo(-0.40824829046386307).Within(1e-12));
            Assert.That(alpha, Is.EqualTo(norm).Within(1e-12));

            (var c, var d) = rotation * (a, b);
            Assert.That(d.Magnitude, Is.EqualTo(0).Within(1e-12));
            Assert.That(c.Magnitude, Is.EqualTo(norm).Within(1e-12));
            Assert.That(c.Imaginary, Is.EqualTo(0).Within(1e-12));
            Assert.That(c.Real, Is.GreaterThan(0));
        }

        [Test]
        public void CreateUnitaryMatrix()
        {
            Complex[] vector = new Complex[]
            {
                new Complex(1, 2),
                new Complex(3, 4),
                new Complex(5, 6),
                new Complex(7, 8)
            };
            var copy = new Complex[vector.Length];
            Array.Copy(vector, copy, vector.Length);

            var norm = Math.Sqrt(vector.Select(x => x.MagnitudeSquared()).Sum());
            var matrix = AmvwSolver.UnitaryMatrix.Create(vector);
            Assert.That(matrix.Rotations.Length, Is.EqualTo(vector.Length - 1));

            // apply transformations to original vector
            for (var i = matrix.Rotations.Length - 1; i >= 0; --i)
            {
                var rotation = matrix.Rotations[i];
                (copy[i], copy[i + 1]) = rotation * (copy[i], copy[i + 1]);
            }

            Assert.Multiple(() =>
            {
                Assert.That(copy[0].Magnitude, Is.EqualTo(norm).Within(1e-12));
                Assert.That(copy[1].Magnitude, Is.EqualTo(0).Within(1e-12));
                Assert.That(copy[2].Magnitude, Is.EqualTo(0).Within(1e-12));
                Assert.That(copy[3].Magnitude, Is.EqualTo(0).Within(1e-12));
            });
        }


        [Test]
        public void CreateFactorizationTestQ()
        {
            var factorization = AmvwFactorization.Create(GetTestCoefficients());
            var matrix = factorization.Q.GetMatrix(4);
            Console.WriteLine(matrix);
            Assert.Multiple(() =>
            {
                Assert.That(matrix[0, 0], Is.EqualTo(Complex.Zero));
                Assert.That(matrix[0, 1], Is.EqualTo(Complex.Zero));

                Assert.That(matrix[1, 0], Is.EqualTo(Complex.One));
                Assert.That(matrix[1, 1], Is.EqualTo(Complex.Zero));
                Assert.That(matrix[2, 2], Is.EqualTo(Complex.Zero));

                Assert.That(matrix[2, 1], Is.EqualTo(Complex.One));

                Assert.That(matrix[0, 2], Is.EqualTo(Complex.One));

                Assert.That(matrix[3, 3], Is.EqualTo(Complex.One));
            });
        }

        /// <summary>
        /// Checks if the C transformations applied to the vector x is a multiple of e1*alpha,
        /// where alpha is the 2-norm of the vector x on page 946 of [1].
        /// </summary>
        [Test]
        public void CreateFactorizationTestC()
        {
            var coefficients = GetTestCoefficients();
            var factorization = AmvwFactorization.Create(coefficients);
            int n = coefficients.Length - 1;
            var x = Vector<Complex>.Build.DenseOfArray(AmvwFactorization.CreateVectorX(coefficients));
            x[n] = x[n].Magnitude;
            var diagonalCorrection = new Complex(0.447213595499958, 0.894427190999916);

            Assert.That(factorization.C.Rotations[0].Cosine.Real, Is.EqualTo(0.512877644532173).Within(1e-12));
            Assert.That(factorization.C.Rotations[0].Cosine.Imaginary, Is.EqualTo(0.093250480824031).Within(1e-12));
            Assert.That(factorization.C.Rotations[0].Sine, Is.EqualTo(-0.853382018538718).Within(1e-12));

            Assert.That(factorization.C.Rotations[1].Cosine.Real, Is.EqualTo(0.928809220003860).Within(1e-12));
            Assert.That(factorization.C.Rotations[1].Cosine.Imaginary, Is.EqualTo(0.218543345883261).Within(1e-12));
            Assert.That(factorization.C.Rotations[1].Sine, Is.EqualTo(-0.299252800832290).Within(1e-12));

            Assert.That(factorization.C.Rotations[2].Cosine.Real, Is.EqualTo(0.408248290463863).Within(1e-12));
            Assert.That(factorization.C.Rotations[2].Cosine.Imaginary, Is.EqualTo(0.816496580927726).Within(1e-12));
            Assert.That(factorization.C.Rotations[2].Sine, Is.EqualTo(0.408248290463863).Within(1e-12));
        }

        [Test]
        public void GetDiagonal()
        {
            var coefficients = GetTestCoefficients();
            var factorization = AmvwFactorization.Create(coefficients);
            var C = factorization.C.GetMatrix(4);
            Console.WriteLine(C);
            for (var i = 0; i < 4; ++i)
            {
                Assert.Multiple(() =>
                {
                    Assert.That(C[i, i].Real, Is.EqualTo(factorization.C.GetDiagonal(i).Real).Within(1e-12));
                    Assert.That(C[i, i].Imaginary, Is.EqualTo(factorization.C.GetDiagonal(i).Imaginary).Within(1e-12));
                });
            }
        }


        [Test]
        public void GetSuperDiagonal()
        {
            var coefficients = GetTestCoefficients();
            var factorization = AmvwFactorization.Create(coefficients);
            var C = factorization.C.GetMatrix(5);
            Console.WriteLine(C);
            for (var i = 1; i < 5; ++i)
            {
                Assert.Multiple(() =>
                {
                    Assert.That(C[i - 1, i].Real, Is.EqualTo(factorization.C.GetSuperdiagonal(i).Real).Within(1e-12));
                    Assert.That(C[i - 1, i].Imaginary, Is.EqualTo(factorization.C.GetSuperdiagonal(i).Imaginary).Within(1e-12));
                });
            }
        }

        /// <summary>
        /// Tests the conversion of a rotation R when a matrix
        /// Z = [0, 1;
        ///      1, 0]
        ///  is applied to it. The new rotation should satisfy
        ///  Z * rotation = diag([d, e]) * newRotation
        /// </summary>
        [TestCase(1, 2, 3)]
        [TestCase(1, -2, 3)]
        [TestCase(1, -2, -3)]
        [TestCase(-1, -2, -3)]
        [TestCase(-1, -2, 3)]
        public void GetReflectedRotation(double cr, double ci, double sine)
        {
            var rotation = new GivensRotation(new Complex(cr, ci), sine);
            rotation.Normalize();

            var r1 = rotation.GetMatrix(2, 0);
            var Z = Matrix<Complex>.Build.Dense(2, 2);
            Z[0, 1] = 1;
            Z[1, 0] = 1;

            var expectedMatrix = r1 * Z;

            (var d, var e, var newRotation) = GivensRotation.GetReflectedRotation(rotation);

            var r2 = newRotation.GetMatrix(2, 0);
            var diagonal = Matrix<Complex>.Build.Sparse(2, 2);
            diagonal[0, 0] = d;
            diagonal[1, 1] = e;

            var actualMatrix = diagonal * r2;
            Assert.That((expectedMatrix - actualMatrix).FrobeniusNorm(), Is.LessThan(1e-12));
        }

        [Test]
        public void CreateFactorizationTestB()
        {
            var coefficients = GetTestCoefficients();
            var factorization = AmvwFactorization.Create(coefficients);

            Assert.That(factorization.B.Rotations[0].Cosine.Real, Is.EqualTo(0.512877644532173).Within(1e-12));
            Assert.That(factorization.B.Rotations[0].Cosine.Imaginary, Is.EqualTo(0.093250480824031).Within(1e-12));
            Assert.That(factorization.B.Rotations[0].Sine, Is.EqualTo(-0.853382018538718).Within(1e-12));

            Assert.That(factorization.B.Rotations[1].Cosine.Real, Is.EqualTo(0.928809220003860).Within(1e-12));
            Assert.That(factorization.B.Rotations[1].Cosine.Imaginary, Is.EqualTo(0.218543345883261).Within(1e-12));
            Assert.That(factorization.B.Rotations[1].Sine, Is.EqualTo(-0.299252800832290).Within(1e-12));

            Assert.That(factorization.B.Rotations[2].Cosine.Real, Is.EqualTo(0.182574185835055).Within(1e-12));
            Assert.That(factorization.B.Rotations[2].Cosine.Imaginary, Is.EqualTo(0.365148371670111).Within(1e-12));
            Assert.That(factorization.B.Rotations[2].Sine, Is.EqualTo(-0.912870929175277).Within(1e-12));
        }

        [Test]
        public void GetTopLeftR()
        {
            var coefficients = GetTestCoefficients();
            var factorization = AmvwFactorization.Create(coefficients);
            (var a, var b, var c, var d) = factorization.GetTopLeftR();
            Assert.Multiple(() =>
            {
                Assert.That(d, Is.EqualTo(Complex.One).Within(1e-12));
                Assert.That(a, Is.EqualTo(Complex.One).Within(1e-12));
                Assert.That(b, Is.EqualTo(Complex.Zero).Within(1e-12));
                Assert.That(c, Is.EqualTo(Complex.Zero).Within(1e-12));
            });
        }

        [Test]
        public void GetRBlock()
        {
            var coefficients = GetTestCoefficients();
            var factorization = AmvwFactorization.Create(coefficients);
            (Complex R00, Complex R01, Complex R10, Complex R11, Complex R20, Complex R21) = factorization.GetRBlock(1);

            Assert.That(R00, Is.EqualTo(Complex.Zero).Within(1e-12));
            Assert.That(R10, Is.EqualTo(Complex.One).Within(1e-12));
            Assert.That(R20, Is.EqualTo(Complex.Zero).Within(1e-12));

            // coeffienct[0] is 1+2j
            Assert.That(R21.Real, Is.EqualTo(-1).Within(1e-8));
            Assert.That(R21.Imaginary, Is.EqualTo(-2).Within(1e-8));
            Assert.That(R11.Real, Is.EqualTo(-5).Within(1e-1));
            Assert.That(R11.Imaginary, Is.EqualTo(-6).Within(1e-8));
            Assert.That(R01.Real, Is.EqualTo(-3).Within(1e-1));
            Assert.That(R01.Imaginary, Is.EqualTo(-4).Within(1e-8));
        }


        [Test]
        public void CreateFactorizationTestD()
        {
            var coefficients = GetTestCoefficients();
            var factorization = AmvwFactorization.Create(coefficients);
            Assert.Multiple(() =>
            {
                Assert.That(factorization.D[0].Real, Is.EqualTo(1).Within(1e-12));
                Assert.That(factorization.D[0].Imaginary, Is.EqualTo(0).Within(1e-12));
                Assert.That(factorization.D[1].Real, Is.EqualTo(1).Within(1e-12));
                Assert.That(factorization.D[1].Imaginary, Is.EqualTo(0).Within(1e-12));

                Assert.That(factorization.D[2].Real, Is.EqualTo(0.447213595499958).Within(1e-12));
                Assert.That(factorization.D[2].Imaginary, Is.EqualTo(0.894427190999916).Within(1e-12));

                Assert.That(factorization.D[3].Real, Is.EqualTo(-0.447213595499958).Within(1e-12));
                Assert.That(factorization.D[3].Imaginary, Is.EqualTo(0.894427190999916).Within(1e-12));
            });
        }

        [Test]
        public void CreateFactorization()
        {
            // A is a 3x3 matrix.
            // [0   0   -1-2j;
            //  1   0   -3-4j;
            //  0   1   -5-6j]
            var coefficients = GetTestCoefficients();
            var factorization = AmvwFactorization.Create(coefficients);

            var xNorm = Math.Sqrt(1 + 2 * 2 + 3 * 3 + 4 * 4 + 5 * 5 * 6 * 6);

            var R00 = factorization.GetDiagonalOfR(0);
            Assert.That(R00.Real, Is.EqualTo(1).Within(1e-12));
            Assert.That(R00.Imaginary, Is.EqualTo(0).Within(1e-12));

            var R11 = factorization.GetDiagonalOfR(1);
            Assert.That(R11.Real, Is.EqualTo(1).Within(1e-12));
            Assert.That(R11.Imaginary, Is.EqualTo(0).Within(1e-12));

            var R22 = factorization.D[2] * factorization.GetDiagonalOfR(2);
            Assert.That(R22.Real, Is.EqualTo(-coefficients[0].Real).Within(1e-12));
            Assert.That(R22.Imaginary, Is.EqualTo(-coefficients[0].Imaginary).Within(1e-12));

            // The factorization should be equal to the companion matrix.
            (var A00, var A01, var A10, var A11) = factorization.GetTopLeft();
            Assert.Multiple(() =>
            {
                Assert.That(A00.Magnitude, Is.EqualTo(0).Within(1e-12));
                Assert.That(A10.Real, Is.EqualTo(1).Within(1e-12));
                Assert.That(A10.Imaginary, Is.EqualTo(0).Within(1e-12));
                Assert.That(A01.Magnitude, Is.EqualTo(0).Within(1e-12));
                Assert.That(A11.Magnitude, Is.EqualTo(0).Within(1e-12));
            });
            (var A11_, var A12, var A21, var A22) = factorization.GetPart(1);
            Assert.Multiple(() =>
            {
                Assert.That(A11_.Magnitude, Is.EqualTo(0).Within(1e-12));
                Assert.That(A21.Real, Is.EqualTo(1).Within(1e-12));
                Assert.That(A21.Imaginary, Is.EqualTo(0).Within(1e-12));
                Assert.That(A12.Real, Is.EqualTo(-3).Within(1e-12));
                Assert.That(A12.Imaginary, Is.EqualTo(-4).Within(1e-12));
                Assert.That(A22.Real, Is.EqualTo(-5).Within(1e-12));
                Assert.That(A22.Imaginary, Is.EqualTo(-6).Within(1e-12));
            });
            int i = 2;
            var test = factorization.D[i - 1] * (factorization.GetDiagonalOfH(i)) / factorization.C.Rotations[i - 1].Sine;
        }

        [Test]
        public void CreateFactorization2()
        {
            // A is a 3x3 matrix.
            // [0   0   0 -1-2j;
            //  1   0   0 -3-4j;
            //  0   1   0 -5-6j;
            //  0   0   1 -7-8j]
            var testCoefficients = new Complex[]
            {
                new Complex(1, 2),
                new Complex(3, 4),
                new Complex(5, 6),
                new Complex(7, 8)
            };
            var factorization = AmvwFactorization.Create(testCoefficients);

            // The factorization should be equal to the companion matrix.
            (var A00, var A01, var A10, var A11) = factorization.GetTopLeft();
            Assert.Multiple(() =>
            {
                Assert.That(A00.Magnitude, Is.EqualTo(0).Within(1e-12));
                Assert.That(A10.Real, Is.EqualTo(1).Within(1e-12));
                Assert.That(A10.Imaginary, Is.EqualTo(0).Within(1e-12));
                Assert.That(A01.Magnitude, Is.EqualTo(0).Within(1e-12));
                Assert.That(A11.Magnitude, Is.EqualTo(0).Within(1e-12));
            });
            (var A11_, var A12, var A21, var A22) = factorization.GetPart(1);
            Assert.Multiple(() =>
            {
                Assert.That(A11_.Magnitude, Is.EqualTo(0).Within(1e-12));
                Assert.That(A21.Real, Is.EqualTo(1).Within(1e-12));
                Assert.That(A21.Imaginary, Is.EqualTo(0).Within(1e-12));
                Assert.That(A12.Real, Is.EqualTo(0).Within(1e-12));
                Assert.That(A12.Imaginary, Is.EqualTo(0).Within(1e-12));
                Assert.That(A22.Real, Is.EqualTo(0).Within(1e-12));
                Assert.That(A22.Imaginary, Is.EqualTo(0).Within(1e-12));
            });
            (var A22_, var A23, var A32, var A33) = factorization.GetPart(2);
            Assert.That(A22_.Magnitude, Is.EqualTo(0).Within(1e-12));
            Assert.That(A23.Real, Is.EqualTo(-testCoefficients[2].Real).Within(1e-12));
            Assert.That(A23.Imaginary, Is.EqualTo(-testCoefficients[2].Imaginary).Within(1e-12));
            Assert.That(A32.Real, Is.EqualTo(1).Within(1e-12));
            Assert.That(A32.Imaginary, Is.EqualTo(0).Within(1e-12));
            Assert.That(A33.Real, Is.EqualTo(-testCoefficients[3].Real).Within(1e-12));
            Assert.That(A33.Imaginary, Is.EqualTo(-testCoefficients[3].Imaginary).Within(1e-12));
        }

        [Test]
        public void CreateFactorizationCompareMatrix()
        {
            // A is a 3x3 matrix.
            // [0   0   0 -1-2j;
            //  1   0   0 -3-4j;
            //  0   1   0 -5-6j;
            //  0   0   1 -7-8j]
            var testCoefficients = new Complex[]
            {
                new Complex(1, 2),
                new Complex(3, 4),
                new Complex(5, 6),
                new Complex(6, 7)
            };
            int n = testCoefficients.Length + 1;
            var factorization = AmvwFactorization.Create(testCoefficients);
            var actualMatrix = factorization.GetMatrix().SubMatrix(0, n - 1, 0, n - 1);
            var expectedMatrix = CompanionMatrix(testCoefficients);
            Console.WriteLine(expectedMatrix);
            Console.WriteLine(actualMatrix);
            Assert.That((expectedMatrix - actualMatrix).FrobeniusNorm(), Is.LessThan(1e-12));
        }

        private Matrix<Complex> CompanionMatrix(Complex[] coefficients)
        {
            int n = coefficients.Length;
            var expectedMatrix = Matrix<Complex>.Build.Dense(n, n, 0);
            for (int i = 0; i < n - 1; ++i)
            {
                expectedMatrix[i + 1, i] = 1;
                expectedMatrix[i, n - 1] = -coefficients[i];
            }
            expectedMatrix[n - 1, n - 1] = -coefficients[n - 1];
            return expectedMatrix;
        }

        [Test]
        public void CreateFactorizationCompareMatrix2()
        {
            // A is a 3x3 matrix.
            // [0   0   -1-2j;
            //  1   0   -3-4j;
            //  0   1   -5-6j;]
            var testCoefficients = new Complex[]
            {
                new Complex(1, 2),
                new Complex(3, 4),
                new Complex(5, 6),
            };
            int n = testCoefficients.Length + 1;
            var factorization = AmvwFactorization.Create(testCoefficients);
            var actualMatrix = factorization.GetMatrix().SubMatrix(0, n - 1, 0, n - 1);
            var expectedMatrix = CompanionMatrix(testCoefficients);
            Assert.That((expectedMatrix - actualMatrix).FrobeniusNorm(), Is.LessThan(1e-12));
        }


        [Test]
        public void CreateFactorizationCheckRoots()
        {
            var coefficients = GetTestCoefficients();
            int n = coefficients.Length;
            var factorization = AmvwFactorization.Create(coefficients);

            var expectedRoots = GetTestSolution().Append(Complex.Zero);
            var actualRoots = factorization.GetMatrix().Evd().EigenValues;
            AssertEqualityWithTolerance(expectedRoots, actualRoots, 12);
        }
        private HashSet<Complex> GetRoundedSet(IEnumerable<Complex> set, int digits)
            => set.OrderBy(x => x.MagnitudeSquared()).ThenBy(x => x.Real).ThenBy(x => x.Imaginary)
                .Select(x => new Complex(x.Real.Round(digits), x.Imaginary.Round(digits))).ToHashSet();

        private void AssertEqualityWithTolerance(IEnumerable<Complex> expectedSet, IEnumerable<Complex> actualSet, int digits)
            => Assert.That(GetRoundedSet(actualSet, digits), Is.EquivalentTo(GetRoundedSet(expectedSet, digits)));


        [Test]
        public void ApplyFrancisStepIdentityCheckSimilarity()
        {
            var coefficients = GetTestCoefficients();
            var factorization = AmvwFactorization.Create(coefficients);
            int n = coefficients.Length + 1;
            var expectedMatrix = factorization.GetMatrix();
            var u1 = new GivensRotation(1, 0);
            AmvwSolver.ApplyFrancisStep(u1, factorization);

            // The process does not keep track of y, but on this example, y will be multiplied with -1.
            factorization.y *= -1;
            var actualMatrix = factorization.GetMatrix();

            // Check similarity
            var actualEigenValues = actualMatrix.Evd().EigenValues;
            var expectedEigenValues = expectedMatrix.Evd().EigenValues;

            AssertEqualityWithTolerance(expectedEigenValues, actualEigenValues, 12);
        }


        [Test]
        public void ApplyTransformationLeftOnFactorization()
        {
            var coefficients = GetTestCoefficients();
            var factorization = AmvwFactorization.Create(coefficients);
            int n = coefficients.Length;
            (var u1, _) = GivensRotation.Create(new Complex(1, 2), new Complex(3, 4));

            var expectedMatrix = u1.GetMatrix(n + 1, 0) * factorization.GetMatrix();

            factorization.ApplyTransformationLeft(u1);

            var actualMatrix = factorization.GetMatrix();
            Assert.That((expectedMatrix - actualMatrix).FrobeniusNorm(), Is.LessThan(1e-12));
        }


        [TestCase(0)]
        [TestCase(1)]
        public void PassthroughRow(int i)
        {
            var testCoefficients = new Complex[]
            {
                new Complex(1, 2),
                new Complex(3, 4),
                new Complex(5, 6),
                new Complex(6, 7)
            };
            var factorization = AmvwFactorization.Create(testCoefficients);
            int n = testCoefficients.Length;

            (var u1, _) = GivensRotation.Create(new Complex(1, 2), new Complex(3, 4));
            var expectedMatrix = factorization.GetMatrix() * u1.GetMatrix(n + 1, i);

            (var u2, var j) = factorization.Passthrough(i, u1);

            Assert.That(j, Is.EqualTo(i + 1));

            var actualMatrix = u2.GetMatrix(n + 1, j) * factorization.GetMatrix();
            Assert.That((expectedMatrix - actualMatrix).FrobeniusNorm(), Is.LessThan(1e-12));
        }

        [Test]
        public void PassthroughFirstRowCheckSimilarity()
        {
            var coefficients = GetTestCoefficients();
            var factorization = AmvwFactorization.Create(coefficients);
            var originalMatrix = factorization.GetMatrix();
            Console.WriteLine(originalMatrix);
            var originalEigenValues = originalMatrix.Evd().EigenValues;
            int n = coefficients.Length;

            (var u1, _) = GivensRotation.Create(new Complex(1, 2), new Complex(3, 4));
            u1.Normalize();
            var u1H = (GivensRotation)u1.Clone();
            u1H.Invert();

            factorization.ApplyTransformationLeft(u1H);

            (var u2, var j) = factorization.Passthrough(0, u1);
            // C = u2 * factorization = u1H * A * u1
            // u2^H * C * u2 = factorization * u2.

            var actualMatrix = u2.GetMatrix(n + 1, j) * factorization.GetMatrix();
            var actualEigenValues = actualMatrix.Evd().EigenValues;

            Console.WriteLine(originalEigenValues);
            Console.WriteLine(actualEigenValues);

        }

        [Test]
        public void PassthroughSecondRow()
        {
            var coefficients = GetTestCoefficients();
            var factorization = AmvwFactorization.Create(coefficients);
            int n = coefficients.Length;

            (var u2, _) = GivensRotation.Create(new Complex(1, 2), new Complex(3, 4));

            // only considering the part Q*D*C'*B because the transformation is not applied to y.
            var expectedMatrix = factorization.Q.GetMatrix(n + 1) *
                Matrix<Complex>.Build.DiagonalOfDiagonalArray(factorization.D) *
                factorization.C.GetMatrix(n + 1).ConjugateTranspose() *
                factorization.B.GetMatrix(n + 1) * u2.GetMatrix(n + 1, 1);

            (u2, var j) = factorization.Passthrough(1, u2);

            Assert.That(j, Is.EqualTo(1));
            Assert.That(u2.Cosine, Is.EqualTo(Complex.One));
            Assert.That(u2.Sine, Is.EqualTo(0));

            var actualMatrix = factorization.Q.GetMatrix(n + 1) *
                Matrix<Complex>.Build.DiagonalOfDiagonalArray(factorization.D) *
                factorization.C.GetMatrix(n + 1).ConjugateTranspose() *
                factorization.B.GetMatrix(n + 1);
            Assert.That((expectedMatrix - actualMatrix).FrobeniusNorm(), Is.LessThan(1e-12));
        }

        [Test]
        public void AbsorbTransformationInQ()
        {
            var coefficients = GetTestCoefficients();
            var factorization = AmvwFactorization.Create(coefficients);
            int n = coefficients.Length;

            // u2_cosine = 0.4016632088371218-0.07302967433402213
            // u2_sine = -0.9128709291752769
            (var u2, _) = GivensRotation.Create(new Complex(1, 2), new Complex(3, 4));
            int i = n - 2;

            var expectedMatrix = factorization.Q.GetMatrix(n + 1) * u2.GetMatrix(n + 1, i) *
                Matrix<Complex>.Build.DiagonalOfDiagonalArray(factorization.D) *
                factorization.C.GetMatrix(n + 1).ConjugateTranspose() *
                (factorization.B.GetMatrix(n + 1) + factorization.GetYMatrix());

            factorization.AbsorbTransformationInQ(u2, n - 2);

            var actualMatrix = factorization.Q.GetMatrix(n + 1) *
                Matrix<Complex>.Build.DiagonalOfDiagonalArray(factorization.D) *
                factorization.C.GetMatrix(n + 1).ConjugateTranspose() *
                (factorization.B.GetMatrix(n + 1) + factorization.GetYMatrix());
            Console.WriteLine(expectedMatrix);
            Console.WriteLine(actualMatrix);
            Assert.That((expectedMatrix - actualMatrix).FrobeniusNorm(), Is.LessThan(1e-12));
        }


        [TestCase(0)]
        [TestCase(2)]
        public void AbsorbTransformationInQShouldThrowException(int i)
        {
            var coefficients = GetTestCoefficients();
            var factorization = AmvwFactorization.Create(coefficients);
            int n = coefficients.Length;

            (var u2, _) = GivensRotation.Create(new Complex(1, 2), new Complex(3, 4));

            Assert.Throws<ArgumentException>(() => factorization.AbsorbTransformationInQ(u2, i));
        }

        [Test]
        public void CreateVectorX()
        {
            var testCoefficients = new Complex[]
            {
                new Complex(1, 2),
                new Complex(3, 4),
                new Complex(5, 6),
                new Complex(7, 8)
            };
            var x = AmvwFactorization.CreateVectorX(testCoefficients);
            Assert.Multiple(() =>
            {
                Assert.That(x[0], Is.EqualTo(-testCoefficients[1]));
                Assert.That(x[1], Is.EqualTo(-testCoefficients[2]));
                Assert.That(x[2], Is.EqualTo(-testCoefficients[3]));
                Assert.That(x[3], Is.EqualTo(-testCoefficients[0]));
                Assert.That(x[4], Is.EqualTo(-Complex.One));
            });
        }

        [Test]
        public void Turnover()
        {
            // rotation1 = ({(0,40166320883712187, -0,07302967433402213)}, 0.9128709291752769)
            (GivensRotation rotation1, _) = GivensRotation.Create(new Complex(1, 2), new Complex(-3, 4));
            // rotation2 = ({(0,5701081850582409, -0,019325701188414927)}, 0.82134230050763524)
            (GivensRotation rotation2, _) = GivensRotation.Create(new Complex(4, -5), new Complex(6, 7));
            // rotation3 = ({(0,619902067513713, -0,008669958986205754)}, 0.784631288251623)
            (GivensRotation rotation3, _) = GivensRotation.Create(new Complex(7, 8), new Complex(-9, 10));

            var expectedMatrix = rotation1.GetMatrix(3, 0) * rotation2.GetMatrix(3, 1) * rotation3.GetMatrix(3, 0);

            (var rotation4, var rotation5, var rotation6) = GivensRotation.Turnover(rotation1, rotation2, rotation3);

            var actualMatrix = rotation4.GetMatrix(3, 1) * rotation5.GetMatrix(3, 0) * rotation6.GetMatrix(3, 1);

            Assert.That((expectedMatrix - actualMatrix).FrobeniusNorm(), Is.LessThan(1e-12));
        }

        [Test]
        public void TurnoverIdentity()
        {
            (var rotation1, _) = GivensRotation.Create(new Complex(1, 2), new Complex(3, 4));
            (var rotation2, _) = GivensRotation.Create(new Complex(4, 5), new Complex(6, 7));
            var rotation3 = new GivensRotation(-Complex.One, 0);

            (var rotation4, var rotation5, var rotation6) = GivensRotation.Turnover(rotation1, rotation2, rotation3);

            var expectedMatrix = rotation1.GetMatrix(3, 0) * rotation2.GetMatrix(3, 1) * rotation3.GetMatrix(3, 0);
            var actualMatrix = rotation4.GetMatrix(3, 1) * rotation5.GetMatrix(3, 0) * rotation6.GetMatrix(3, 1);
            Console.WriteLine(expectedMatrix);
            Console.WriteLine(actualMatrix);
            Assert.That((expectedMatrix - actualMatrix).FrobeniusNorm(), Is.LessThan(1e-12));
        }

        [Test]
        public void TurnoverIdentity2()
        {
            var coefficients = GetTestCoefficients();
            var factorization = AmvwFactorization.Create(coefficients);
            var rotation1 = factorization.B.Rotations[0];
            rotation1.Normalize();
            var rotation2 = factorization.B.Rotations[1];
            rotation2.Normalize();
            var rotation3 = new GivensRotation(Complex.One, 0);

            (var rotation4, var rotation5, var rotation6) = GivensRotation.Turnover(rotation1, rotation2, rotation3);

            var expectedMatrix = rotation1.GetMatrix(3, 0) * rotation2.GetMatrix(3, 1) * rotation3.GetMatrix(3, 0);
            var actualMatrix = rotation4.GetMatrix(3, 1) * rotation5.GetMatrix(3, 0) * rotation6.GetMatrix(3, 1);
            Console.WriteLine(expectedMatrix);
            Console.WriteLine(actualMatrix);
            Assert.That((expectedMatrix - actualMatrix).FrobeniusNorm(), Is.LessThan(1e-12));
        }

        [Test]
        public void TurnoverIdentity3()
        {
            var rotation1 = new GivensRotation(new Complex(0.312771621085, -0.417028828), -0.853382018538);
            rotation1.Normalize();
            var rotation2 = new GivensRotation(new Complex(0.61084722178, -0.73301666613), -0.299252800);
            rotation2.Normalize();
            var rotation3 = new GivensRotation(Complex.One, 0);

            (var rotation4, var rotation5, var rotation6) = GivensRotation.Turnover(rotation1, rotation2, rotation3);

            var expectedMatrix = rotation1.GetMatrix(3, 0) * rotation2.GetMatrix(3, 1) * rotation3.GetMatrix(3, 0);
            var actualMatrix = rotation4.GetMatrix(3, 1) * rotation5.GetMatrix(3, 0) * rotation6.GetMatrix(3, 1);
            Console.WriteLine(expectedMatrix);
            Console.WriteLine(actualMatrix);
            Assert.That((expectedMatrix - actualMatrix).FrobeniusNorm(), Is.LessThan(1e-12));
        }

        /// <summary>
        /// Checks if the turnover operation is executed on the correct transformations of a unitary matrix
        /// described by core transformations.
        /// </summary>
        [TestCase(0)]
        [TestCase(1)]
        [TestCase(2)]
        [TestCase(3)]
        public void ApplyTransformation(int i)
        {
            (var transformation, _) = GivensRotation.Create(new Complex(1, 2), new Complex(3, 4));
            var B = UnitaryMatrix.Create(new Complex[]
            {
                new Complex(4, 5),
                new Complex(9, -12),
                new Complex(7, 3),
                new Complex(2, -3),
                new Complex(9, 12),
                new Complex(-8, 7)
            });
            int n = B.Rotations.Length + 2;
            var expectedMatrix = B.GetMatrix(n) * transformation.GetMatrix(n, i);

            (transformation, var j) = B.ApplyTransformationRight(transformation, i);

            // transformation should move down one row.
            Assert.That(j, Is.EqualTo(i + 1));

            var actualMatrix = transformation.GetMatrix(n, j) * B.GetMatrix(n);

            Console.WriteLine(expectedMatrix);
            Console.WriteLine(actualMatrix);
            Assert.That((expectedMatrix - actualMatrix).FrobeniusNorm(), Is.LessThan(1e-12));
        }

        /// <summary>
        /// Checks if the turnover operation is executed on the correct complex conjugate 
        /// transformations of a unitary matrix described by core transformations.
        /// </summary>
        [TestCase(1)]
        [TestCase(2)]
        [TestCase(3)]
        [TestCase(4)]
        public void ApplyTransformationComplexConjugate(int i)
        {
            (var transformation, _) = GivensRotation.Create(new Complex(1, -2), new Complex(-3, 4));
            var B = UnitaryMatrix.Create(new Complex[]
            {
                new Complex(4, 5),
                new Complex(9, -12),
                new Complex(7, 3),
                new Complex(2, -3),
                new Complex(9, 12),
                new Complex(-8, 7)
            });
            int n = B.Rotations.Length + 2;
            var expectedMatrix = B.GetMatrix(n).ConjugateTranspose() * transformation.GetMatrix(n, i);

            (transformation, var j) = B.ApplyTransformationConjugate(transformation, i);

            // transformation should move down one row.
            Assert.That(j, Is.EqualTo(i - 1));

            var actualMatrix = transformation.GetMatrix(n, j) * B.GetMatrix(n).ConjugateTranspose();

            Console.WriteLine(expectedMatrix);
            Console.WriteLine(actualMatrix);
            Assert.That((expectedMatrix - actualMatrix).FrobeniusNorm(), Is.LessThan(1e-12));
        }

        [Test]
        public void ApplyTransformationShouldRaiseException()
        {
            var coefficients = GetTestCoefficients();
            var B = UnitaryMatrix.Create(coefficients);
            Assert.Throws<ArgumentException>(() => B.ApplyTransformationRight(new GivensRotation(), coefficients.Length - 1));
        }


        [Test]
        public void TurnoverComplexConjugate()
        {
            var rotation1 = new GivensRotation(new Complex(0.312771621085, -0.417028828), -0.853382018538);
            rotation1.Normalize();
            var rotation2 = new GivensRotation(new Complex(0.61084722178, -0.73301666613), -0.299252800);
            rotation2.Normalize();
            (var rotation3, _) = GivensRotation.Create(new Complex(1, 2), new Complex(3, 4));

            // Executing turnover on rotation1', rotation2' and rotation3.
            rotation1.Invert();
            rotation2.Invert();

            (var rotation4, var rotation5, var rotation6) = GivensRotation.Turnover(rotation1, rotation2, rotation3);

            rotation1.Invert();
            rotation2.Invert();
            rotation5.Invert();
            rotation6.Invert();

            var expectedMatrix = rotation4.GetMatrix(3, 0) * rotation5.GetMatrix(3, 1).ConjugateTranspose() * rotation6.GetMatrix(3, 0).ConjugateTranspose();

            (rotation4, rotation5, rotation6) = GivensRotation.TurnoverComplexConjugate(rotation1, rotation2, rotation3);

            var actualMatrix = rotation4.GetMatrix(3, 0) * rotation5.GetMatrix(3, 1).ConjugateTranspose() * rotation6.GetMatrix(3, 0).ConjugateTranspose();
            Console.WriteLine(expectedMatrix);
            Console.WriteLine(actualMatrix);
            Assert.That((expectedMatrix - actualMatrix).FrobeniusNorm(), Is.LessThan(1e-12));
        }

        [Test]
        public void TurnoverComplexConjugateIdentity()
        {
            var rotation1 = new GivensRotation(new Complex(0.312771621085, -0.417028828), -0.853382018538);
            rotation1.Normalize();
            var rotation2 = new GivensRotation(new Complex(0.61084722178, -0.73301666613), -0.299252800);
            rotation2.Normalize();
            var rotation3 = new GivensRotation(1, 0);

            // Executing turnover on rotation1', rotation2' and rotation3.
            rotation1.Invert();
            rotation2.Invert();

            (var rotation4, var rotation5, var rotation6) = GivensRotation.Turnover(rotation1, rotation2, rotation3);

            rotation1.Invert();
            rotation2.Invert();
            rotation5.Invert();
            rotation6.Invert();

            var expectedMatrix = rotation4.GetMatrix(3, 0) * rotation5.GetMatrix(3, 1).ConjugateTranspose() * rotation6.GetMatrix(3, 0).ConjugateTranspose();

            (rotation4, rotation5, rotation6) = GivensRotation.TurnoverComplexConjugate(rotation1, rotation2, rotation3);

            var actualMatrix = rotation4.GetMatrix(3, 0) * rotation5.GetMatrix(3, 1).ConjugateTranspose() * rotation6.GetMatrix(3, 0).ConjugateTranspose();
            Console.WriteLine(expectedMatrix);
            Console.WriteLine(actualMatrix);
            Assert.That((expectedMatrix - actualMatrix).FrobeniusNorm(), Is.LessThan(1e-12));
        }


        [Test]
        public void ApplyTransformationLeftReal()
        {
            var u1 = new GivensRotation(new Complex(2, 0), 5);
            u1.Normalize();
            var matrix = UnitaryMatrix.Create(new Complex[] { 1, 2, 3, 4 });
            int n = 4;
            int i = 0;

            var expectedMatrix = u1.GetMatrix(n, i) * matrix.GetMatrix(n);
            var diagonal = matrix.ApplyTransformationLeft(u1, i);
            var actualMatrix = matrix.GetMatrix(n) * Matrix<Complex>.Build.DiagonalOfDiagonalArray(diagonal);

            Assert.That((expectedMatrix - actualMatrix).FrobeniusNorm(), Is.LessThan(1e-12));

        }

        [TestCase(0)]
        //[TestCase(1)] would require a turnover.
        public void ApplyTransformationLeft(int i)
        {
            var coefficients = GetTestCoefficients();
            (var u1, _) = GivensRotation.Create(new Complex(1, 2), new Complex(1, 4));
            var matrix = UnitaryMatrix.Create(coefficients);
            var n = coefficients.Length;

            var expectedMatrix = u1.GetMatrix(n, i) * matrix.GetMatrix(n);
            var diagonal = matrix.ApplyTransformationLeft(u1, i);
            var actualMatrix = matrix.GetMatrix(n) * Matrix<Complex>.Build.DiagonalOfDiagonalArray(diagonal);

            Console.WriteLine(expectedMatrix);
            Console.WriteLine(actualMatrix);

            Assert.That((expectedMatrix - actualMatrix).FrobeniusNorm(), Is.LessThan(1e-12));

        }

        [Test]
        public void Fusion()
        {
            (var g1, _) = GivensRotation.Create(new Complex(1, -2), new Complex(3, -4));
            (var g2, _) = GivensRotation.Create(new Complex(-5, 6), new Complex(-7, 8));

            var expectedMatrix = g1.GetMatrix(2, 0) * g2.GetMatrix(2, 0);
            (var d, var e) = g1.Fusion(g2);
            var actualMatrix = g1.GetMatrix(2, 0) * Matrix<Complex>.Build.DiagonalOfDiagonalArray(new Complex[] { d, e });

            Assert.That((expectedMatrix - actualMatrix).FrobeniusNorm(), Is.LessThan(1e-12));
        }

        [Test]
        public void Fusion2()
        {
            (var g1, _) = GivensRotation.Create(new Complex(1, -2), new Complex(3, -4));
            (var g2, _) = GivensRotation.Create(new Complex(-5, 6), new Complex(-7, 8));

            var expectedMatrix = g1.GetMatrix(2, 0) * g2.GetMatrix(2, 0);
            (var d, var e) = g1.Fusion(g2);
            var actualMatrix = g1.GetMatrix(2, 0) * Matrix<Complex>.Build.DiagonalOfDiagonalArray(new Complex[] { d, e });

            Assert.That((expectedMatrix - actualMatrix).FrobeniusNorm(), Is.LessThan(1e-12));
        }

        [Test]
        public void FusionOnUnitaryMatrix()
        {
            var coefficients = GetTestCoefficients();
            var factorization = AmvwFactorization.Create(coefficients);
            var transformation = new GivensRotation(new Complex(-0.6518827942845281, -0.7094018643684572), -0.2679511473129994);
            var expectedMatrix = factorization.Q.GetMatrix(4) * transformation.GetMatrix(4, 1);
            (var d, var e) = factorization.Q.Rotations[1].Fusion(transformation);

            var actualMatrix = factorization.Q.GetMatrix(4) * Matrix<Complex>.Build.SparseOfDiagonalArray(new Complex[]
            {
                1, d, e, 1
            });

            Assert.That((expectedMatrix - actualMatrix).FrobeniusNorm(), Is.LessThan(1e-12));
        }

        [Test]
        public void PassthroughDiagonalFromRight()
        {
            var d = new Complex(1, 2);
            var e = new Complex(3, 4);
            d /= d.Magnitude;
            e /= e.Magnitude;
            (var rotation, _) = GivensRotation.Create(new Complex(5, 4), new Complex(3, 1));
            var expectedMatrix = Matrix<Complex>.Build.DiagonalOfDiagonalArray(new Complex[] { d, e }) * rotation.GetMatrix(2, 0);

            (var passedRotation, var f, var g) = GivensRotation.PassthroughRotationToLeft(d, e, rotation);

            // checking passedRotation * diag([f, g]) == diag([d, e]) * rotation
            var actualMatrix = passedRotation.GetMatrix(2, 0) * Matrix<Complex>.Build.DiagonalOfDiagonalArray(new Complex[] { f, g });
            Assert.That((expectedMatrix - actualMatrix).FrobeniusNorm(), Is.LessThan(1e-12));
        }

        [Test]
        public void PassthroughDiagonalFromLeft()
        {
            var d = new Complex(0.447213595499958, -0.894427190999916);
            var e = new Complex(-0.447213595499958, -0.894427190999916);
            (var rotation, _) = GivensRotation.Create(new Complex(5, 4), new Complex(3, 1));
            var expectedMatrix = rotation.GetMatrix(2, 0) * Matrix<Complex>.Build.DiagonalOfDiagonalArray(new Complex[] { d, e });

            (var passedRotation, var f, var g) = GivensRotation.PassthroughRotationToLeft(d, e, rotation);

            // checking passedRotation * diag([f, g]) == diag([d, e]) * rotation
            var actualMatrix = Matrix<Complex>.Build.DiagonalOfDiagonalArray(new Complex[] { f, g }) * passedRotation.GetMatrix(2, 0);
            Assert.That((expectedMatrix - actualMatrix).FrobeniusNorm(), Is.LessThan(1e-12));
        }

        [Test]
        public void PassthroughDiagonalFromLeftHermitian()
        {
            var d = new Complex(0.447213595499958, -0.894427190999916);
            var e = new Complex(-0.447213595499958, -0.894427190999916);
            (var rotation, _) = GivensRotation.Create(new Complex(5, 4), new Complex(3, 1));
            var expectedMatrix = rotation.GetMatrix(2, 0).ConjugateTranspose() * Matrix<Complex>.Build.DiagonalOfDiagonalArray(new Complex[] { d, e });

            (var passedRotation, var f, var g) = GivensRotation.PassthroughRotationToLeftHermitian(d, e, rotation);

            var actualMatrix = Matrix<Complex>.Build.DiagonalOfDiagonalArray(new Complex[] { f, g }) * passedRotation.GetMatrix(2, 0).ConjugateTranspose();
            Assert.That((expectedMatrix - actualMatrix).FrobeniusNorm(), Is.LessThan(1e-12));
        }

        [Test]
        public void PassthroughDiagonalToLeft()
        {
            var diagonal = new Complex[] {
                new Complex(0.447213595499958, -0.894427190999916),
                new Complex(0.599999877331049, -0.800000092001698),
                new Complex(0.894427190999916, 0.447213595499958),
                new Complex(-0.6, 0.8)
            };
            var testCoefficients = new Complex[] { new Complex(1, 2), new Complex(3, 4), new Complex(7, 9), new Complex(9, 1) };
            int n = testCoefficients.Length;
            var matrix = UnitaryMatrix.Create(testCoefficients);
            var expectedProduct = matrix.GetMatrix(n) * Matrix<Complex>.Build.DiagonalOfDiagonalArray(diagonal);

            var newDiagonal = matrix.PassthroughDiagonalToLeft(diagonal);
            var actualProduct = Matrix<Complex>.Build.DiagonalOfDiagonalArray(newDiagonal) * matrix.GetMatrix(n);

            Console.WriteLine(expectedProduct);
            Console.WriteLine(actualProduct);
            Assert.That((expectedProduct - actualProduct).FrobeniusNorm(), Is.LessThan(1e-12));
        }


        [Test]
        public void PassthroughDiagonalToRight()
        {
            var diagonal = new Complex[] {
                new Complex(0.447213595499958, -0.894427190999916),
                new Complex(0.447213595499958, 0.894427190999916),
                new Complex(0.599999877331049, -0.800000092001698),
                new Complex(0.443015775381405, 0.896513816270119)
            };
            var testCoefficients = new Complex[] { new Complex(1, 2), new Complex(3, 4), new Complex(7, 9), new Complex(9, 1) };
            int n = testCoefficients.Length;
            var matrix = UnitaryMatrix.Create(testCoefficients);
            var expectedProduct = Matrix<Complex>.Build.DiagonalOfDiagonalArray(diagonal) * matrix.GetMatrix(n);

            var newDiagonal = matrix.PassthroughDiagonalToRight(diagonal, 0);
            var actualProduct = matrix.GetMatrix(n) * Matrix<Complex>.Build.DiagonalOfDiagonalArray(newDiagonal);

            Console.WriteLine(expectedProduct);
            Console.WriteLine(actualProduct);
            Assert.That((expectedProduct - actualProduct).FrobeniusNorm(), Is.LessThan(1e-12));
        }


        [Test]
        public void PassthroughDiagonalHermitian()
        {
            var diagonal = new Complex[] {
                new Complex(0.447213595499958, -0.894427190999916),
                new Complex(0.447213595499958, 0.894427190999916),
                new Complex(0.599999877331049, -0.800000092001698),
                new Complex(0.443015775381405, 0.896513816270119)
            };
            var testCoefficients = new Complex[] { new Complex(1, 2), new Complex(3, 4), new Complex(7, 9), new Complex(9, 1) };
            int n = testCoefficients.Length;
            var matrix = UnitaryMatrix.Create(testCoefficients);
            var expectedProduct = matrix.GetMatrix(n).ConjugateTranspose() * Matrix<Complex>.Build.DiagonalOfDiagonalArray(diagonal);

            var newDiagonal = matrix.PassthroughDiagonalHermitian(diagonal);
            var actualProduct = Matrix<Complex>.Build.DiagonalOfDiagonalArray(newDiagonal) * matrix.GetMatrix(n).ConjugateTranspose();

            Console.WriteLine(expectedProduct);
            Console.WriteLine(actualProduct);
            Assert.That((expectedProduct - actualProduct).FrobeniusNorm(), Is.LessThan(1e-12));
        }

        [Test]
        public void ScaleCoefficients()
        {
            var alpha = new Complex(4, 0);
            // 3 + 2x + 5x^2 + x^3 = (x + 4.7106)*(x+0.1449+0.78581i)*(x+0.1449-0.78581i)
            var coefficients = new Complex[] { 3, 2, 5 };
            int n = 3;
            var scaledCoefficients = AmvwSolver.ScaleCoefficients(coefficients, alpha);
            
            // should be: 3/alpha^n + 2x/alpha^(n-1) + 5x^2/alpha^(n-2) + x^3/(alpha^(n-3))

            Assert.That(scaledCoefficients[0], Is.EqualTo(coefficients[0] / alpha.Power(n)));
            Assert.That(scaledCoefficients[1], Is.EqualTo(coefficients[1] / alpha.Power(n - 1)));
            Assert.That(scaledCoefficients[2], Is.EqualTo(coefficients[2] / alpha.Power(n - 2)));
        }

        [Test]
        public void Solve3()
        {
            var coefficients = GetTestCoefficients();
            var expectedRoots = GetTestSolution();
            var actualRoots = AmvwSolver.Solve(coefficients, 1e-4, ShiftStrategy.Wilkinson);
            AssertEqualityWithTolerance(expectedRoots, actualRoots.Item1, 1);
        }


        [Test]
        public void Solve()
        {
            // (x - (1+ 2i))(x - (1- 2i))(x - 5) = -25 + 15 x - 7 x^2 + x^3
            var testCoefficients = new Complex[] { -25, 15, -7 }; // leading coeffient is not provided.
            var expectedRoots = new Complex[] { 5, new Complex(1, 2), new Complex(1, -2) };
            var factorization = AmvwFactorization.Create(testCoefficients);
            (var actualRoots, _) = AmvwSolver.Solve(factorization, 1e-12);

            AssertEqualityWithTolerance(expectedRoots, actualRoots, 8);
        }


        [Test]
        public void Solve2()
        {
            // (x-5)*(x-10)*(x-15)
            var testCoefficients = new Complex[] { -750, 275, -30 }; // leading coeffient is not provided.
            var expectedRoots = new Complex[] { 5, 10, 15 };

            var actualRoots = AmvwSolver.Solve(testCoefficients, 1e-8, ShiftStrategy.Wilkinson);
            AssertEqualityWithTolerance(expectedRoots, actualRoots.Item1, 6);
        }

        [Test]
        public void Solve4()
        {
            // (x-5)*(x-10)*(x-15)*(x-20)
            // = x^4 - 50 x^3 + 875 x^2 - 6250 x + 15000
            var testCoefficients = new Complex[] { 15000, -6250, 875, -50 }; // leading coeffient is not provided.
            var expectedRoots = new Complex[] { 5, 10, 15, 20 };

            var actualRoots = AmvwSolver.Solve(testCoefficients, 1e-10, ShiftStrategy.Wilkinson);
            AssertEqualityWithTolerance(expectedRoots, actualRoots.Item1, 2);
        }
    }
}
