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
using NUnit.Framework;

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
        /// <summary>
        /// Checks if the method GivensRotation.Create calculates a rotation,
        /// which zeros the second component of the vector (1+2j, 3+4j).
        /// </summary>
        [TestCase(1, 2, 3, 4)]
        [TestCase(1, 2, 0, 0)]
        public void CalculateGivensRotationTest(double ar, double ai, double br, double bi)
        {
            (var a, var b) = (new Complex(ar, ai), new Complex(br, bi));
            (var rotation, var norm) = AmvwSolver.GivensRotation.Create(ref a, ref b);
            Assert.That(norm, Is.EqualTo(Math.Sqrt(ar*ar + ai*ai + br*br + bi*bi)));
            (var c, var d) = rotation * (new Complex(ar, ai), new Complex(br, bi));
            Assert.Multiple(() =>
            {
                Assert.That(c.Magnitude, Is.EqualTo(norm).Within(1e-12));
                Assert.That(a.Magnitude, Is.EqualTo(norm).Within(1e-12));
                Assert.That(d.Magnitude, Is.EqualTo(0).Within(1e-12));
                Assert.That(b.Magnitude, Is.EqualTo(0).Within(1e-12));
            });
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
            
            // apply transformations to original vector
            for(var i = matrix.CoreTransformations.Length-1; i >= 0; --i)
            {
                var rotation = matrix.CoreTransformations[i];
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
        public void CreateFactorization()
        {
            var coefficients = new Complex[]
            {
                new Complex(1, 2),
                new Complex(3, 4),
                new Complex(5, 6)
            };
            var factorization = AmvwFactorization.Create(coefficients);

            // The factorization should be equal to the companion matrix.
            (var a, var b, var c, var d) = factorization.GetPart(0);
            Assert.Multiple(() =>
            {
                Assert.That(a.Magnitude, Is.EqualTo(0).Within(1e-12));
                Assert.That(c.Real, Is.EqualTo(1).Within(1e-12));
                Assert.That(c.Imaginary, Is.EqualTo(0).Within(1e-12));
                Assert.That(b.Magnitude, Is.EqualTo(0).Within(1e-12));
                Assert.That(d.Magnitude, Is.EqualTo(0).Within(1e-12));
            });
        }

        [Test]
        public void CreateVectorX()
        {
            var coefficients = new Complex[]
            {
                new Complex(1, 2),
                new Complex(3, 4),
                new Complex(5, 6),
                new Complex(7, 8)
            };
            var x = AmvwFactorization.CreateVectorX(coefficients);
            Assert.Multiple(() =>
            {
                Assert.That(x[0], Is.EqualTo(-coefficients[1]));
                Assert.That(x[1], Is.EqualTo(-coefficients[2]));
                Assert.That(x[2], Is.EqualTo(-coefficients[3]));
                Assert.That(x[3], Is.EqualTo(-coefficients[0]));
                Assert.That(x[4], Is.EqualTo(Complex.One));
            });
        }

        [Test]
        public void Turnover()
        {
            (GivensRotation rotation1, _) = GivensRotation.Create(new Complex(1, 2), new Complex(3, 4));
            (GivensRotation rotation2, _) = GivensRotation.Create(new Complex(4, 5), new Complex(6, 7));
            (GivensRotation rotation3, _) = GivensRotation.Create(new Complex(7, 8), new Complex(9, 10));
            var a = new Complex(3, 5);
            var b = new Complex(7, 2);
            var c = new Complex(9, 12);

            Complex d = a;
            Complex e = b;
            Complex f = c;
            
            (d, e) = rotation3 * (d, e);
            (e, f) = rotation2 * (e, f);
            (d, e) = rotation1 * (d, e);

            (var rotation4, var rotation5, var rotation6) = GivensRotation.Turnover(rotation1, rotation2, rotation3);

            (b, c) = rotation6 * (b, c);
            (a, b) = rotation5 * (a, b);
            (b, c) = rotation4 * (b, c);

            Assert.AreEqual(d.Real, a.Real, 1e-12);
            Assert.AreEqual(d.Imaginary, a.Imaginary, 1e-12);
        }
    }
}
