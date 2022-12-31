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
            Assert.That(norm, Is.EqualTo(Math.Sqrt(ar*ar + ai*ai + br*br + bi*bi)));
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
            for(var i = matrix.Rotations.Length-1; i >= 0; --i)
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
        public void CreateFactorization()
        {
            var coefficients = new Complex[]
            {
                new Complex(1, 2),
                new Complex(3, 4),
                new Complex(5, 6)
            };
            // A is a 3x3 matrix.
            var factorization = AmvwFactorization.Create(coefficients);

            // The factorization should be equal to the companion matrix.
            (var A00, var A01, var A10, var A11) = factorization.GetPart(0);
                Assert.That(A00.Magnitude, Is.EqualTo(0).Within(1e-12));
                Assert.That(A10.Real, Is.EqualTo(1).Within(1e-12));
                Assert.That(A10.Imaginary, Is.EqualTo(0).Within(1e-12));
                Assert.That(A01.Magnitude, Is.EqualTo(0).Within(1e-12));
                Assert.That(A11.Magnitude, Is.EqualTo(0).Within(1e-12));

            (var A11_, var A12, var A21, var A22) = factorization.GetPart(1);
                Assert.That(A11_.Magnitude, Is.EqualTo(0).Within(1e-12));
                Assert.That(A21.Real, Is.EqualTo(1).Within(1e-12));
                Assert.That(A21.Imaginary, Is.EqualTo(0).Within(1e-12));
                Assert.That(A12.Magnitude, Is.EqualTo(0).Within(1e-12));
                Assert.That(A22.Magnitude, Is.EqualTo(0).Within(1e-12));
        }

        [Test]
        public void GetPart()
        {

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
            // rotation1 = ({(0,40166320883712187, -0,07302967433402213)}, 0.9128709291752769)
            (GivensRotation rotation1, _) = GivensRotation.Create(new Complex(1, 2), new Complex(3, 4));
            // rotation2 = ({(0,5701081850582409, -0,019325701188414927)}, 0.82134230050763524)
            (GivensRotation rotation2, _) = GivensRotation.Create(new Complex(4, 5), new Complex(6, 7));
            // rotation3 = ({(0,619902067513713, -0,008669958986205754)}, 0.784631288251623)
            (GivensRotation rotation3, _) = GivensRotation.Create(new Complex(7, 8), new Complex(9, 10));
            var a = new Complex(3, 0);
            var b = new Complex(7, 0);
            var c = new Complex(9, 0);

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

            Assert.That(a.Real, Is.EqualTo(d.Real).Within(1e-12));
            Assert.That(a.Imaginary, Is.EqualTo(d.Imaginary).Within(1e-12));

            Assert.That(b.Real, Is.EqualTo(e.Real).Within(1e-12));
            Assert.That(b.Imaginary, Is.EqualTo(e.Imaginary).Within(1e-12));

            Assert.That(c.Real, Is.EqualTo(f.Real).Within(1e-12));
            Assert.That(c.Imaginary, Is.EqualTo(f.Imaginary).Within(1e-12));
        }

        [TestCase(3, 4, 7, 8)]
        public void Fusion(double ar, double ai, double br, double bi)
        {

            // c1=0.40166320883712187-0.07302967433402213j, s1=-0.9128709291752769
            (var g1, _) = GivensRotation.Create(new Complex(1, 2), new Complex(3, 4));
            // c2=0.5919216793966253-0.01426317299750905j, s2=-0.805869274359261 
            (var g2, _) = GivensRotation.Create(new Complex(5, 6), new Complex(7, 8));

            var a = new Complex(ar, ai);
            var b = new Complex(br, bi);

            // calculate g1 * g2 * (a,b).
            (var c, var d) = g1 * (g2 * (a, b));

            (Complex f, Complex g) = g1.Fusion(g2);
            
            // calculate g1 * diag([f, g]) * (a,b).
            (a, b) = g1 * (f * a, g * b);

            Assert.Multiple(() =>
            {
                Assert.That(c.Real, Is.EqualTo(a.Real).Within(1e-12));
                Assert.That(c.Imaginary, Is.EqualTo(a.Imaginary).Within(1e-12));
                Assert.That(d.Real, Is.EqualTo(b.Real).Within(1e-12));
                Assert.That(d.Imaginary, Is.EqualTo(b.Imaginary).Within(1e-12));
            });
        }

        [Test]
        public void PassthroughDiagonal()
        {
            var d = new Complex(3, 4) / 5;
            var e = new Complex(3, -4) / 5;
            (var rotation, _) = GivensRotation.Create(new Complex(5, 4), new Complex(3, 1));

            (var passedRotation, var f, var g) = AmvwSolver.PassthroughDiagonal(d, e, rotation);

            // checking passedRotation * diag([f, g]) == diag([d, e]) * rotation

            var expected11 = d * rotation.Cosine;
            var expected12 = d * (-rotation.Sine);
            var expected21 = e * rotation.Sine;
            var expected22 = e * rotation.Cosine.Conjugate();

            var actual11 = passedRotation.Cosine * f;
            var actual12 = -passedRotation.Sine * g;
            var actual21 = passedRotation.Sine * f;
            var actual22 = passedRotation.Cosine.Conjugate() * g;

            Assert.Multiple(() =>
            {
                Assert.That(expected11.Real, Is.EqualTo(actual11.Real).Within(1e-12));
                Assert.That(expected11.Imaginary, Is.EqualTo(actual11.Imaginary).Within(1e-12));

                Assert.That(expected12.Real, Is.EqualTo(actual12.Real).Within(1e-12));
                Assert.That(expected12.Imaginary, Is.EqualTo(actual12.Imaginary).Within(1e-12));

                Assert.That(expected21.Real, Is.EqualTo(actual21.Real).Within(1e-12));
                Assert.That(expected21.Imaginary, Is.EqualTo(actual21.Imaginary).Within(1e-12));

                Assert.That(expected22.Real, Is.EqualTo(actual22.Real).Within(1e-12));
                Assert.That(expected22.Imaginary, Is.EqualTo(actual22.Imaginary).Within(1e-12));
            });
        }

        [Test]
        public void Solve()
        {
            // (x - (1+ 2i))(x - (1- 2i))(x - 5) = -25 + 15 x - 7 x^2 + x^3
            var coefficients = new Complex[] { -25, 15, -7 }; // leading coeffient is not provided.

            AmvwSolver.Solve(coefficients);
        }
    }
}
