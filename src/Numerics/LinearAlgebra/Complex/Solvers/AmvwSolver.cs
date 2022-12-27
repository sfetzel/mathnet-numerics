// <copyright file="AmvwSolver.cs" company="Math.NET">
// Math.NET Numerics, part of the Math.NET Project
// http://numerics.mathdotnet.com
// http://github.com/mathnet/mathnet-numerics
//
// Copyright (c) 2009-2014 Math.NET
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
using MathNet.Numerics.LinearAlgebra.Solvers;
using System.Collections.Generic;
using System.Linq;
using Complex = System.Numerics.Complex;
using System.Reflection.Metadata.Ecma335;

namespace MathNet.Numerics.LinearAlgebra.Complex.Solvers
{
    using Complex = System.Numerics.Complex;
    /// <summary>
    /// Calculates the roots of a complex polynomial p(z) = z^n + a_{n-1}z^{n-1} + ... + a_1 z + a_0
    /// using the algorithm described by [1]. The coefficient vector is therefore of length n.
    /// </summary>
    public class AmvwSolver
    {
        /// <summary>
        /// Represents a Givens Rotator, which corresponds to the matrix [c, s; -s, c'].
        /// </summary>
        public class GivensRotation
        {
            public Complex Cosine { get; set; }
            public double Sine { get; set; }

            public static (Complex a, Complex b) operator *(GivensRotation rotation, (Complex a, Complex b) tuple)
            {
                return (tuple.a * rotation.Cosine + tuple.b * rotation.Sine, tuple.b * rotation.Cosine.Conjugate() - tuple.a * rotation.Sine);
            }

            public GivensRotation(Complex cosine, double sine)
            {
                Cosine = cosine;
                Sine = sine;
            }

            public static (GivensRotation Rotation, double norm) Create(Complex a, Complex b)
            {
                var bNorm = b.Magnitude;
                var cosine = Complex.One;
                var sine = 0.0;
                var alpha = Math.Sqrt(a.MagnitudeSquared() + bNorm * bNorm);

                if (bNorm != 0 && alpha != 0)
                {
                    // the rotation must be described with a complex cosine
                    // and a real sine.
                    // This is done by factoring out the phase of the sine part.

                    cosine = (a * (b / bNorm).Conjugate() / alpha).Conjugate();
                    sine = bNorm / alpha;
                }
                return (new GivensRotation(cosine, sine), alpha);
            }

            public static (GivensRotation Rotation, double norm) Create(ref Complex a, ref Complex b)
            {
                (var rotation, var alpha) = Create(a, b);
                a = rotation.Cosine * a + b * rotation.Sine;
                b = 0;
                
                return (rotation, alpha);
            }

            public GivensRotation()
            {
                Cosine = 1;
                Sine = 0;
            }

            public static (GivensRotation first, GivensRotation second, GivensRotation third) Turnover(GivensRotation leftRotation, GivensRotation middleRotation, GivensRotation rightRotation)
            {
                // calculate a turnover on the second and third row of first column

                // t21 is the content of first column, second row.
                var t21 = leftRotation.Sine * rightRotation.Cosine + leftRotation.Cosine.Conjugate() * middleRotation.Cosine * rightRotation.Sine;
                var t31 = new Complex(middleRotation.Sine * rightRotation.Sine, 0);

                (var firstRotation, var norm1) = GivensRotation.Create(ref t21, ref t31);

                var t11 = leftRotation.Cosine * rightRotation.Cosine - leftRotation.Sine * middleRotation.Cosine * rightRotation.Sine;
                // t21 contains now norm1 because of the call to GivensRotation.Create.

                (var secondRotation, _) = GivensRotation.Create(ref t11, ref t21);

                var t23 = middleRotation.Cosine * firstRotation.Cosine + leftRotation.Cosine * middleRotation.Sine * firstRotation.Sine;
                // ensuring zero imaginary part by explicitly calculating the real part.
                var t33 = new Complex(leftRotation.Sine * middleRotation.Sine * secondRotation.Sine
                    - firstRotation.Sine * (secondRotation.Cosine.Real * middleRotation.Cosine.Real + secondRotation.Cosine.Imaginary * middleRotation.Cosine.Imaginary)
                    + middleRotation.Sine * (
                        leftRotation.Cosine.Real * (firstRotation.Cosine.Real * secondRotation.Cosine.Real + firstRotation.Cosine.Imaginary * secondRotation.Cosine.Imaginary)
                        + firstRotation.Cosine.Imaginary * (firstRotation.Cosine.Real * secondRotation.Cosine.Imaginary - firstRotation.Cosine.Imaginary * secondRotation.Cosine.Real)
                    ), 0);

                (var thirdRotation, _) = GivensRotation.Create(ref t23, ref t33);

                return (firstRotation, secondRotation, thirdRotation);
            }
        }

        public class UnitaryMatrix
        {
            public GivensRotation[] CoreTransformations { get; set; }

            public UnitaryMatrix(GivensRotation[] rotations)
            {
                CoreTransformations = rotations;
            }

            /// <summary>
            /// Creates a unitary matrix Q of givens rotations,
            /// which transforms the vector v to a multiple of its norm multiplied with e_1.
            /// </summary>
            /// <param name="vector"></param>
            /// <returns></returns>
            public static UnitaryMatrix Create(Complex[] vector)
            {
                var n = vector.Length;
                GivensRotation[] rotations = new GivensRotation[n-1];
                for(int i = n-1; i > 0; --i)
                {
                    (rotations[i - 1], _) = GivensRotation.Create(ref vector[i - 1], ref vector[i]);
                }
                return new UnitaryMatrix(rotations);
            }
        }

        /// <summary>
        /// Represents the factorization A = Q * C^H * (B + e1*y^T)
        /// </summary>
        public class AmvwFactorization
        {
            /// <summary>
            /// Represents the Q part of the factorization. The i-th transformation
            /// is active on rows i and i+1.
            /// </summary>
            public UnitaryMatrix Q { get; set; }

            public Complex[] D { get; set; }

            /// <summary>
            /// Represents the hermitian C matrix of the factorization. The i-th transformation
            /// is active on rows (n-(i+1), n-i).
            /// </summary>
            public UnitaryMatrix C { get; set; }

            /// <summary>
            /// Represents the B part of the factorization. The i-th transformation
            /// is active on rows (i, i+1).
            /// </summary>
            public UnitaryMatrix B { get; set; }

            public int n { get; set; }

            private AmvwFactorization() { }

            /// <summary>
            /// 
            /// </summary>
            /// <param name="coefficients">The coefficient vector starting with a_0.</param>
            /// <returns></returns>
            public static AmvwFactorization Create(Complex[] coefficients)
            {
                var m = coefficients.Length;
                var Q = new GivensRotation[m];
                var D = new Complex[m];

                for (var i = 0; i < m - 1; ++i)
                {
                    Q[i] = new GivensRotation(0, -1);

                    D[i] = Complex.One;
                }
                Q[m - 1] = new GivensRotation(1, 0);
                Complex[] x = CreateVectorX(coefficients);

                var B = UnitaryMatrix.Create(x);
                var C = new UnitaryMatrix(B.CoreTransformations.
                    Select(x => new GivensRotation(x.Cosine.Conjugate(), x.Sine)).Reverse().ToArray());

                return new AmvwFactorization() { B = B, C = C, D = D, Q = new(Q), n = m };
            }

            public static Complex[] CreateVectorX(Complex[] coefficients)
            {
                var m = coefficients.Length;
                // create vector x = (-a_1, -a_2, ..., -a_{n-1}, -a_0, 1)
                var x = new Complex[m + 1];
                for (var i = 0; i < m - 1; ++i)
                {
                    x[i] = -coefficients[i + 1];
                }
                x[m - 1] = -coefficients[0];
                x[m] = 1;
                return x;
            }

            public (Complex a, Complex b, Complex c, Complex d) GetPart(int i)
            {
                // Calculate R(i:i+1, i:i+1) = [a, b; c, d]

                // Using equation 4.14 from the paper.
                // H = B + e_1 * y therefore H corresponds to B.
                // h_{{i+1},i} = c_{{i+1},i} * r_{i,i}

                // r_{i,i}
                Complex a = B.CoreTransformations[i].Sine / C.CoreTransformations[n - 1 - i].Sine;
                // r_{i+1,i+1}
                Complex d = B.CoreTransformations[i + 1].Sine / C.CoreTransformations[n - 1 - (i + 1)].Sine;

                // h_{i,i} = c_{i,i-1} * r_{i-1,i} + c_{i,i} * r_{i,i}
                // r_{i-1,i} = (h_{i,i} - c_{i,i} * r_{i,i}) / c_{i,i-1}
                // r_{i,i+1} = (h_{i+1, i+1} - c_{i+1, i+1} * r_{i+1, i+1}) / c_{i+1, i}
                Complex b = (B.CoreTransformations[i].Cosine - C.CoreTransformations[n - 1 - i].Cosine.Conjugate() * d) / C.CoreTransformations[n - 1 - i].Sine;

                Complex c = 0;

                d = Q.CoreTransformations[i + 1].Cosine * d;

                (a, c) = Q.CoreTransformations[i] * (a, c);
                (b, d) = Q.CoreTransformations[i] * (b, d);

                return (a, b, c, d);
            }
        }

        public List<Complex> Solve(Polynomial polynomial)
        {
            return Solve(polynomial.Coefficients);
        }

        public List<Complex> Solve(double[] coefficients)
        {
            return Solve(coefficients.Select(x => new Complex(x, 0)).ToArray());
        }

        public List<Complex> Solve(Complex[] coefficients)
        {
            var m = coefficients.Length;

            // The companion matrix A is factorized into:
            // A = Q * C' * (B' + e1 * y')
            // where Q, C, B are unitary and stored using Givens rotations.
            var factorization = AmvwFactorization.Create(coefficients);


            return null;
        }
    }
}
