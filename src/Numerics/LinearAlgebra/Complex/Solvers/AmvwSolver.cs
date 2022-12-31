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
using System.Security.Cryptography;

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
        /// Represents a Givens Rotator, which corresponds to the matrix
        /// [ c, -s;
        ///  s, c'].
        /// </summary>
        public class GivensRotation : ICloneable
        {
            public Complex Cosine { get; set; }
            public double Sine { get; set; }

            public static (Complex a, Complex b) operator *(GivensRotation rotation, (Complex a, Complex b) tuple)
            {
                return (tuple.a * rotation.Cosine - tuple.b * rotation.Sine, tuple.b * rotation.Cosine.Conjugate() + tuple.a * rotation.Sine);
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
                    var originalCosine = (a / alpha).Conjugate();
                    var originalSine = -b / alpha;

                    cosine = (originalCosine.Conjugate() * (b.Conjugate() / bNorm)).Conjugate();
                    sine = -bNorm / alpha;
                }
                return (new GivensRotation(cosine, sine), alpha);
            }

            public static (GivensRotation Rotation, double norm) Create(ref Complex a, ref Complex b)
            {
                (var rotation, var alpha) = Create(a, b);
                a = rotation.Cosine * a - b * rotation.Sine;
                b = 0;
                
                return (rotation, alpha);
            }

            public GivensRotation()
            {
                Cosine = 1;
                Sine = 0;
            }

            public void Invert()
            {
                Cosine = Cosine.Conjugate();
                Sine = -Sine;
            }

            public static (GivensRotation first, GivensRotation second, GivensRotation third) Turnover(GivensRotation g1, GivensRotation g2, GivensRotation g3)
            {
                // calculate a turnover on the second and third row of first column
                
                // Assuming T = g1*g2*g3
                // Doing the turnover operation is the same as calculating the QR factorization
                // of T. First rotation is on the second row and third row of the first column, therefore
                // zeroing the last row in the first column.

                // t21 is the content of first column, second row.
                var t11 = g1.Cosine * g3.Cosine - g1.Sine * g2.Cosine * g3.Sine;
                var t21 = (g1.Cosine.Conjugate() * g2.Cosine * g3.Sine) + g1.Sine * g3.Cosine;

                var t31 = new Complex(g2.Sine * g3.Sine, 0);

                (var r1, var norm1) = GivensRotation.Create(ref t21, ref t31);

                // t21 contains now norm1 because of the call to GivensRotation.Create.

                (var r2, _) = GivensRotation.Create(t11, t21);

                var t12 = -g1.Cosine * g3.Sine - g2.Cosine * g1.Sine * g3.Cosine.Conjugate();
                var t22 = g2.Cosine * (g1.Cosine * g3.Cosine).Conjugate() - g1.Sine * g3.Sine;
                var t32 = g2.Sine * g3.Cosine.Conjugate();

                (t22, t32) = r1 * (t22, t32);
                (_, t22) = r2 * (t12, t22);

                (var r3, _) = GivensRotation.Create(t22, t32);

                // r3*r2*r1 * (g1*g2*g3) = I
                // Therefore (r3*r2*r1)^{-1} = g1*g2*g3
                // r1^{-1} r2^{-1} r3^{-1} = g1*g2*g3

                r3.Invert();
                r2.Invert();
                r1.Invert();

                return (r1, r2, r3);
            }

            /// <summary>
            /// Fuses this rotation with the provided rotation from the right.
            /// Assuming this rotation is g1 and the provided rotation is g2, we
            /// compute g1 = g1*g2. Makes the sine part real again using the 
            /// resulting complex numbers f and g (diagonal matrix entries).
            /// In accordance with page 957 of [1].
            /// </summary>
            /// <param name="rotation">The right rotation.</param>
            /// <returns></returns>
            public (Complex f, Complex g) Fusion(GivensRotation rotation)
            {
                var cosine = Cosine * rotation.Cosine - Sine * rotation.Sine;
                var complexSine = Sine * rotation.Cosine + Cosine.Conjugate() * rotation.Sine;

                Sine = complexSine.Magnitude;
                Complex phi = Complex.One;
                if (Sine != 0)
                {
                    // f is equal to phi.
                    phi = complexSine / Sine;
                    // the formula in [1] is incorrect. 
                    Cosine = cosine / phi;
                }
                else
                {
                    Cosine = cosine;
                }
                var norm = Math.Sqrt(Cosine.MagnitudeSquared() + Sine * Sine);

                Sine /= norm;
                Cosine /= norm;

                return (phi, phi.Conjugate());
            }

            public object Clone()
            {
                return new GivensRotation(Cosine, Sine);
            }
        }

        public class UnitaryMatrix
        {
            public GivensRotation[] Rotations { get; set; }

            public UnitaryMatrix(GivensRotation[] rotations)
            {
                Rotations = rotations;
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
            /// <param name="coefficients">The coefficient vector starting with a_0 without the coefficient a_{n-1}=1.</param>
            /// <returns></returns>
            public static AmvwFactorization Create(Complex[] coefficients)
            {
                var n = coefficients.Length;
                var Q = new GivensRotation[n - 1];
                var D = new Complex[n + 1];

                for (var i = 0; i < n - 1; ++i)
                {
                    Q[i] = new GivensRotation(0, 1);

                    D[i] = Complex.One;
                }
                D[n - 1] = Complex.One;
                D[n] = Complex.One;

                Complex[] x = CreateVectorX(coefficients);

                var B = UnitaryMatrix.Create(x);
                var C = new UnitaryMatrix(B.Rotations.
                    Select(x => new GivensRotation(x.Cosine.Conjugate(), x.Sine)).Reverse().ToArray());

                return new AmvwFactorization() { B = B, C = C, D = D, Q = new(Q), n = n };
            }

            public static Complex[] CreateVectorX(Complex[] coefficients)
            {
                var n = coefficients.Length;
                // create vector x = (-a_1, -a_2, ..., -a_{n-1}, -a_0, 1)
                var x = new Complex[n + 1];
                for (var i = 0; i < n - 1; ++i)
                {
                    x[i] = -coefficients[i + 1];
                }
                x[n - 1] = -coefficients[0];
                x[n] = 1;
                return x;
            }

            public (Complex a, Complex b, Complex c, Complex d) GetPart(int i)
            {
                // Calculate R(i:i+1, i:i+1) = [a, b; c, d]

                // Using equation 4.14 from the paper.
                // H = B + e_1 * y therefore H corresponds to B.
                // h_{{i+1},i} = c_{{i+1},i} * r_{i,i}

                // r_{i,i}
                Complex a = B.Rotations[i].Sine / C.Rotations[n - 1 - i].Sine;
                // r_{i+1,i+1}
                Complex d = B.Rotations[i + 1].Sine / C.Rotations[n - 1 - (i + 1)].Sine;

                // h_{i,i} = c_{i,i-1} * r_{i-1,i} + c_{i,i} * r_{i,i}
                // r_{i-1,i} = (h_{i,i} - c_{i,i} * r_{i,i}) / c_{i,i-1}
                // r_{i,i+1} = (h_{i+1, i+1} - c_{i+1, i+1} * r_{i+1, i+1}) / c_{i+1, i}
                Complex b = (B.Rotations[i].Cosine - C.Rotations[n - 1 - i].Cosine.Conjugate() * d) / C.Rotations[n - 1 - i].Sine;

                Complex c = 0;

                // Apply diagonal:
                a *= D[i];
                b *= D[i];
                d *= D[i + 1];

                // Apply Q
                if (i + 1 < n - 1)
                {
                    d = Q.Rotations[i + 1].Cosine * d;
                }

                (a, c) = Q.Rotations[i] * (a, c);
                (b, d) = Q.Rotations[i] * (b, d);

                return (a, b, c, d);
            }

            public (Complex e1, Complex e2) GetWilkinsonShift()
            {
                (var AUpperLeft, var AUpperRight, var ALowerLeft, var ALowerRight) = GetPart(n - 1);
                return AmvwSolver.GetWilkinsonShift(AUpperLeft, AUpperRight, ALowerLeft, ALowerRight);
            }
        }

        internal static (Complex a, Complex b) GetWilkinsonShift(Complex AUpperLeft, Complex AUpperRight, Complex ALowerLeft, Complex ALowerRight)
        {
            var trace = AUpperLeft + ALowerRight;
            var determinant = AUpperLeft * ALowerRight - ALowerLeft * AUpperRight;

            var disc = (trace * trace - 4 * determinant).SquareRoot();

            var traceDiscSum = trace + disc;
            var traceDiscDiff = trace - disc;

            var e1 = Complex.Zero;
            var e2 = Complex.Zero;

            if (traceDiscSum.Magnitude > traceDiscDiff.Magnitude)
            {
                if (traceDiscSum.Magnitude > 0)
                {
                    e1 = traceDiscSum / 2;
                    e2 = determinant / e1;
                }
            }
            else
            {
                if (traceDiscDiff.Magnitude > 0)
                {
                    e1 = traceDiscDiff / 2;
                    e2 = determinant / e1;
                }
            }

            return (e1, e2);
        }

        public static (GivensRotation, Complex f, Complex g) PassthroughDiagonal(Complex d, Complex e, GivensRotation rotation)
        {
            var f = e;
            var g = d;
            var s2 = rotation.Sine;
            var c2 = rotation.Cosine * d * e.Conjugate();
            return (new GivensRotation(c2, s2), f, g);
        }

        public List<Complex> Solve(Polynomial polynomial)
        {
            return Solve(polynomial.Coefficients);
        }

        public List<Complex> Solve(double[] coefficients)
        {
            return Solve(coefficients.Select(x => new Complex(x, 0)).ToArray());
        }

        public static List<Complex> Solve(Complex[] coefficients)
        {
            var m = coefficients.Length;

            // The companion matrix A is factorized into:
            // A = Q * C' * (B' + e1 * y')
            // where Q, C, B are unitary and stored using Givens rotations.
            var factorization = AmvwFactorization.Create(coefficients);
            var n = factorization.n;
            Complex Aul; Complex Aur; Complex All; Complex Alr;

            for (int j = 0; j < 30; ++j)
            {
                (var A11, var A12, var A21, var A22) = factorization.GetPart(0);
                (Aul, Aur, All, Alr) = factorization.GetPart(n-2);
                
                (var shift1, var shift2) = GetWilkinsonShift(Aul, Aur, All, Alr);
                Console.WriteLine("Aul={0}, Aur={1}, All={2}, Alr={3}", Aul, Aur, All, Alr);

                //var shift = (shift1 - Alr).MagnitudeSquared() < (shift2 - Alr).MagnitudeSquared() ? shift1 : shift2;
                var shift = 0;

                (var u1, _) = GivensRotation.Create(A11 - shift, A21);

                var u1H = (GivensRotation)u1.Clone();
                u1H.Invert();

                // fuse u1 with the first Q core transformation.
                (factorization.D[0], factorization.D[1]) = u1H.Fusion(factorization.Q.Rotations[0]);

                factorization.Q.Rotations[0] = u1H;

                var transformation = u1;
                // i is equal to the active rows (i,i+1) the current rotation works on.
                for (int i = 0; i < n;)
                {
                    // Turnover on B:
                    (transformation, factorization.B.Rotations[i], factorization.B.Rotations[i + 1]) =
                        GivensRotation.Turnover(factorization.B.Rotations[i], factorization.B.Rotations[i + 1], transformation);
                    // the transformation moved down one row:
                    i += 1;

                    // Turnover on C:
                    (transformation, factorization.C.Rotations[n - 1 - i], factorization.C.Rotations[n - 1 - (i - 1)]) =
                        GivensRotation.Turnover(factorization.C.Rotations[n - 1 - i], factorization.C.Rotations[n - 1 - (i - 1)], transformation);

                    i -= 1;

                    // pass through diagonal
                    (transformation, factorization.D[i], factorization.D[i + 1]) =
                        PassthroughDiagonal(factorization.D[i], factorization.D[i + 1], transformation);

                    if (i == n - 2)
                    {
                        // Fuse the transformation with the last rotation in Q.
                        factorization.Q.Rotations[i].Fusion(transformation);
                        break;
                    }
                    else
                    {
                        // Turnover on Q:
                        (transformation, factorization.Q.Rotations[i], factorization.Q.Rotations[i + 1]) =
                            GivensRotation.Turnover(factorization.Q.Rotations[i], factorization.Q.Rotations[i + 1], transformation);
                        i += 1;
                    }
                }
            }

            (Aul, Aur, All, Alr) = factorization.GetPart(n - 2);
            Console.WriteLine("Aul={0}, Alr={1}", Aul, Alr);

            return null;
        }
    }
}
