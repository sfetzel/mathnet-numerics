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
using System.Globalization;

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

                if (alpha != 0)
                {
                    var originalCosine = (a / alpha).Conjugate();
                    if (bNorm == 0)
                    {
                        // If b is zero, then a*a' will always be a positive real number.
                        cosine = originalCosine;
                    }
                    else
                    {
                        // the rotation must be described with a complex cosine
                        // and a real sine.
                        // This is done by factoring out the phase of the sine part.
                        var originalSine = -b / alpha;
                        cosine = (originalCosine.Conjugate() * (b.Conjugate() / bNorm)).Conjugate();
                        sine = -bNorm / alpha;
                    }
                }
                return (new GivensRotation(cosine, sine), alpha);
            }

            public static (GivensRotation Rotation, double norm) CreatePositive(Complex a, Complex b)
            {
                (var rotation, var alpha) = Create(a, b);
                if (b.Real < 0)
                {
                    rotation.Cosine *= -1;
                    rotation.Sine *= -1;
                }
                return (rotation, alpha);
            }


            public static (GivensRotation Rotation, double norm) CreatePositive(ref Complex a, ref Complex b)
            {
                (var rotation, var alpha) = CreatePositive(a, b);
                a = rotation.Cosine * a - b * rotation.Sine;
                b = 0;

                return (rotation, alpha);
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

            public double Norm() => Math.Sqrt(Cosine.MagnitudeSquared() + Sine * Sine);

            public void Normalize()
            {
                var norm = Norm();
                Cosine /= norm;
                Sine /= norm;
            }

            public static GivensRotation PassthroughRotation(GivensRotation rotation, Complex cosineCorrectionTerm)
            {
                var s2 = rotation.Sine;
                var c2 = rotation.Cosine * cosineCorrectionTerm;
                var newRotation = new GivensRotation(c2, s2);
                //newRotation.Normalize();
                return newRotation;
            }

            /// <summary>
            /// Passes a rotation R through a diagonal D=diag([d,e]) from "left to right", such that:
            /// diag([d, e]) * R == G * diag([d,e])
            /// where G is the new rotation and d,e are the new diagonal elements.
            /// </summary>
            /// <param name="d"></param>
            /// <param name="e"></param>
            /// <param name="rotation"></param>
            /// <returns></returns>
            public static (GivensRotation, Complex f, Complex g) PassthroughRotationToLeft(Complex d, Complex e, GivensRotation rotation)
                => (PassthroughRotation(rotation, d * e.Conjugate()), e, d);

            public static (GivensRotation, Complex f, Complex g) PassthroughRotationToLeftHermitian(Complex d, Complex e, GivensRotation rotation)
                => (PassthroughRotation(rotation, d.Conjugate() * e), e, d);


            /// <summary>
            /// Passes a rotation R through a diagonal D=diag([d,e]), such that:
            /// R * diag([d, e]) == diag([f,g]) * G,
            /// where G is the new rotation and d,e are the new diagonal elements.
            /// </summary>
            /// <param name="d"></param>
            /// <param name="e"></param>
            /// <param name="rotation"></param>
            /// <returns></returns>
            public static (Complex f, Complex g, GivensRotation) PassthroughRotationToRight(GivensRotation rotation, Complex d, Complex e)
            {
                var f = e;
                var g = d;
                var s2 = rotation.Sine;
                var c2 = rotation.Cosine / (d * e.Conjugate());
                return (f, g, new GivensRotation(c2, s2));
            }

            public static (GivensRotation first, GivensRotation second, GivensRotation third) Turnover(GivensRotation g1, GivensRotation g2, GivensRotation g3)
                => Turnover(g1.Cosine, g1.Sine, g2.Cosine, g2.Sine, g3.Cosine, g3.Sine);

            public static (GivensRotation r1, GivensRotation r2, GivensRotation r3) TurnoverComplexConjugate(GivensRotation g1, GivensRotation g2, GivensRotation g3)
                => TurnoverComplexConjugate(g1.Cosine.Conjugate(), -g1.Sine, g2.Cosine.Conjugate(), -g2.Sine, g3.Cosine, g3.Sine);

            public static (GivensRotation first, GivensRotation second, GivensRotation third) Turnover(Complex g1Cosine, double g1Sine, Complex g2Cosine, double g2Sine, Complex g3Cosine, double g3Sine)
            {
                (var r1, var r2, var r3) = TurnoverComplexConjugate(g1Cosine, g1Sine, g2Cosine, g2Sine, g3Cosine, g3Sine);

                r3.Invert();
                r2.Invert();

                return (r1, r2, r3);
            }

            public static (GivensRotation r1, GivensRotation r2, GivensRotation r3) TurnoverComplexConjugate(Complex g1Cosine, double g1Sine, Complex g2Cosine, double g2Sine, Complex g3Cosine, double g3Sine)
            {
                // calculate a turnover on the second and third row of first column

                // Assuming T = g1*g2*g3
                // Doing the turnover operation is the same as calculating the QR factorization
                // of T. First rotation is on the second row and third row of the first column, therefore
                // zeroing the last row in the first column.

                // t21 is the content of first column, second row.
                var t11 = g1Cosine * g3Cosine - g1Sine * g2Cosine * g3Sine;
                var t21 = (g1Cosine.Conjugate() * g2Cosine * g3Sine) + g1Sine * g3Cosine;

                var t31 = new Complex(g2Sine * g3Sine, 0);

                (var r1, var norm1) = GivensRotation.CreatePositive(ref t21, ref t31);

                if (t21.Real < -1e-15)
                {
                    throw new Exception();
                }

                // t21 contains now norm1 because of the call to GivensRotation.Create.

                (var r2, _) = GivensRotation.CreatePositive(ref t11, ref t21);
                if (t11.Real < -1e-15)
                {
                    throw new Exception();
                }

                var t12 = -g1Cosine * g3Sine - g2Cosine * g1Sine * g3Cosine.Conjugate();
                var t22 = g2Cosine * (g1Cosine * g3Cosine).Conjugate() - g1Sine * g3Sine;
                var t32 = g2Sine * g3Cosine.Conjugate();

                (t22, t32) = r1 * (t22, t32);
                (_, t22) = r2 * (t12, t22);

                (var r3, _) = GivensRotation.CreatePositive(ref t22, ref t32);

                if (t22.Real < -1e-15)
                {
                    throw new Exception();
                }

                // r3*r2*r1 * (g1*g2*g3) = I
                // Therefore (r3*r2*r1)^{-1} = g1*g2*g3
                // r1^{-1} r2^{-1} r3^{-1} = g1*g2*g3

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

                (var phi, Cosine, Sine) = GetRotation(cosine, complexSine);

                return (phi, phi.Conjugate());
            }

            public static (Complex phi, Complex Cosine, double sine) GetRotation(Complex cosine, Complex complexSine)
            {
                double Sine = complexSine.Magnitude;
                Complex Cosine;
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
                return (phi, Cosine, Sine);
            }

            public Matrix<Complex> GetMatrix(int n, int i)
            {
                var matrix = Matrix.Build.DenseIdentity(n, n);
                matrix[i, i] = Cosine;
                matrix[i, i + 1] = -Sine;
                matrix[i + 1, i] = Sine;
                matrix[i + 1, i + 1] = Cosine.Conjugate();
                return matrix;
            }

            public object Clone()
            {
                return new GivensRotation(Cosine, Sine);
            }

            public GivensRotation CloneInverse()
            {
                return new GivensRotation(Cosine.Conjugate(), -Sine);
            }

            public static (Complex d, Complex e, GivensRotation) GetReflectedRotation(GivensRotation rotation)
            {
                // This swaps cosine and sine part.
                // We must make the sine real by using the diagonal matrix D.
                // we have the matrix:
                // [ -s,  c;
                //    c', s ]
                // factor for row one: c'/(cosineNorm) = phase.
                // factor for row two: c/cosineNorm = -phase'.
                var cosineNorm = rotation.Cosine.Magnitude;
                var sine = rotation.Sine;
                var phase = rotation.Cosine.Conjugate() / cosineNorm;

                var newRotation = new GivensRotation(-sine * phase, -(rotation.Cosine * phase).Real);

                // The elements d,e of the diagonal matrix are the inverse of the factors.
                return (1 / phase, -1 / phase.Conjugate(), newRotation);
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
                GivensRotation[] rotations = new GivensRotation[n - 1];
                for (int i = n - 1; i > 0; --i)
                {
                    (rotations[i - 1], _) = GivensRotation.CreatePositive(ref vector[i - 1], ref vector[i]);
                }
                return new UnitaryMatrix(rotations);
            }

            public Matrix<Complex> GetMatrix(int n)
                => Rotations.Select((x, index) => x.GetMatrix(n, index)).Aggregate((a, b) => a * b);

            /// <summary>
            /// Gets a diagonal element of this matrix.
            /// </summary>
            /// <param name="i">The index of the diagonal matrix, starting from 0.</param>
            /// <returns></returns>
            public Complex GetDiagonal(int i) =>
                (i < Rotations.Length ? Rotations[i].Cosine : 1) * (i > 0 ? Rotations[i - 1].Cosine.Conjugate() : 1);

            public Complex GetSuperdiagonal(int i) =>
                (i > Rotations.Length ? 0 : 1) *
                -(i < Rotations.Length ? Rotations[i].Cosine : 1) *
                 (i - 1 >= 0 && i - 1 < Rotations.Length ? Rotations[i - 1].Sine : 1) *
                (i - 2 >= 0 && i - 2 < Rotations.Length ? Rotations[i - 2].Cosine.Conjugate() : 1);

            /// <summary>
            /// Applys a transformation T from the right hand side to this unitary matrix U and
            /// returns the transformation S such that S * U = U * T, where
            /// S moved one row up.
            /// </summary>
            /// <param name="transformation">The transformation T.</param>
            /// <param name="i">The first row on which T is active.</param>
            /// <returns>The transformation S and the first row on which S is active.</returns>
            public (GivensRotation, int) ApplyTransformationRight(GivensRotation transformation, int i)
            {
                if (i < Rotations.Length - 1)
                {
                    // Using a turnover operation to pass the transformation through the matrix.
                    // Input as pictogram is: Up down up.
                    (transformation, Rotations[i], Rotations[i + 1]) =
                        GivensRotation.Turnover(Rotations[i], Rotations[i + 1], transformation);
                    // result is Down Up Down.
                    i += 1;
                }
                else
                {
                    throw new ArgumentException($"Cannot apply a transformation active on rows {i}:{i + 1}.");
                }
                return (transformation, i);
            }

            /// <summary>
            /// Applys a transformation T from the left hand side to this unitary matrix U and
            /// returns the resulting diagonal matrix D on the right hand side:
            /// T*U = U' * D
            /// </summary>
            /// <param name="transformation"></param>
            /// <param name="i"></param>
            /// <returns>the diagonal matrix D.</returns>
            public Complex[] ApplyTransformationLeft(GivensRotation transformation, int i)
            {
                if (i > 0)
                {
                    throw new ArgumentException("currently not supported");
                }

                var D = Enumerable.Repeat(Complex.One, Rotations.Length + 1).ToArray();
                var newTransformation = (GivensRotation)transformation.Clone();
                (D[i], D[i + 1]) = newTransformation.Fusion(Rotations[i]);
                Rotations[i] = newTransformation;
                PassthroughDiagonalToRight(D, i + 1);
                return D;
            }

            /// <summary>
            /// Applys a transformation T from the right hand side to complex conjugate of this unitary matrix U and
            /// returns the transformation S such that S * U' = U' * T, where
            /// U moved one row up.
            /// </summary>
            /// <param name="transformation">The transformation T.</param>
            /// <param name="i">The first row on which T is active.</param>
            /// <returns>The transformation S and the first row on which S is active.</returns>
            public (GivensRotation, int) ApplyTransformationConjugate(GivensRotation transformation, int i)
            {
                if (i < Rotations.Length && i >= 1)
                {
                    // Using a turnover operation to pass the transformation through the matrix.
                    // Input as pictogram is: Down Up Down.
                    (transformation, Rotations[i], Rotations[i - 1]) =
                        GivensRotation.TurnoverComplexConjugate(Rotations[i], Rotations[i - 1], transformation);
                    // result is Up Down Up.
                    i -= 1;
                }
                else
                {
                    throw new ArgumentException($"Cannot apply a transformation active on rows {i}:{i + 1}.");
                }
                return (transformation, i);
            }

            public Complex[] PassthroughDiagonalToLeft(Complex[] diagonal) => PassthroughDiagonalToLeft(diagonal, diagonal.Length - 2);

            public Complex[] PassthroughDiagonalToLeft(Complex[] diagonal, int start)
            {
                if (diagonal.Length != Rotations.Length + 1)
                {
                    throw new ArgumentException(nameof(diagonal));
                }
                var n = diagonal.Length - 1;
                for (int i = start; i >= 0; --i)
                {
                    // Rotation is active on rows i and i+1.
                    (Rotations[i], diagonal[i], diagonal[i + 1]) = GivensRotation.PassthroughRotationToLeft(diagonal[i], diagonal[i + 1], Rotations[i]);
                }
                return diagonal;
            }


            public Complex[] PassthroughDiagonalToRight(Complex[] diagonal, int start)
            {
                if (diagonal.Length != Rotations.Length + 1)
                {
                    throw new ArgumentException(nameof(diagonal));
                }
                for (int i = start; i < Rotations.Length; ++i)
                {
                    // Rotation is active on rows i and i+1.
                    (Rotations[i], diagonal[i], diagonal[i + 1]) = GivensRotation.PassthroughRotationToLeft(diagonal[i], diagonal[i + 1], Rotations[i]);
                    diagonal[i] /= diagonal[i].Magnitude;
                    diagonal[i + 1] /= diagonal[i + 1].Magnitude;
                }
                return diagonal;
            }

            public Complex[] PassthroughDiagonalHermitian(Complex[] diagonal)
            {
                if (diagonal.Length != Rotations.Length + 1)
                {
                    throw new ArgumentException(nameof(diagonal));
                }
                var n = diagonal.Length - 1;
                for (var i = 0; i < Rotations.Length; ++i)
                {
                    // Rotation is active on rows i and i+1.
                    (Rotations[i], diagonal[i], diagonal[i + 1]) = GivensRotation.PassthroughRotationToLeftHermitian(diagonal[i], diagonal[i + 1], Rotations[i]);
                }
                return diagonal;
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
            /// Represents the C matrix of the factorization. The i-th transformation
            /// is active on rows (n-(i+1), n-i).
            /// </summary>
            public UnitaryMatrix C { get; set; }

            /// <summary>
            /// Represents the B part of the factorization. The i-th transformation
            /// is active on rows (i, i+1).
            /// </summary>
            public UnitaryMatrix B { get; set; }

            public int n { get; set; }

            public Complex y { get; set; }

            public AmvwFactorization() { }

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

                // correct negative sign which is introduced from Q_1,..,Q_{n-1}.
                Complex[] x = CreateVectorX(coefficients);
                // B_i = C_i for i=1,.., n-1
                var C = UnitaryMatrix.Create(x);
                Complex yNorm = x[0].Real;

                var B = new UnitaryMatrix(C.Rotations.
                    Select(x => new GivensRotation(x.Cosine, x.Sine)).ToArray());

                // B*Z_n is almost a diagonal matrix, but the last two rows/columns are [0, 1; 1, 0].
                // This swaps cosine and sine part.
                // We must make the sine real by using the diagonal matrix D.
                // factor for row one: c/(cosineNorm) = phase.
                // factor for row two: -c'/cosineNorm = -phase'.

                (D[n - 1], D[n], B.Rotations[n - 1]) = GivensRotation.GetReflectedRotation(B.Rotations[n - 1]);

                // correct negative sign which is introduced from Q_1,..,Q_{n-1}:
                // Z_n becomes [0, -1; 1, 0];
                //D[n-1] *= n % 2 == 0 ? -1 : 1;

                D = B.PassthroughDiagonalToLeft(D, n - 2);
                yNorm /= D[0];
                D = C.PassthroughDiagonalHermitian(D);

                // correct negative sign which is introduced from Q_1,..,Q_{n-1}.
                // this happens when n is an even number.
                if (n % 2 == 0)
                {
                    D[n - 1] *= -1;
                }


                return new AmvwFactorization() { B = B, C = C, D = D, Q = new(Q), n = n, y = yNorm };
            }

            public Matrix<Complex> GetDiagonalMatrix() => Matrix<Complex>.Build.DiagonalOfDiagonalArray(D);

            public Matrix<Complex> GetYMatrix()
            {
                var yMatrix = Matrix<Complex>.Build.Sparse(n + 1, n + 1, 0);
                yMatrix[0, n - 1] = y;
                return yMatrix;
            }

            public Vector<Complex> GetDefaultY()
            {
                var yVector = Vector<Complex>.Build.Sparse(n + 1);
                yVector[0] = y;
                return yVector;
            }

            public Matrix<Complex> GetMatrix() => GetMatrix(GetDefaultY());

            public Matrix<Complex> GetMatrix(Vector<Complex> y)
            {
                var Cstar = C.GetMatrix(n + 1).ConjugateTranspose();
                var H = B.GetMatrix(n + 1);
                H.SetColumn(n - 1, H.Column(n - 1) + y);
                var R = GetDiagonalMatrix() * Cstar * H;
                var A = Q.GetMatrix(n + 1) * R;
                return A;
            }

            public static Complex[] CreateVectorX(Complex[] coefficients)
            {
                var n = coefficients.Length;
                // create vector x = (-a_1, -a_2, ..., -a_{n-1}, -a_0, -1)
                var x = new Complex[n + 1];
                for (var i = 0; i < n - 1; ++i)
                {
                    x[i] = -coefficients[i + 1];
                }
                // compensate the sign from the Q core transformations.
                x[n - 1] = -coefficients[0];
                x[n] = -1;
                return x;
            }

            public Complex GetBestWilkinsonShift(int n)
            {
                (var Aul, var Aur, var All, var Alr) = GetPart(n - 2);
                (var shift1, var shift2) = AmvwSolver.GetWilkinsonShift(Aul, Aur, All, Alr);
                // select the shift which is nearer to the bottom right entry of A.
                var shift = (shift1 - Alr).MagnitudeSquared() < (shift2 - Alr).MagnitudeSquared() ? shift1 : shift2;
                return shift;
            }
            public Complex GetBestRayleighShift(int n)
            {
                (_, _, _, var Alr) = GetPart(n - 2);
                return Alr;
            }

            public GivensRotation GetFirstTransformation(int n, ShiftStrategy strategy)
                => strategy switch
                {
                    ShiftStrategy.Wilkinson => GetFirstTransformation(n, GetBestWilkinsonShift(n)),
                    ShiftStrategy.Rayleigh => GetFirstTransformation(n, GetBestRayleighShift(n)),
                    ShiftStrategy.Random => GetFirstTransformation(n, new Complex(RandomNumberGenerator.GetInt32(n), 0)),
                    _ => GetFirstTransformation(n, Complex.Zero),
                };

            public GivensRotation GetFirstTransformation(int n) => GetFirstTransformation(n, GetBestWilkinsonShift(n));

            public GivensRotation GetFirstTransformation(int n, Complex shift)
            {
                (var A11, _, var A21, _) = GetTopLeft();

                (var u1, _) = GivensRotation.Create(A11 - shift, A21);
                return u1;
            }

            /// <summary>
            /// Gets the entry R_{i,i} of the R matrix (A=Q*R).
            /// Using equation 4.12 of [1]: h_{i+1, i}/c_{i+1, i}.
            /// </summary>
            /// <param name="i">diagonal entry (starting from 0).</param>
            /// <returns></returns>
            public Complex GetDiagonalOfR(int i) => B.Rotations[i].Sine / C.Rotations[i].Sine;

            public Complex GetDiagonalOfH(int i) => B.Rotations[i].Cosine * (i > 0 ? B.Rotations[i - 1].Cosine.Conjugate() : 1);

            /// <summary>
            /// Gets the hermitian Givens rotation of the C matrix, which is active on rows i and i+1.
            /// </summary>
            /// <param name="i"></param>
            /// <returns></returns>
            public GivensRotation GetCRotation(int i) => C.Rotations[n - 1 - i];

            public Complex GetDiagonalOfC(int i) => GetCRotation(i).Cosine.Conjugate() * (i == n ? 1.0 : GetCRotation(i - 1).Cosine);

            public (Complex a, Complex b, Complex c, Complex d) GetTopLeftR()
            {
                // Calculate R(0:1, 0:1) = [a, b; c, d]

                // Using equation 4.14 from the paper.
                // H = B + e_1 * y therefore H corresponds to B.
                // h_{{i+1},i} = c_{{i+1},i} * r_{i,i}

                // r_{i,i}
                Complex a = D[0] * GetDiagonalOfR(0);
                // r_{i+1,i+1}
                Complex d = GetDiagonalOfR(1);

                // h_{i,i} = c_{i,i-1} * r_{i-1,i} + c_{i,i} * r_{i,i}
                // r_{i-1,i} = (h_{i,i} - c_{i,i} * r_{i,i}) / c_{i,i-1}

                // r_{i,i+1} = (h_{i+1, i+1} - c_{i+1, i+1} * r_{i+1, i+1}) / c_{i+1, i}
                Complex b = D[0] * (GetDiagonalOfH(1) - C.GetDiagonal(1) * d) / C.Rotations[0].Sine;

                Complex c = 0;
                d *= D[1];
                return (a, b, c, d);
            }

            public (Complex a, Complex b, Complex c, Complex d) GetTopLeft()
            {
                (var a, var b, var c, var d) = GetTopLeftR();

                // Apply Q
                c = Q.Rotations[1].Cosine * c;
                d = Q.Rotations[1].Cosine * d;

                (a, c) = Q.Rotations[0] * (a, c);
                (b, d) = Q.Rotations[0] * (b, d);

                return (a, b, c, d);
            }

            /// <summary>
            /// Returns the matrix (octave/matlab indexing but i from parameter)
            /// R(i-1:i+1, i:i+1) = [ R00, R01;
            ///                     R10, R11;
            ///                     R20, R21; ].
            /// The diagonal elements are R10, R21. The diagonal is already applied to this matrix.
            /// R00 refers to R(i-1, i)
            /// R10 refers to R(i, i)
            /// R(i+1, i) is zero and will not be calculated.
            /// R01 refers to R(i-1, i+1)
            /// R11 refers to R(i, i+1)
            /// R21 refers to R(i+1, i+1)
            /// </summary>
            /// <param name="i">The zero-based index of the part of R to calculate. Must not be zero (use function GetTopLeftR).</param>
            /// <returns></returns>
            public (Complex R00, Complex R01, Complex R10, Complex R11, Complex R20, Complex R21) GetRBlock(int i)
            {
                if (i == 0)
                {
                    throw new ArgumentOutOfRangeException(nameof(i));
                }
                // R(i,i)
                var R10 = GetDiagonalOfR(i);
                // R(i+1, i+1)
                var R21 = GetDiagonalOfR(i + 1);

                // r_{i-1,i} = (h_{i,i} - c_{i,i} * r_{i,i}) / c_{i,i-1}
                // R(i-1, i)
                var R00 = D[i - 1] * (GetDiagonalOfH(i) - C.GetDiagonal(i) * R10) / C.Rotations[i - 1].Sine;

                // r_{i,i+1} = (h_{i+1,i+1} - c_{i+1,i+1} * r_{i+1,i+1}) / c_{i+1,i}
                // R(i, i+1)
                var R11 = (GetDiagonalOfH(i + 1) - C.GetDiagonal(i + 1) * R21) / C.Rotations[i].Sine;

                // Using r_{i-2, i} = (h_{{i-1}, i} - c_{i-1, i-1} * r_{i-1, i} - c_{i-1, i} * r_{i,i}) / c_{i-1, i-2}
                // R(i-1, i+1)
                var R01 = D[i - 1] * (B.GetSuperdiagonal(i + 1) - C.GetDiagonal(i) * R11 - C.GetSuperdiagonal(i + 1) * R21) / C.Rotations[i - 1].Sine;

                R10 *= D[i];
                R21 *= D[i + 1];
                R11 *= D[i];

                return (R00, R01, R10, R11, 0, R21);
            }

            public (Complex a, Complex b, Complex c, Complex d) GetPart(int i)
            {
                if (i == 0)
                {
                    return GetTopLeft();
                }
                (var R00, var R01, var R10, var R11, var R20, var R21) = GetRBlock(i);

                if (i + 1 < Q.Rotations.Length)
                {
                    R21 = Q.Rotations[i + 1].Cosine * R21;
                }

                // Apply transformation which is active on rows i, i+1.
                (R10, R20) = Q.Rotations[i] * (R10, R20);
                (R11, R21) = Q.Rotations[i] * (R11, R21);

                // Apply transformation which is active on rows i-1, i.
                (_, R10) = Q.Rotations[i - 1] * (R00, R10);
                (_, R11) = Q.Rotations[i - 1] * (R01, R11);
                return (R10, R11, R20, R21);
            }

            public (Complex e1, Complex e2) GetWilkinsonShift()
            {
                (var AUpperLeft, var AUpperRight, var ALowerLeft, var ALowerRight) = GetPart(n - 1);
                return AmvwSolver.GetWilkinsonShift(AUpperLeft, AUpperRight, ALowerLeft, ALowerRight);
            }

            public void ApplyTransformationLeft(GivensRotation rotation)
            {
                var diagonal = Q.ApplyTransformationLeft(rotation, 0);
                for (var i = 0; i < D.Length && i < diagonal.Length; ++i)
                {
                    D[i] *= diagonal[i];
                }
                rotation.Cosine = 1;
                rotation.Sine = 0;
            }

            public void AbsorbTransformationInQ(GivensRotation transformation, int i)
            {
                if (i != n - 2)
                {
                    throw new ArgumentException("transformation only can be absorbed into last Q transformation.");
                }
                (var d, var e) = Q.Rotations[n - 2].Fusion(transformation);
                D[n - 2] *= d;
                D[n - 1] *= e;
                transformation.Cosine = 1;
                transformation.Sine = 0;
            }

            public (GivensRotation transformation, int i) Passthrough(int i, GivensRotation transformation)
            {
                // Turnover on B:
                (transformation, i) = B.ApplyTransformationRight(transformation, i);

                // Turnover on C:
                (transformation, i) = C.ApplyTransformationConjugate(transformation, i);

                // pass through diagonal
                (transformation, D[i], D[i + 1]) = GivensRotation.PassthroughRotationToLeft(D[i], D[i + 1], transformation);

                if (i == n - 2)
                {
                    // Fuse the transformation with the last rotation in Q.
                    AbsorbTransformationInQ(transformation, i);
                }
                else
                {
                    // Turnover on Q:
                    (transformation, i) = Q.ApplyTransformationRight(transformation, i);
                }
                return (transformation, i);
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

        public List<Complex> Solve(Polynomial polynomial, double precision)
        {
            return Solve(polynomial.Coefficients, precision);
        }

        public List<Complex> Solve(double[] coefficients, double precision, ShiftStrategy shiftStrategy = ShiftStrategy.Wilkinson)
        {
            return Solve(coefficients.Select(x => new Complex(x, 0)).ToArray(), precision, shiftStrategy).Item1;
        }

        public static void ApplyFrancisStep(AmvwFactorization factorization)
        {
            int n = factorization.n;
            (var A11, var A12, var A21, var A22) = factorization.GetTopLeft();
            (var Aul, var Aur, var All, var Alr) = factorization.GetPart(n - 2);

            (var shift1, var shift2) = GetWilkinsonShift(Aul, Aur, All, Alr);

            var shift = (shift1 - Alr).MagnitudeSquared() < (shift2 - Alr).MagnitudeSquared() ? shift1 : shift2;
            //var shift = 0;
            Console.WriteLine("Using shift: {0}", shift1);

            (var u1, _) = GivensRotation.Create(A11 - shift, A21);
            ApplyFrancisStep(u1, factorization);
        }

        /// <summary>
        /// Applies a similarity transformation on the matrix A using u1:
        /// un^H * ... * u1^H * A * u1 * ... * un .
        /// </summary>
        /// <param name="u1"></param>
        /// <param name="factorization"></param>
        public static void ApplyFrancisStep(GivensRotation u1, AmvwFactorization factorization)
        {
            var n = factorization.n;

            // fuse u1 with the first Q core transformation.
            //(var d, var e) = u1H.Fusion(factorization.Q.Rotations[0]);
            //factorization.D[0] *= d;
            //factorization.D[1] *= e;
            //factorization.Q.Rotations[0] = u1H;
            factorization.ApplyTransformationLeft(u1.CloneInverse());

            var transformation = u1;
            // i is equal to the active rows (i,i+1) the current rotation works on.
            for (int i = 0; i < n - 1;)
            {
                (transformation, i) = factorization.Passthrough(i, transformation);
                if (transformation.Cosine == Complex.One && transformation.Sine == 0)
                {
                    break;
                }
            }
        }

        public enum ShiftStrategy
        {
            Wilkinson,
            Rayleigh,
            None,
            Random
        }
        public static Complex[] ScaleCoefficients(Complex[] coefficients, Complex scaling)
        {
            int n = coefficients.Length;
            var scaledCoefficients = coefficients.Select((x, i) => x / scaling.Power(n - i));
            return scaledCoefficients.ToArray();
        }

        public static (List<Complex>, int) SolveScaled(Complex[] coefficients, double precision, ShiftStrategy shiftStrategy)
        {
            var scaling = coefficients[0].Power(1.0 / coefficients.Length);
            var scaledCoefficients = ScaleCoefficients(coefficients, scaling);
            (var roots, int iterations) = Solve(scaledCoefficients.ToArray(), precision, shiftStrategy);
            return (roots.Select(x => x * scaling).ToList(), iterations);
        }

        public static (List<Complex>, int) Solve(Complex[] coefficients, double precision, ShiftStrategy shiftStrategy)
        {
            // The companion matrix A is factorized into:
            // A = Q * C' * (B' + e1 * y')
            // where Q, C, B are unitary and stored using Givens rotations.
            var factorization = AmvwFactorization.Create(coefficients);
            return Solve(factorization, precision, shiftStrategy);
        }

        public static (List<Complex>, int) Solve(AmvwFactorization factorization, double precision = 1e-7, ShiftStrategy shiftStrategy = ShiftStrategy.Wilkinson)
        {
            var n = factorization.n;
            var roots = new List<Complex>();
            int j = 0;
            while(true)
            {

                if (factorization.Q.Rotations[n - 2].Sine < precision)
                {
                    (_, _, _, var Alr) = factorization.GetPart(n - 2);
                    roots.Add(Alr);
                    Console.WriteLine($"deflate root {Alr}");
                    n -= 1;

                    if (n <= 2)
                    {
                        (var A11, var A12, var A21, var A22) = factorization.GetTopLeft();
                        (var root1, var root2) = GetWilkinsonShift(A11, A12, A21, A22);
                        roots.Add(root1);
                        roots.Add(root2);
                        break;
                    }
                }
                var u1 = factorization.GetFirstTransformation(n, shiftStrategy);
                //Console.WriteLine($"sin: {factorization.Q.Rotations[n - 2].Sine}");
                ApplyFrancisStep(u1, factorization);
                ++j;
            }


            return (roots, j);
        }
    }
}
