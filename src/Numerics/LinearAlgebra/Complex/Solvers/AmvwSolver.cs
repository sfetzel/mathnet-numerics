// <copyright file="BiCgStab.cs" company="Math.NET">
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

namespace MathNet.Numerics.LinearAlgebra.Complex.Solvers
{
    using Complex = System.Numerics.Complex;
    public class AmvwSolver
    {
        public class GivensRotation
        {
            public Complex Cosine { get; set; }
            public double Sine { get; set; }

            public static (Complex a, Complex b) operator *(GivensRotation rotation, (Complex a, Complex b) tuple)
            {
                return (tuple.a * rotation.Cosine + tuple.b * rotation.Sine, tuple.b * rotation.Cosine - tuple.a * rotation.Sine);
            }

            public GivensRotation(Complex cosine, double sine)
            {
                Cosine = cosine;
                Sine = sine;
            }

            public static (GivensRotation Rotation, double norm) Create(Complex a, Complex b)
            {
                var bNorm = b.Magnitude();
                var alpha = Math.Sqrt(a.MagnitudeSquared() + bNorm * bNorm);
                // the rotation must be described with a complex cosine
                // and a real sine. Basically this is done by factoring out the sine part.

                // calculating cosine and sine.
                var scaledA = a / alpha;
                var scaledB = b / alpha;

                // rotation would be given by (scaledA, scaledB). Factor out scaledB'.
                var cosine = scaledA * (b / bNorm).Conjugate() / alpha;
                var sine = bNorm / alpha;
                return (new GivensRotation(cosine, sine), alpha);
            }
            
            public GivensRotation()
            {
                Cosine = 1;
                Sine = 0;
            }
        }

        public class UnitaryMatrix
        {
            public Dictionary<int, GivensRotation> CoreTransformations = new Dictionary<int, GivensRotation>();
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
            int n = coefficients.Length;

            // The companion matrix A is factorized into:
            // A = Q * C * (B' + e1 * y')
            // where Q, C, B are unitary and stored using Givens rotations.

            GivensRotation[] Q = new GivensRotation[n];

            for(int i = 0; i < n-1; ++i)
            {
                Q[i] = new GivensRotation(0, -1);
            }
            Q[n-1] = new GivensRotation(1, 0);



            return null;
        }
    }
}