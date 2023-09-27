using System.Diagnostics.CodeAnalysis;
using System.Linq.Expressions;
using System.Numerics;
using System.Runtime.Intrinsics;

namespace TOOP_3sem
{
    interface IVector : IList<double> { }
    interface IParametricFunction
    {
        IFunction Bind(IVector parameters);
    }
    interface IFunction
    {
        double Value(IVector point);
    }
    interface IDifferentiableFunction : IFunction
    {
        // По параметрам исходной IParametricFunction
        IVector Gradient(IVector point);
    }
    interface IFunctional
    {
        double Value(IFunction function);
    }
    interface IDifferentiableFunctional : IFunctional
    {
        IVector Gradient(IFunction function);
    }
    interface IMatrix : IList<IList<double>>
    {

    }
    interface ILeastSquaresFunctional : IFunctional
    {
        IVector Residual(IFunction function);
        IMatrix Jacobian(IFunction function);
    }
    interface IOptimizator
    {
        IVector Minimize(IFunctional objective, IParametricFunction function, IVector initialParameters, IVector minimumParameters = default, IVector maximumParameters = default);
    }
    public class Vector : List<double>, IVector
    {
    }
    class LineFunction : IParametricFunction
    {
        class InternalLineFunction : IFunction
        {
            public double a, b;
            public double Value(IVector point) => a * point[0] + b;
        }
        public IFunction Bind(IVector parameters) => new InternalLineFunction() { a = parameters[0], b = parameters[1] };
    }
    class MyFunctional : IFunctional
    {
        public List<(double x, double y)> points;
        public double Value(IFunction function)
        {
            double sum = 0;
            foreach (var point in points)
            {
                var param = new Vector();
                param.Add(point.x);
                var s = function.Value(param) - point.y;
                sum += s * s;
            }
            return sum;
        }
    }
    class MyLineNFunction : IParametricFunction
    {
        class InternalLineFunction : IDifferentiableFunction
        {
            public IVector a;
            public double Value(IVector point)
            {
                double sum = 0;
                for (int i = 0; i < point.Count; i++)
                {
                    sum += a[i] * point[i];
                }
                sum += a.Last();
                return sum;
            }
            public IVector Gradient(IVector point)
            {
                IVector grad = new Vector();

                for (int i = 0; i < point.Count; i++)
                    grad.Add(point[i]);
                grad.Add(0);
                return grad;
            }
        }
        public IFunction Bind(IVector parameters) => new InternalLineFunction() { a = parameters };

    }
    class MyPolynom : IParametricFunction
    {
        class InternalLineFunction : IFunction
        {
            public IVector a;
            public double Value(IVector point)
            {
                double sum = 0;
                for (int i = a.Count - 1; i >= 0; i--)
                    sum += a[i] * Math.Pow(point[0], i);
                return sum;
            }
        }
        public IFunction Bind(IVector parameters) => new InternalLineFunction() { a = parameters };

    }
    class MyL1Functional : IDifferentiableFunctional
    {
        public List<(IVector point, double f)> pointsAndF;
        public double Value(IFunction function)
        {
            double sum = 0;
            foreach (var pointAndF in pointsAndF)
                sum += Math.Abs(function.Value(pointAndF.point) - pointAndF.f);
            return sum;
        }
        public IVector Gradient(IFunction function)
        {
            Vector grad = new Vector();
            IDifferentiableFunction fun;
            if (function is IDifferentiableFunction)
            {
                double sign = 0;

                fun = function as IDifferentiableFunction;

                for (int i = 0; i < pointsAndF[0].point.Count; i++)
                    grad.Add(0);
                grad.Add(1);

                for (int i = 0; i < pointsAndF.Count; i++)
                {
                    sign = Math.Sign(function.Value(pointsAndF[i].point) - pointsAndF[i].f);
                    for (int j = 0; j <= pointsAndF[i].point.Count; j++)
                        grad[j] += sign * fun.Gradient(pointsAndF[i].point)[j];
                }
                return grad;
            }
            else
            {
                Console.WriteLine("ERROR!!! This Function is not Differentiable!");
                return grad;
            }
        }
    }
    class MyL2Functional : IDifferentiableFunctional
    {
        public List<(IVector point, double f)> pointsAndF;
        public double Value(IFunction function)
        {
            double sum = 0, s = 0;
            foreach (var pointAndF in pointsAndF)
            {
                s = function.Value(pointAndF.point) - pointAndF.f;
                sum += s * s;
            }
            return sum;
        }
        public IVector Gradient(IFunction function)
        {
            Vector grad = new Vector();
            IDifferentiableFunction fun;
            if (function is IDifferentiableFunction)
            {
                fun = function as IDifferentiableFunction;

                for (int i = 0; i < pointsAndF[0].point.Count; i++)
                    grad.Add(0);
                grad.Add(1);
                for (int i = 0; i < pointsAndF.Count; i++)
                {
                    for (int j = 0; j <= pointsAndF[i].point.Count; j++)
                        grad[j] += 2 * (fun.Value(pointsAndF[i].point) - pointsAndF[i].f) * fun.Gradient(pointsAndF[i].point)[j];
                }
                return grad;
            }
            else
            {
                Console.WriteLine("ERROR!!! This Function is not Differentiable!");
                return grad;
            }
        }
    }
    class MinimizerGradient : IOptimizator
    {
        public int MaxIter = 100000;
        private double lambda = 1e-3, eps = 1e-8;
        private double MakeSimplefx(double x, IDifferentiableFunctional objective, IParametricFunction function, IVector parameters)
        {
            Vector buffer = new();
            int n = parameters.Count;

            for (int i = 0; i < n; i++)
            {
                buffer.Add(parameters[i] - x * objective.Gradient(function.Bind(parameters))[i]);
            }

            return objective.Value(function.Bind(buffer));
        }
        double GoldenSelection(double a, double b, double eps, IDifferentiableFunctional objective, IParametricFunction function, IVector parameters)
        {
            const double fi = 1.6180339887;
            double x1, x2;
            double y1, y2;

            x1 = b - ((b - a) / fi);
            x2 = a + ((b - a) / fi);

            y1 = MakeSimplefx(x1, objective, function, parameters);
            y2 = MakeSimplefx(x2, objective, function, parameters);
            while (Math.Abs(b - a) > eps)
            {
                if (y1 <= y2)
                {
                    b = x2;
                    x2 = x1;
                    x1 = b - ((b - a) / fi);
                    y2 = y1;
                    y1 = MakeSimplefx(x1, objective, function, parameters);
                }
                else
                {
                    a = x1;
                    x1 = x2;
                    x2 = a + ((b - a) / fi);
                    y1 = y2;
                    y2 = MakeSimplefx(x2, objective, function, parameters);
                }
            }

            return (a + b) / 2;
        }
        public IVector Minimize(IFunctional objective, IParametricFunction function, IVector initialParameters, IVector minimumParameters = null, IVector maximumParameters = null)
        {
            var param = new Vector();
            var minparam = new Vector();
            foreach (var p in initialParameters) param.Add(p);
            foreach (var p in initialParameters) minparam.Add(p);

            IDifferentiableFunctional obj;
            if (objective is IDifferentiableFunctional)
            {
                obj = objective as IDifferentiableFunctional;

                var f = objective.Value(function.Bind(param));
                int i;
                for (i = 0; i < MaxIter; i++)
                {
                    lambda = GoldenSelection(0, 2, 1e-9, obj, function, minparam);

                    for (int j = 0; j < param.Count; j++)
                        param[j] = minparam[j] - lambda * obj.Gradient(function.Bind(minparam))[j];

                    f = objective.Value(function.Bind(param));

                    if (f < eps)
                        break;
                    else
                        for (int j = 0; j < param.Count; j++)
                            minparam[j] = param[j];
                    //if (i == MaxIter - 1)
                    //    Console.WriteLine("Maxiter! f = " + f.ToString());
                }
                Console.WriteLine("Gradient: Functional = " + f.ToString() + " Iter = " + i.ToString());
            }
            return param;
        }
    }
    class MinimizerMonteCarlo : IOptimizator
    {
        public int MaxIter = 100000;

        public IVector Minimize(IFunctional objective, IParametricFunction function, IVector initialParameters, IVector minimumParameters = null, IVector maximumParameters = null)
        {
            var param = new Vector();
            var minparam = new Vector();
            foreach (var p in initialParameters) param.Add(p);
            foreach (var p in initialParameters) minparam.Add(p);
            var fun = function.Bind(param);
            var currentmin = objective.Value(fun);
            var rand = new Random(0);
            int i;
            for (i = 0; i < MaxIter; i++)
            {
                for (int j = 0; j < param.Count; j++) param[j] = (rand.NextDouble() - 0.5) * 20;
                var f = objective.Value(function.Bind(param));
                if (f < currentmin)
                {
                    currentmin = f;
                    for (int j = 0; j < param.Count; j++) minparam[j] = param[j];
                }
            }
            Console.WriteLine("MonteCarlo: Functional = " + currentmin.ToString());

            return minparam;
        }
    }
    class Program
    {
        static void Main(string[] args)
        {
            var optimizer = new MinimizerMonteCarlo();
            //var optimizer = new MinimizerGradient();
            var initial = new Vector();
            Console.Write("Пространство: ");
            int n = int.Parse(Console.ReadLine());

            Console.Write("Количество точек: ");
            int m = int.Parse(Console.ReadLine());
            
            for (int i = 0; i <= n; i++)
                initial.Add(1);

            List<(IVector x, double y)> points = new();
            
            IVector point;

            for (int i = 0; i < m; i++)
            {
                point = new Vector();
                var str = Console.ReadLine().Split();
                for (int j = 0; j < n; j++)
                { 
                    point.Add(double.Parse(str[j]));
                }

                points.Add((point, double.Parse(str.Last())));
            }
            //var functinal = new MyFunctional() { points = points };
            var functinal = new MyL1Functional() { pointsAndF = points };
            
            //var fun = new LineFunction();
            var fun = new MyLineNFunction();

            var res = optimizer.Minimize(functinal, fun, initial);

            for (int i = 0; i <=n; i++)
                Console.WriteLine($"{res[i]} ");
        }
    }
}