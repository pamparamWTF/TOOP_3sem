using System.Diagnostics.CodeAnalysis;
using System.Linq.Expressions;

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

                for(int i = 0; i < point.Count; i++)
                    grad.Add(point[i]);
                grad.Add(1);
                return grad;
            }
        }
        public IFunction Bind(IVector parameters) => new InternalLineFunction() { a = parameters };

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
    class MyL1NormDifferentiableFunctional : IDifferentiableFunctional
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

                for (int i = 0; i < pointsAndF.Count; i++)
                {
                    sign = Math.Sign(function.Value(pointsAndF[i].point) - pointsAndF[i].f);
                    for (int j = 0; j < pointsAndF[i].point.Count; j++)
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
    class MinimizerGradient : IOptimizator
    {
        public int MaxIter = 100000;
        private double lambda = 0.1;
        public IVector Minimize(IFunctional objective, IParametricFunction function, IVector initialParameters, IVector minimumParameters = null, IVector maximumParameters = null)
        {
            var param = new Vector();
            var minparam = new Vector();
            foreach (var p in initialParameters) param.Add(p);
            foreach (var p in initialParameters) minparam.Add(p);
            var fun = function.Bind(param);
            
            //var currentmin = objective.Value(fun);

            // засунуть в тру катч
            IDifferentiableFunctional obj;
            if (objective is IDifferentiableFunctional)
            {
                obj = objective as IDifferentiableFunctional;

                for (int i = 0; i < MaxIter; i++)
                {
                    for (int j = 0; j < param.Count; j++)
                    {
                        param[j] = minparam[j] - lambda * obj.Gradient(function.Bind(minparam))[j];
                    }
                    var f = objective.Value(function.Bind(param));

                }
                return param;
            }
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
            for (int i = 0; i < MaxIter; i++)
            {
                for (int j = 0; j < param.Count; j++) param[j] = rand.NextDouble();
                var f = objective.Value(function.Bind(param));
                if (f < currentmin)
                {
                    currentmin = f;
                    for (int j = 0; j < param.Count; j++) minparam[j] = param[j];
                }
            }
            return minparam;
        }
    }

    class Program
    {
        static void Main(string[] args)
        {
            var optimizer = new MinimizerMonteCarlo();
            var initial = new Vector();
            initial.Add(1);
            initial.Add(1);
            int n = int.Parse(Console.ReadLine());
            List<(double x, double y)> points = new();
            for (int i = 0; i < n; i++)
            {
                var str = Console.ReadLine().Split();
                points.Add((double.Parse(str[0]), double.Parse(str[1])));
            }
            var functinal = new MyFunctional() { points = points };
            var fun = new LineFunction();

            var res = optimizer.Minimize(functinal, fun, initial);
            Console.WriteLine($"a={res[0]},b={res[1]}");
        }
    }
}