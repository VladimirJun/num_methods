package org.example.task2;

public class RichardsonSolver {
    private static final int MAX_ITER = 1000;
    private static final double EPSILON = 1e-10;

    public static double[] solveChebyshev(double[][] A, double[] b, int s) {
        int n = b.length;

        double alphaMax = powerIteration(A);
        double alphaMin = inversePowerIteration(A);

        double[] tau = new double[s];
        for (int j = 0; j < s; j++) {
            tau[j] = 2.0 / (alphaMax + alphaMin + (alphaMax - alphaMin) * Math.cos(Math.PI * (2 * j + 1) / (2.0 * s)));
        }

        double[] x = new double[n];

        for (int iter = 0; iter < MAX_ITER; iter++) {
            for (int j = 0; j < s; j++) {
                // Вычисляем T_j * x и прибавляем τ_j * b
                double[] Tx = new double[n];
                for (int i = 0; i < n; i++) {
                    Tx[i] = (1.0 - tau[j] * A[i][i]) * x[i];
                    for (int k = 0; k < n; k++) {
                        if (i != k) {
                            Tx[i] -= tau[j] * A[i][k] * x[k];
                        }
                    }
                    Tx[i] += tau[j] * b[i];
                }

                // Обновляем x
                System.arraycopy(Tx, 0, x, 0, n);

                // Вычисляем невязку r = A*x - b
                double[] Ax = multiply(A, x);
                double[] r = new double[n];
                for (int i = 0; i < n; i++) {
                    r[i] = Ax[i] - b[i];
                }
                if (norm(r) < EPSILON) return x;
            }
        }

        return x;
    }

    private static double powerIteration(double[][] A) {
        int n = A.length;
        double[] b = new double[n];
        for (int i = 0; i < n; i++) b[i] = 1;

        double lambdaOld = 0;
        for (int iter = 0; iter < 1000; iter++) {
            double[] Ab = multiply(A, b);
            double lambda = dot(Ab, b) / dot(b, b);
            double norm = Math.sqrt(dot(Ab, Ab));
            for (int i = 0; i < n; i++) b[i] = Ab[i] / norm;

            if (Math.abs(lambda - lambdaOld) < 1e-8) break;
            lambdaOld = lambda;
        }
        return lambdaOld;
    }

    private static double inversePowerIteration(double[][] A) {
        int n = A.length;
        double[] b = new double[n];
        for (int i = 0; i < n; i++) b[i] = 1;

        double lambdaOld = 0;
        for (int iter = 0; iter < 1000; iter++) {
            double[] x = solveGaussian(A, b);
            double norm = Math.sqrt(dot(x, x));
            for (int i = 0; i < n; i++) b[i] = x[i] / norm;
            double[] Ax = multiply(A, b);
            double lambda = dot(Ax, b) / dot(b, b);

            if (Math.abs(lambda - lambdaOld) < 1e-8) break;
            lambdaOld = lambda;
        }
        return lambdaOld;
    }

    private static double[] multiply(double[][] A, double[] x) {
        int n = x.length;
        double[] result = new double[n];
        for (int i = 0; i < n; i++) {
            result[i] = 0;
            for (int j = 0; j < n; j++) {
                result[i] += A[i][j] * x[j];
            }
        }
        return result;
    }

    private static void multiply(double[][] A, double[] x, double[] result) {
        int n = x.length;
        for (int i = 0; i < n; i++) {
            result[i] = 0;
            for (int j = 0; j < n; j++) {
                result[i] += A[i][j] * x[j];
            }
        }
    }

    private static double[] solveGaussian(double[][] A, double[] b) {
        int n = b.length;
        double[][] a = new double[n][n];
        double[] rhs = b.clone();

        // Deep copy A to a
        for (int i = 0; i < n; i++)
            a[i] = A[i].clone();

        for (int i = 0; i < n; i++) {
            int maxRow = i;
            for (int k = i + 1; k < n; k++)
                if (Math.abs(a[k][i]) > Math.abs(a[maxRow][i]))
                    maxRow = k;

            double[] tmp = a[i];
            a[i] = a[maxRow];
            a[maxRow] = tmp;
            double t = rhs[i];
            rhs[i] = rhs[maxRow];
            rhs[maxRow] = t;

            for (int k = i + 1; k < n; k++) {
                double factor = a[k][i] / a[i][i];
                for (int j = i; j < n; j++)
                    a[k][j] -= factor * a[i][j];
                rhs[k] -= factor * rhs[i];
            }
        }

        double[] x = new double[n];
        for (int i = n - 1; i >= 0; i--) {
            x[i] = rhs[i];
            for (int j = i + 1; j < n; j++)
                x[i] -= a[i][j] * x[j];
            x[i] /= a[i][i];
        }
        return x;
    }

    private static double dot(double[] a, double[] b) {
        double sum = 0;
        for (int i = 0; i < a.length; i++) sum += a[i] * b[i];
        return sum;
    }

    private static double norm(double[] x) {
        return Math.sqrt(dot(x, x));
    }

    public static void main(String[] args) {
        double[][] A = {
                {4, 1, 2},
                {1, 3, 1},
                {2, 1, 5}
        };
        // Вектор правой части b
        double[] b = {7, 8, 10};

        // Решение методом Ричардсона с Чебышевскими параметрами
        double[] x = solveChebyshev(A, b, 5);

        // Вывод результатов
        System.out.println("Решение системы A*x = b:");
        for (int i = 0; i < x.length; i++) {
            System.out.printf("x[%d] = %.10f%n", i, x[i]);
        }
    }
}
