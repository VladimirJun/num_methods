package org.example;

import org.example.task1.Generator;
import org.example.task2.RichardsonSolver;

public class MatrixErrorMeasurement {
    public static void main(String[] args) {
        Generator generator = new Generator();

        int n = 100;
        double alpha, beta, aNorm, aInvNorm, vA, errorNorm, dzeta, rNorm;
        int cycles = 15;

        System.out.println("-".repeat(121));
        System.out.println("|    alpha     |     beta     |   norm A     |  norm A^-1   |    obusl     |   norm err   |    dzeta     |   nevyazka   |");
        System.out.println("-".repeat(121));

        for (int i = 0; i < cycles; i++) {
            alpha = Math.pow(10, -(i + 1));
            beta = 1;

            double[][] A = new double[n][n];
            double[][] Ainv = new double[n][n];
            generator.myGen(A, Ainv, n, alpha, beta, 1, 2, 0, 1);

            double[] xExpected = new double[n];
            for (int j = 0; j < n; j++) xExpected[j] = 1.0;

            double[] b = new double[n];
            for (int j = 0; j < n; j++) {
                for (int k = 0; k < n; k++) {
                    b[j] += A[j][k] * xExpected[k];
                }
            }

            double[] xComputed = RichardsonSolver.solveChebyshev(A, b, 10);

            aNorm = generator.matrixInfNorm(A, n);
            aInvNorm = generator.matrixInfNorm(Ainv, n);
            vA = aNorm * aInvNorm;
            errorNorm = vectorErrorNorm(xExpected, xComputed);
            dzeta = errorNorm / vectorNorm(xExpected);

            double[] Ax = new double[n];
            for (int j = 0; j < n; j++)
                for (int k = 0; k < n; k++)
                    Ax[j] += A[j][k] * xComputed[k];
            double[] r = new double[n];
            for (int j = 0; j < n; j++) r[j] = Ax[j] - b[j];
            rNorm = vectorNorm(r);

            System.out.printf("| %e | %2.10f | %2.10f | %e | %e | %e | %e | %e |\n",
                    alpha, beta, aNorm, aInvNorm, vA, errorNorm, dzeta, rNorm);
        }

        System.out.println("-".repeat(121));
    }

    private static double vectorNorm(double[] v) {
        double sum = 0;
        for (double val : v) sum += Math.abs(val);
        return sum;
    }

    private static double vectorErrorNorm(double[] expected, double[] actual) {
        double max = 0;
        for (int i = 0; i < expected.length; i++) {
            max = Math.max(max, Math.abs(expected[i] - actual[i]));
        }
        return max;
    }
}
