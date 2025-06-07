package org.example.task3;

import java.util.Random;

public class DegMethodTask3 {

    // Метод для нахождения минимального по модулю собственного значения
    public static double findMinimalEigenvalue(double[][] a, double eps, int maxIterations) {
        int n = a.length;

        // Создаем копию матрицы A, так как Gauss.inverse изменяет входную матрицу
        double[][] aCopy = new double[n][n];
        for (int i = 0; i < n; i++) {
            System.arraycopy(a[i], 0, aCopy[i], 0, n);
        }

        // Шаг 1: Обратная матрица
        double[][] inverse = Gauss.inverse(aCopy);

        // Шаг 2: Степенной метод для нахождения максимального собственного значения A⁻¹
        double[] x = new double[n];
        for (int i = 0; i < n; i++) x[i] = 1.0;

        double lambdaOld = 0;
        for (int iter = 0; iter < maxIterations; iter++) {
            // Умножение A⁻¹ * x
            double[] y = new double[n];
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < n; j++) {
                    y[i] += inverse[i][j] * x[j];
                }
            }

            // Нормализация и расчет собственного значения
            double norm = 0;
            for (double v : y) norm += v * v;
            norm = Math.sqrt(norm);
            for (int i = 0; i < n; i++) x[i] = y[i] / norm;

            double lambda = 0;
            for (int i = 0; i < n; i++) lambda += x[i] * y[i];

            if (Math.abs(lambda - lambdaOld) < eps) break;
            lambdaOld = lambda;
        }

        // Шаг 3: Минимальное собственное значение исходной матрицы
        return 1.0 / lambdaOld;
    }
}


