package org.example.task3;

public class Gauss {

    public static double[][] matrixMultiplication(double[][] a, double[][] b) {
        double[][] result = new double[a.length][a.length];
        for (int i = 0; i < a.length; ++i) {
            for (int j = 0; j < a[i].length; ++j) {
                for (int k = 0; k < a[i].length; ++k) {
                    result[i][j] += a[i][k] * b[k][j];
                }
            }
        }
        return result;
    }

    public static double[][] matrixDifference(double[][] a, double[][] b) {
        double[][] result = new double[a.length][a.length];
        for (int i = 0; i < a.length; i++) {
            for (int j = 0; j < a[i].length; j++) {
                result[i][j] = a[i][j] - b[i][j];
            }
        }
        return result;
    }


    public static void printMatrix(double[][] matrix) {
        for (double[] row : matrix) {
            for (double val : row) {
                System.out.printf("%10.32f ", val);
            }
            System.out.println();
        }
        System.out.println();
    }

    public static double[][] createIdentityMatrix(int n) {
        double[][] identityMatrix = new double[n][n];
        for (int i = 0; i < n; ++i) {
            identityMatrix[i][i] = 1.0;
        }
        return identityMatrix;
    }

    public static double[][] inverse(double[][] a) {
        int n = a.length;
        double[][] identityMatrix = createIdentityMatrix(n);

        for (int k = 0; k < n; ++k) {
            // Поиск строки с максимальным элементом в текущем столбце
            int maxRow = k;
            for (int i = k + 1; i < n; ++i) {
                if (Math.abs(a[i][k]) > Math.abs(a[maxRow][k])) {
                    maxRow = i;
                }
            }

            // Перестановка текущей строки с найденной строкой в обеих матрицах
            double[] tempA = a[k];
            a[k] = a[maxRow];
            a[maxRow] = tempA;

            double[] tempI = identityMatrix[k];
            identityMatrix[k] = identityMatrix[maxRow];
            identityMatrix[maxRow] = tempI;

            // Нормализация текущей строки
            double pivot = a[k][k];
            for (int j = 0; j < n; ++j) {
                a[k][j] /= pivot;
                identityMatrix[k][j] /= pivot;
            }

            // Обнуление текущего столбца в других строках
            for (int i = 0; i < n; ++i) {
                if (i != k) {
                    double factor = a[i][k];
                    for (int j = 0; j < n; ++j) {
                        a[i][j] -= factor * a[k][j];
                        identityMatrix[i][j] -= factor * identityMatrix[k][j];
                    }
                }
            }
        }
        printMatrix(identityMatrix);
        return identityMatrix;
    }
}