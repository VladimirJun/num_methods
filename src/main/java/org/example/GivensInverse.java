package org.example;


public class GivensInverse {
    public static boolean inverse(double[][] matrix) {
        int n = matrix.length;

        // Проверка, что матрица квадратная
        if (n == 0 || matrix.length != matrix[0].length) {
            System.out.println("Матрица должна быть квадратной.");
            return false;
        }

        // Прямое преобразование в верхнетреугольную матрицу
        for (int j = 0; j < n; j++) {
            for (int i = j + 1; i < n; i++) {
                if (matrix[i][j] != 0) {
                    double r = Math.sqrt(matrix[j][j] * matrix[j][j] + matrix[i][j] * matrix[i][j]);
                    double c = matrix[j][j] / r;
                    double s = -matrix[i][j] / r;

                    // Вращение строк
                    for (int k = 0; k < n; k++) {
                        double temp = matrix[j][k];
                        matrix[j][k] = c * temp - s * matrix[i][k];
                        matrix[i][k] = s * temp + c * matrix[i][k];
                    }
                }
            }
        }

        // Обращение верхнетреугольной матрицы (матрица остается верхнетреугольной)
        for (int i = n - 1; i >= 0; i--) {
            if (Math.abs(matrix[i][i]) < 1e-50) {
                System.out.println("Матрица вырожденная.");
                return false;
            }
            matrix[i][i] = 1.0 / matrix[i][i];
            for (int j = i - 1; j >= 0; j--) {
                double sum = 0;
                for (int k = j + 1; k < n; k++) {
                    sum += matrix[j][k] * matrix[k][i];
                }
                matrix[j][i] = -sum / matrix[j][j];
            }
        }

        // Обратное применение вращений Гивенса
        for (int j = n - 1; j > 0; j--) {
            for (int i = j - 1; i >= 0; i--) {
                double r = Math.sqrt(matrix[j][j] * matrix[j][j] + matrix[i][j] * matrix[i][j]);
                double c = matrix[j][j] / r;
                double s = matrix[i][j] / r;

                for (int k = 0; k < n; k++) {
                    double temp = matrix[j][k];
                    matrix[k][j] = c * temp - s * matrix[k][i];
                    matrix[k][i] = s * temp + c * matrix[k][i];
                }
            }
        }

        return true;
    }

    public void printMatrix(double[][] matrix) {
        for (double[] row : matrix) {
            for (double elem : row) {
                System.out.printf("%8.4f ", elem);
            }
            System.out.println();
        }
    }


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


    public static double[][] createIdentityMatrix(int n) {
        double[][] identityMatrix = new double[n][n];
        for (int i = 0; i < n; ++i) {
            identityMatrix[i][i] = 1.0;
        }
        return identityMatrix;
    }

    public void main(String[] args) {
        double[][] matrix = {
                {4, 7},
                {2, 6},
        };

        System.out.println("Исходная матрица:");
        printMatrix(matrix);

        if (inverse(matrix)) {
            System.out.println("Обратная матрица:");
            printMatrix(matrix);
        } else {
            System.out.println("Матрица необратима.");
        }
    }
}