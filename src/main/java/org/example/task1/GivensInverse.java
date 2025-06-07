package org.example.task1;

public class GivensInverse {
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
    public static double[][] createIdentityMatrix(int n) {
        double[][] identityMatrix = new double[n][n];
        for (int i = 0; i < n; ++i) {
            identityMatrix[i][i] = 1.0;
        }
        return identityMatrix;
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

//    public static double[][] inverse(double[][] input) {
//        int n = input.length;
//        double[][] A = copyMatrix(input);           // Исходная матрица, будем модифицировать
//        double[][] Inv = identityMatrix(n);         // Начинаем с единичной, превращаем в обратную
//
//        // Прямой ход: вращения Гивенса
//        for (int j = 0; j < n; j++) {
//            for (int i = n - 1; i > j; i--) {
//                if (A[i][j] != 0.0) {
//                    double a = A[i - 1][j];
//                    double b = A[i][j];
//                    double r = Math.hypot(a, b);
//                    double c = a / r;
//                    double s = -b / r;
//
//                    for (int k = 0; k < n; k++) {
//                        // Вращение строк в A
//                        double t1 = A[i - 1][k];
//                        double t2 = A[i][k];
//                        A[i - 1][k] = c * t1 - s * t2;
//                        A[i][k]     = s * t1 + c * t2;
//
//                        // То же вращение применяем к Inv
//                        t1 = Inv[i - 1][k];
//                        t2 = Inv[i][k];
//                        Inv[i - 1][k] = c * t1 - s * t2;
//                        Inv[i][k]     = s * t1 + c * t2;
//                    }
//                }
//            }
//        }
//
//        // Обратный ход: решаем A (верхнетреугольная) * X = Inv
//        for (int col = 0; col < n; col++) {
//            for (int i = n - 1; i >= 0; i--) {
//                double sum = Inv[i][col];
//                for (int j = i + 1; j < n; j++) {
//                    sum -= A[i][j] * Inv[j][col];
//                }
//                Inv[i][col] = sum / A[i][i];
//            }
//        }
//
//        return Inv;
//    }

    public static double[][] inverse(double[][] a) {
        int n = a.length;
        double[][] A = new double[n][2 * n];

        // Формируем расширенную матрицу: A | I
        for (int i = 0; i < n; i++) {
            System.arraycopy(a[i], 0, A[i], 0, n);      // Копия исходной матрицы
            A[i][n + i] = 1.0;                          // Единичная справа
        }

        double[] tempRow = new double[2 * n];

        // Прямой ход — зануляем ниже диагонали
        for (int j = 0; j < n; j++) {
            for (int i = n - 1; i > j; i--) {
                if (Math.abs(A[i][j]) > 1e-10) {
                    double a_ = A[i - 1][j];
                    double b = A[i][j];
                    double r = Math.hypot(a_, b);
                    double c = a_ / r;
                    double s = -b / r;

                    System.arraycopy(A[i - 1], 0, tempRow, 0, 2 * n);

                    for (int k = 0; k < 2 * n; k++) {
                        double t1 = tempRow[k];
                        double t2 = A[i][k];
                        A[i - 1][k] = c * t1 - s * t2;
                        A[i][k] = s * t1 + c * t2;
                    }
                }
            }
        }

        // Обратный ход — приводим к единичной
        for (int j = n - 1; j >= 0; j--) {
            // Нормируем строку
            double diag = A[j][j];
            for (int k = 0; k < 2 * n; k++) {
                A[j][k] /= diag;
            }

            // Зануляем выше диагонали
            for (int i = j - 1; i >= 0; i--) {
                double factor = A[i][j];
                for (int k = 0; k < 2 * n; k++) {
                    A[i][k] -= factor * A[j][k];
                }
            }
        }

        // Копируем правую часть в `a` (обратная матрица)
        for (int i = 0; i < n; i++) {
            System.arraycopy(A[i], n, a[i], 0, n);
        }

        return a;
    }

    // Печать
    private static void printMatrix(double[][] matrix) {
        for (double[] row : matrix) {
            for (double val : row) {
                System.out.printf("%10.5f", val);
            }
            System.out.println();
        }
    }

    // Тест
    public static void main(String[] args) {
        double[][] A = {
                {2, -1, 0},
                {-1, 2, -1},
                {0, -1, 2}
        };

        System.out.println("Original matrix A:");
        printMatrix(A);

        double[][] inv = inverse(A);
        System.out.println("\nInverse of A:");
        printMatrix(inv);

        double[][] product = matrixMultiplication(A, inv);
        System.out.println("\nA * inverse(A) ≈ I:");
        printMatrix(product);
    }
}