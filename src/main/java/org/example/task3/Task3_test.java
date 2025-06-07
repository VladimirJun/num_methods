package org.example.task3;

public class Task3_test {
    public static void main(String[] args) {
        test3x3Matrix();
        test4x4Matrix();
    }

    public static void test3x3Matrix() {
        double[][] matrix = {
                {4, 1, 1},
                {1, 3, 0},
                {1, 0, 2}
        };

        double expectedMinEigenvalue = 1.46; //из калькулятора
        double actual = DegMethodTask3.findMinimalEigenvalue(matrix, 1e-10, 1000);
        System.out.printf("3x3 matrix test:\nExpected ≈ %.6f, Actual = %.6f\n\n", expectedMinEigenvalue, actual);
    }

    public static void test4x4Matrix() {
        double[][] matrix = {
                {5, 2, 0, 0},
                {2, 6, 2, 0},
                {0, 2, 7, 2},
                {0, 0, 2, 8}
        };

        double expectedMinEigenvalue = 3.0; // примерно
        double actual = DegMethodTask3.findMinimalEigenvalue(matrix, 1e-10, 1000);
        System.out.printf("4x4 matrix test:\nExpected ≈ %.6f, Actual = %.6f\n\n", expectedMinEigenvalue, actual);
    }
}
