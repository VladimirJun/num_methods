package org.example;
import static org.junit.jupiter.api.Assertions.assertEquals;

import org.junit.jupiter.api.Test;

public class GivensInverseTest {

    private static final double EPSILON = 1e-9;

    @Test
    public void testInverse() {
        double[][] matrix = {
                {4, 7, 1},
                {2, 6, 2},
                {7, 8, 9},
        };
        double[][] expectedInverse = {
                {0.004388261531064058, 0.127502469925538,-0.05689095361581528},
                {-0.14320224437179235, 0.190180751870142,0.13213030970940484},
                {-0.025281671644413806, 0.191217119083307,0.2620461449544067}
        };

        GivensInverse.inverse(matrix);
        assertMatrixEquals(expectedInverse, matrix);
    }

    @Test
    public void testIdentityMatrix() {
        double[][] identity = {
                {1, 0, 0},
                {0, 1, 0},
                {0, 0, 1}
        };

    }

    private void assertMatrixEquals(double[][] expected, double[][] actual) {
        assertEquals(expected.length, actual.length, "Row count mismatch");
        assertEquals(expected[0].length, actual[0].length, "Column count mismatch");

        for (int i = 0; i < expected.length; i++) {
            for (int j = 0; j < expected[i].length; j++) {
                assertEquals(expected[i][j], actual[i][j], GivensInverseTest.EPSILON, "Mismatch at [" + i + "][" + j + "]");
            }
        }
    }
}
