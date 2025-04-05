package org.example;

import org.junit.Test;

/**
 * Unit test for simple App.
 */
public class GenTest
{

	private int N = 3;
	private double ALPHA = 1.;
	private double BETA  = 1.e+3;

	@Test
	public void GenTest()
	{
		int n = N;
		double alpha = ALPHA;
		double beta  = BETA;

		double[][] a = new double[n][];
		for (int i = 0; i < n; i++)	a[i] = new double[n];

		double[][] a_inv = new double[n][];
		for (int i = 0; i < n; i++)	a_inv[i] = new double[n];

		Gen g = new Gen();

		g.myGen ( a, a_inv, n, alpha, beta, 1, 2, 0, 1 ); // симметричная
	//	g.mygen ( a, a_inv, n, alpha, beta, 1, 2, 1, 1 ); //проостой структуры
	//	g.mygen ( a, a_inv, n, alpha, beta, 0, 0, 2, 1 ); //жорданова клетка

		g.printMatrix(a,n);
		g.printMatrix(a_inv,n);
	}

}
