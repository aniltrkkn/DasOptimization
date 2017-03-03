/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package optimization;

import cern.colt.matrix.DoubleMatrix1D;
import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.impl.DenseDoubleMatrix1D;
import cern.colt.matrix.impl.DenseDoubleMatrix2D;
import cern.colt.matrix.impl.SparseDoubleMatrix2D;
import optimization.functionImplementation.ObjectiveFunction;
import optimization.functionImplementation.Options;
import solvers.NonlinearEquationSolver;

/**
 *
 * @author O. Anil Turkkan <turkkan.1@osu.edu> - anilturkkan.com
 */
public class Optimization {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        //
        final int n=100;
        ObjectiveFunction f = new ObjectiveFunction() {
            @Override
            public DoubleMatrix1D getF(DoubleMatrix1D x) {
                DoubleMatrix1D f = new DenseDoubleMatrix1D(n);
                for(int i=0;i<n/2;i++){
                  f.set(2*i, 10 * (x.get(2*i+1) - x.get(2*i) * x.get(2*i))) ; 
                  f.set(2*i+1,1-x.get(2*i));
                }
                return f;
            }

            @Override
            public DoubleMatrix1D getD(DoubleMatrix1D x) {
                return null;
            }

            @Override
            public DoubleMatrix2D getH(DoubleMatrix1D x) {
                DoubleMatrix2D H = new DenseDoubleMatrix2D(n,n);
                for(int i=0;i<n/2;i++){
                  H.set(2*i,2*i, -20*x.get(2*i)) ; 
                  H.set(2*i, 2*i+1,10);
                  H.set(2*i+1,2*i, -1);
                }

                return H;
            }
        };
        double[] initialGuess = new double[n];
        Options options = new Options(n);
        options.setAnalyticalHessian(true);
        NonlinearEquationSolver solver = new NonlinearEquationSolver(f, options, initialGuess);
        long startTime = System.nanoTime();
        for (int i = 0; i < 1e3; i++) {
            System.out.println(i);
            solver.solve(new DenseDoubleMatrix1D(initialGuess));
        }
        long endTime = System.nanoTime();
        System.out.println((endTime - startTime) / 1e9);
    }

}
