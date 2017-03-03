/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package optimization;
import optimization.functionImplementation.ObjectiveFunction;
import optimization.functionImplementation.Options;
import org.ejml.data.DenseMatrix64F;
import org.ejml.ops.CommonOps;
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
        final int n=500;
        ObjectiveFunction f = new ObjectiveFunction() {
            @Override
            public DenseMatrix64F getF(DenseMatrix64F x) {
                DenseMatrix64F f = new DenseMatrix64F(n,1);
                for(int i=0;i<n/2;i++){
                  f.set(2*i,0, 10 * (x.get(2*i+1,0) - x.get(2*i,0) * x.get(2*i,0))) ; 
                  f.set(2*i+1,0,1-x.get(2*i,0));
                }
                return f;
            }

            @Override
            public DenseMatrix64F getD(DenseMatrix64F x) {
                return null;
            }

            @Override
            public DenseMatrix64F getH(DenseMatrix64F x) {
                DenseMatrix64F H = new DenseMatrix64F(n,n);
                for(int i=0;i<n/2;i++){
                  H.set(2*i,2*i, -20*x.get(2*i)) ; 
                  H.set(2*i, 2*i+1,10);
                  H.set(2*i+1,2*i, -1);
                }

                return H;
            }
        };
        DenseMatrix64F initialGuess = new DenseMatrix64F(n,1);
        CommonOps.fill(initialGuess,0.0);
        Options options = new Options(n);
        options.setAnalyticalHessian(true);
        NonlinearEquationSolver solver = new NonlinearEquationSolver(f, options);
        long startTime = System.nanoTime();
        for (int i = 0; i < 1e1; i++) {
            System.out.println(i);
            solver.solve(new DenseMatrix64F(initialGuess));
        }
        long endTime = System.nanoTime();
        System.out.println((endTime - startTime) / 1e9);
        System.out.println(solver.getX());
        //System.out.println(solver.getFx());
        //System.out.println(solver.getJacobian());
    }

}
