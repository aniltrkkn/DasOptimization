/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package solverImplementation;

import optimization.functionImplementation.Options;
import org.ejml.data.DenseMatrix64F;
import org.ejml.ops.CommonOps;
import org.ejml.ops.NormOps;
import solvers.Solver;

/**
 *
 * @author turkkan.1
 */
public class LineSearch{

    //solver status
    private static final int STATUS_ONGOING = 0;
    private static final int STATUS_FAILED = 1;
    private static final int STATUS_INSIDE_FUNCTION = 2;
    private static int solverStatus;
    private static boolean maxStepTaken;

    public static boolean getSolverStatus() {
        return solverStatus != STATUS_FAILED;
    }

    public static boolean isMaxStepTaken() {
        return maxStepTaken;
    }

    //xp=x+lambda*sn such that f(xp) <= f(x) + alpha*lambda*g^T*p
    public static DenseMatrix64F lineSearch(DenseMatrix64F g, DenseMatrix64F sn, DenseMatrix64F x, Options solverOptions, Solver solver) {
        //maximum step taken
        maxStepTaken = false;
        //solver status
        solverStatus = STATUS_INSIDE_FUNCTION;
        //alpha
        double alpha = 1e-4;
        //norm of sn*Sx
        DenseMatrix64F dummySn = new DenseMatrix64F(sn);
        CommonOps.elementMult(dummySn, solverOptions.getTypicalX());
        double newtonLength = NormOps.normP2(dummySn);
        if (newtonLength > solverOptions.getMaxStep()) {
            //newton step xp=x+sn is longer than maximum allowed
            for (int i = 0; i < sn.numRows; i++) {
                sn.set(i, sn.get(i) * (solverOptions.getMaxStep() / newtonLength));
            }
            newtonLength = solverOptions.getMaxStep();
        }
        //initial slope
        double initialSlope = CommonOps.dot(g, sn);
        //relative length of sn as calculated in the stopping routine
        double relativeLength = Double.MIN_VALUE;
        for (int i = 0; i < sn.numRows; i++) {
            relativeLength = Math.max(relativeLength, Math.abs(sn.get(i)) / Math.max(Math.abs(x.get(i)), 1 / solverOptions.getTypicalX().get(i)));
        }
        //minimum allowable step length
        double minLambda = solverOptions.getStepTolerance() / relativeLength;
        /* find the new x values*/
        DenseMatrix64F xPlus = new DenseMatrix64F(x.numRows, x.numCols);
        double lambda = 1.0;
        double initialFunctionNorm = solver.functionNorm(x);
        double lambdaPrevious = 1.0;
        double previousFunctionNorm = 0.0;
        //loop to check whether xp=x+lambda*sn is satisfactory (generate new lambda if required)
        while (solverStatus == STATUS_INSIDE_FUNCTION) {
            //new x values
            for (int i = 0; i < x.numRows; i++) {
                xPlus.set(i, 0, x.get(i, 0) + lambda * sn.get(i, 0));
            }
            //new function norm
            double newFunctionNorm = solver.functionNorm(xPlus);
            if (newFunctionNorm <= initialFunctionNorm + alpha * lambda) {
                //satisfactory new x found
                solverStatus = STATUS_ONGOING;
                if (lambda == 1.0 && newtonLength > 0.99 * solverOptions.getMaxStep()) {
                    maxStepTaken = true;
                }
            } else if (lambda < minLambda) {
                //no satisfactory new x can be found different than initial x
                solverStatus = STATUS_FAILED;
                xPlus = x.copy();
            } else {
                //reduce lambda
                double lambdaTemp;
                if (lambda == 1) {
                    // first backtrack(quadratic fit)
                    lambdaTemp = -initialSlope / (2 * (newFunctionNorm - initialFunctionNorm - initialSlope));
                } else {
                    // cubic fit
                    DenseMatrix64F a = new DenseMatrix64F(new double[][]{{1 / (lambda * lambda), -1 / (lambdaPrevious * lambdaPrevious)}, {-lambdaPrevious / (lambda * lambda), lambda / (lambdaPrevious * lambdaPrevious)}});
                    DenseMatrix64F b = new DenseMatrix64F(2, 1);
                    b.set(0, 0, newFunctionNorm - initialFunctionNorm - lambda * initialSlope);
                    b.set(1, 0, previousFunctionNorm - initialFunctionNorm - lambdaPrevious * initialSlope);
                    DenseMatrix64F multiplication = new DenseMatrix64F(2, 1);
                    CommonOps.mult(a, b, multiplication);
                    multiplication.set(0, multiplication.get(0, 0) / (lambda - lambdaPrevious));
                    multiplication.set(1, multiplication.get(1, 0) / (lambda - lambdaPrevious));
                    double disc = multiplication.get(1) * multiplication.get(1) - 3 * multiplication.get(0) * initialSlope;
                    if (multiplication.get(0) == 0.0) {
                        //cubic is quadratic
                        lambdaTemp = -initialSlope / (2 * multiplication.get(1));
                    } else {
                        //legitimate cubic
                        lambdaTemp = (-multiplication.get(1) + Math.sqrt(Math.abs(disc))) / (3 * multiplication.get(0));
                    }
                    if (lambdaTemp > 0.5 * lambda) {
                        lambdaTemp = 0.5 * lambda;
                    }
                }
                lambdaPrevious = lambda;
                previousFunctionNorm = newFunctionNorm;
                if (lambdaTemp <= 0.1 * lambda) {
                    lambda = 0.1 * lambda;
                } else {
                    lambda = lambdaTemp;
                }
            }
        }
        return xPlus;
    }

}
