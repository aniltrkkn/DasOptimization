/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package solvers;

import optimization.functionImplementation.ObjectiveFunction;
import optimization.functionImplementation.Options;
import org.ejml.alg.dense.linsol.chol.LinearSolverChol_B64;
import org.ejml.alg.dense.linsol.qr.LinearSolverQrHouseTran_D64;
import org.ejml.data.DenseMatrix64F;
import org.ejml.ops.CommonOps;
import org.ejml.ops.NormOps;
import solverImplementation.LineSearch;
import solverImplementation.NewtonStep;
import solverImplementation.TrussRegionDoubleDogleg;

/**
 *
 * @author O. Anil Turkkan <turkkan.1@osu.edu> - anilturkkan.com
 */
public class NonlinearEquationSolver implements Solver {

    private final Options solverOptions;
    private final ObjectiveFunction equations;
    //number of consecutive past steps with length maxStep
    private int consecmax;
    private DenseMatrix64F x;
    //f(x)
    private DenseMatrix64F fx;
    //J(x)
    private DenseMatrix64F jacobian;
    //termination status
    public static final int SOLVER_RUNNING = 0;
    public static final int CONVERGED__FUNCTION_TOLERANCE = 1;
    public static final int CONVERGED__STEP_TOLERANCE = 2;
    public static final int FAILED__CANNOT_DECREASE_F = 3;
    public static final int FAILED__ITERATION_LIMIT_REACHED = 4;
    public static final int FAILED__TOO_MANY_MAXSTEP = 5;
    public static final int FAILED__ANOTHER_LOCAL_MINIMUM = 6;
    private int terminationStatus;

    public NonlinearEquationSolver(ObjectiveFunction equations, Options solverOptions) {
        //deep copy the options
        this.solverOptions = new Options(solverOptions);
        this.equations = equations;
        consecmax = 0;
    }

    private DenseMatrix64F getJacobian(DenseMatrix64F x) {
        //check if analytical jacobian
        if (solverOptions.isAnalyticalHessian()) {
            return equations.getH(x);
        } else {
            //finite difference Jacobian
            DenseMatrix64F finiteDifferenceJacobian = new DenseMatrix64F(solverOptions.getN(), solverOptions.getN());
            DenseMatrix64F fx = equations.getF(x);
            DenseMatrix64F dummyFx;
            double sqrDelta = Math.sqrt(solverOptions.getMachineEpsilon());
            for (int j = 0; j < solverOptions.getN(); j++) {
                //calculate column j of J
                //get the sign of x
                double sign = Math.signum(x.get(j));
                if (sign == 0.0) {
                    sign = 1.0;
                }
                //calculate the step size
                double stepSize = sqrDelta * Math.max(Math.abs(x.get(j)), Math.abs(1 / solverOptions.getTypicalX().get(j))) * sign;
                //save the initial value
                double temp = x.get(j);
                //calculate x+stepSize
                x.set(j, x.get(j) + stepSize);
                //to reduce finite difference precision errors
                stepSize = x.get(j) - temp;
                //get the value at x+stepsize
                dummyFx = equations.getF(x);
                for (int i = 0; i < solverOptions.getN(); i++) {
                    finiteDifferenceJacobian.set(i, j, (dummyFx.get(i) - fx.get(i)) / stepSize);
                }
                //resetQuick the value
                x.set(j, temp);
            }
            return finiteDifferenceJacobian;
        }
    }

    private boolean checkInitialGuess(DenseMatrix64F initialGuess) {
        double maximumValue = Double.MIN_VALUE;
        DenseMatrix64F functionValues = this.equations.getF(initialGuess);
        //calculate the devation of function at initial guess from 0
        for (int i = 0; i < solverOptions.getN(); i++) {
            maximumValue = Math.max(maximumValue, Math.abs(functionValues.get(i) * solverOptions.getTypicalX().get(i)));
        }
        return maximumValue <= 0.01 * solverOptions.getFunctionTolerance();
    }

    @Override
    public double functionNorm(DenseMatrix64F x) {
        fx = equations.getF(x);
        return 0.5 * NormOps.fastNormP2(fx) * NormOps.fastNormP2(solverOptions.getTypicalF());
    }

    private double norm(DenseMatrix64F f) {
        return NormOps.fastNormP2(f);
    }

    private void checkConvergence(int iteration, DenseMatrix64F x, DenseMatrix64F xPlus, DenseMatrix64F fx, DenseMatrix64F g, boolean solverStatus, boolean maxStepTaken) {
        double functionTolerance = Double.MIN_VALUE;
        for (int i = 0; i < fx.numRows; i++) {
            functionTolerance = Math.max(functionTolerance, solverOptions.getTypicalF().get(i) * Math.abs(fx.get(i)));
        }
        double lastStepMagnitude = Double.MIN_VALUE;
        for (int i = 0; i < x.numRows; i++) {
            lastStepMagnitude = Math.max(lastStepMagnitude, Math.abs(xPlus.get(i) - x.get(i)) / Math.max(Math.abs(xPlus.get(i)), 1 / solverOptions.getTypicalX().get(i)));
        }
        double localMinimum = Double.MIN_VALUE;
        double functionNorm = norm(fx);
        for (int i = 0; i < g.numRows; i++) {
            localMinimum = Math.max(localMinimum, Math.abs(g.get(i)) * Math.max(Math.abs(xPlus.get(i)), 1 / solverOptions.getTypicalX().get(i)) / (Math.max(functionNorm, solverOptions.getN() / 2)));
        }

        if (!solverStatus) {
            terminationStatus = FAILED__CANNOT_DECREASE_F;
        } else if (functionTolerance <= solverOptions.getFunctionTolerance()) {
            terminationStatus = CONVERGED__FUNCTION_TOLERANCE;
        } else if (lastStepMagnitude <= solverOptions.getStepTolerance()) {
            terminationStatus = CONVERGED__STEP_TOLERANCE;
        } else if (iteration >= solverOptions.getMaxIterations()) {
            terminationStatus = FAILED__ITERATION_LIMIT_REACHED;
        } else if (maxStepTaken) {
            consecmax += 1;
            if (consecmax == 5) {
                terminationStatus = FAILED__TOO_MANY_MAXSTEP;
            }
        } else {
            consecmax = 0;
            if (localMinimum <= solverOptions.getMinTolerance()) {
                terminationStatus = FAILED__ANOTHER_LOCAL_MINIMUM;
            }
        }
    }

    public void solve(DenseMatrix64F initialGuess) {
        /* check initial guess */
        if (checkInitialGuess(new DenseMatrix64F(initialGuess))) {
            //
        } else {

        }
        /*initialize iteration*/
        int iteration = 0;
        terminationStatus = SOLVER_RUNNING;
        //x
        x = initialGuess;
        //f(x)
        fx = equations.getF(x);
        //J(x)
        jacobian = getJacobian(x);
        //g(x)=J^T*F(x)
        DenseMatrix64F g = new DenseMatrix64F(solverOptions.getN(), 1);
        CommonOps.multTransA(jacobian, fx, g);
        CommonOps.elementMult(g, solverOptions.getTypicalF());
        /*iterate until solver succeeds, fails or maximum iteration number is reached*/
        while (terminationStatus == SOLVER_RUNNING) {
            iteration += 1;
            DenseMatrix64F step;
            DenseMatrix64F lowerTriangle;
            if (!solverOptions.isBFGSHessian()) {
                //get the newton step
                step = NewtonStep.newtonStep(jacobian, fx, g, solverOptions);
                lowerTriangle = NewtonStep.getLowerTriangleR();
            } else {
                step = null;
                lowerTriangle = null;
            }
            //get new x values
            DenseMatrix64F xPlus;
            switch (solverOptions.getAlgorithm()) {
                case Options.LINE_SEARCH:
                    xPlus = LineSearch.lineSearch(g, step, x, solverOptions, this);
                    break;
                case Options.DOGLEG_TRUST_REGION:
                    xPlus = TrussRegionDoubleDogleg.dogDriver(g, step, x, lowerTriangle, solverOptions, this);
                    break;
                default:
                    xPlus = LineSearch.lineSearch(g, step, x, solverOptions, this);
            }
            //new function values
            fx = equations.getF(xPlus);
            //get new jacobian
            jacobian = getJacobian(xPlus);
            //get new gradient g(x)=J^T*F(x)
            CommonOps.multTransA(jacobian, fx, g);
            CommonOps.elementMult(g, solverOptions.getTypicalF());
            //check for convergence
            checkConvergence(iteration, x, xPlus, fx, g, LineSearch.getSolverStatus(), LineSearch.isMaxStepTaken());
            //update x
            x = xPlus;
            System.out.println(iteration);
        }
        //System.out.println(x);
        //System.out.println(terminationStatus);
    }

    public DenseMatrix64F getX() {
        return x;
    }

    public DenseMatrix64F getFx() {
        return fx;
    }

    public DenseMatrix64F getJacobian() {
        return jacobian;
    }

}
