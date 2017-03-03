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
import org.ejml.alg.dense.mult.MatrixVectorMult;
import org.ejml.data.DenseMatrix64F;
import org.ejml.ops.CommonOps;
import org.ejml.ops.MatrixFeatures;
import org.ejml.ops.NormOps;

/**
 *
 * @author O. Anil Turkkan <turkkan.1@osu.edu> - anilturkkan.com
 */
public class NonlinearEquationSolver {

    private final Options solverOptions;
    private final ObjectiveFunction equations;
    //number of consecutive past steps with length maxStep
    private int consecmax;
    private boolean maxStepTaken;
    private DenseMatrix64F x;
    //f(x)
    private DenseMatrix64F fx;
    //J(x)
    private DenseMatrix64F jacobian;
    //solver status
    public static final int STATUS_ONGOING = 0;
    public static final int STATUS_FAILED = 1;
    public static final int STATUS_INSIDE_FUNCTION = 2;
    private int solverStatus;
    //termination status
    public static final int SOLVER_RUNNING = 0;
    public static final int CONVERGED__FUNCTION_TOLERANCE = 1;
    public static final int CONVERGED__STEP_TOLERANCE = 2;
    public static final int FAILED__CANNOT_DECREASE_F = 3;
    public static final int FAILED__ITERATION_LIMIT_REACHED = 4;
    public static final int FAILED__TOO_MANY_MAXSTEP = 5;
    public static final int FAILED__ANOTHER_LOCAL_MINIMUM = 6;
    private int terminationStatus;
    //large or small matrix distinction
    public static final int LARGE_MATRIX = 100000;

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

    private double totalError(DenseMatrix64F fx) {
        return 0.5 * NormOps.fastNormP2(fx) * NormOps.fastNormP2(solverOptions.getTypicalF());
    }

    private double norm(DenseMatrix64F f) {
        return NormOps.fastNormP2(f);
    }

    private DenseMatrix64F newtonStep(DenseMatrix64F jacobian, DenseMatrix64F fx, DenseMatrix64F g) {
        /*Calculate QR decomposition of DfJ*/
        DenseMatrix64F dummyJacobian = jacobian.copy();
        //Df*J
        for (int i = 0; i < dummyJacobian.numRows; i++) {
            for (int j = 0; j < dummyJacobian.numCols; j++) {
                dummyJacobian.set(i, j, dummyJacobian.get(i, j) * solverOptions.getTypicalF().get(i));
            }
        }
        //qr decomposition
        LinearSolverQrHouseTran_D64 qrSolver = new LinearSolverQrHouseTran_D64();
        
        qrSolver.setA(dummyJacobian);
        /*get the condition number*/
        double conditionNumber = qrSolver.quality();
        //get the condition number
        //double conditionNumber = algebra.cond(rMatrix);
        //double conditionNumber = ;
        //if (MatrixFeatures.rank(dummyJacobian) == solverOptions.getN()&& conditionNumber > 1e3*Math.sqrt(solverOptions.getMachineEpsilon())) {
        if (conditionNumber < 1e3*Math.sqrt(solverOptions.getMachineEpsilon())) {
            /*calculate the newton step -J*sn=fx*/
            DenseMatrix64F newtonianStep = new DenseMatrix64F(dummyJacobian.numRows, 1);
            qrSolver.solve(fx, newtonianStep);
            CommonOps.elementMult(newtonianStep, solverOptions.getTypicalF());
            CommonOps.changeSign(newtonianStep);
            return newtonianStep;
        } else {
            //ill-conditioned or singular jacobian
            /* solve -H*sn=g where
                H=J^T*Sf^2*J+sqrt(n*machineEpsilon)*||J^T*Sf^2*J||*Sx^2
             */
            //H=J^T*Sf^2*J
            DenseMatrix64F h = new DenseMatrix64F(dummyJacobian.numRows, 1);
            CommonOps.multInner(jacobian, jacobian);
            for (int i = 0; i < h.numRows; i++) {
                for (int j = 0; j < h.numCols; j++) {
                    h.set(i, j, h.get(i, j) * solverOptions.getTypicalF().get(i) * solverOptions.getTypicalF().get(i));
                }
            }
            double normH = NormOps.normP1(h);
            for (int i = 0; i < h.numRows; i++) {
                h.set(i, i, h.get(i, i) + Math.sqrt(h.numRows * solverOptions.getMachineEpsilon()) * normH * solverOptions.getTypicalX().get(i) * solverOptions.getTypicalX().get(i));
            }
            //cholesky decomposition
            LinearSolverChol_B64 cDecomposition = new LinearSolverChol_B64();
            cDecomposition.setA(h);
            DenseMatrix64F newtonianStep = new DenseMatrix64F(dummyJacobian.numRows, 1);
            cDecomposition.solve(g, newtonianStep);
            CommonOps.changeSign(newtonianStep);
            return newtonianStep;
        }
    }

    //xp=x+lambda*sn such that f(xp) <= f(x) + alpha*lambda*g^T*p
    private DenseMatrix64F lineSearch(DenseMatrix64F g, DenseMatrix64F sn, DenseMatrix64F x) {
        //maximum step taken
        maxStepTaken = false;
        //solver status
        solverStatus = STATUS_INSIDE_FUNCTION;
        //alpha
        double alpha = 1e-4;
        //norm of sn*Sx
        DenseMatrix64F dummySn = new DenseMatrix64F(sn);
        CommonOps.elementMult(dummySn, solverOptions.getTypicalX());
        double newtonLength = NormOps.conditionP2(dummySn);
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
        DenseMatrix64F xPlus = new DenseMatrix64F(x.numRows,x.numCols);
        double lambda = 1.0;
        double initialFunctionNorm = norm(equations.getF(x));
        double lambdaPrevious = 1.0;
        double previousFunctionNorm = 0.0;
        //loop to check whether xp=x+lambda*sn is satisfactory (generate new lambda if required)
        while (solverStatus == STATUS_INSIDE_FUNCTION) {
            //new x values
            for (int i = 0; i < x.numRows; i++) {
                xPlus.set(i, 0, x.get(i,0) + lambda * sn.get(i,0));
            }
            //new function norm
            double newFunctionNorm = norm(equations.getF(xPlus));
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
                        lambdaTemp = (-multiplication.get(1) + Math.sqrt(disc)) / (3 * multiplication.get(0));
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

    private void checkConvergence(int iteration, DenseMatrix64F x, DenseMatrix64F xPlus, DenseMatrix64F fx, DenseMatrix64F g) {
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

        if (solverStatus == STATUS_FAILED) {
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
        CommonOps.multTransA(jacobian, fx, g);
        CommonOps.elementMult(g, solverOptions.getTypicalF());
        /*iterate until solver succeeds, fails or maximum iteration number is reached*/
        while (terminationStatus == SOLVER_RUNNING) {
            iteration += 1;
            //get the newton step
            DenseMatrix64F newtonianStep = newtonStep(jacobian, fx, g);
            //get new x values
            DenseMatrix64F xPlus;
            switch (solverOptions.getAlgorithm()) {
                case Options.LINE_SEARCH:
                    xPlus = lineSearch(g, newtonianStep, x);
                    break;
                default:
                    xPlus = lineSearch(g, newtonianStep, x);
            }
            //new function values
            fx = equations.getF(xPlus);
            //get new jacobian
            jacobian = getJacobian(xPlus);
            //get new gradient g(x)=J^T*F(x)
            CommonOps.multTransA(jacobian, fx, g);
            CommonOps.elementMult(g, solverOptions.getTypicalF());
            //check for convergence
            checkConvergence(iteration, x, xPlus, fx, g);
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
