/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package solvers;

import cern.colt.matrix.DoubleMatrix1D;
import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.impl.DenseDoubleMatrix1D;
import cern.colt.matrix.impl.DenseDoubleMatrix2D;
import cern.colt.matrix.linalg.Algebra;
import cern.colt.matrix.linalg.CholeskyDecomposition;
import cern.colt.matrix.linalg.QRDecomposition;
import optimization.functionImplementation.ObjectiveFunction;
import optimization.functionImplementation.Options;
import org.ejml.data.DenseMatrix64F;
import org.ejml.ops.NormOps;

/**
 *
 * @author O. Anil Turkkan <turkkan.1@osu.edu> - anilturkkan.com
 */
public class NonlinearEquationSolver {

    private final Options solverOptions;
    private final ObjectiveFunction equations;
    private final Algebra algebra;
    //number of consecutive past steps with length maxStep
    private int consecmax;
    private boolean maxStepTaken;
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

    public NonlinearEquationSolver(ObjectiveFunction equations, Options solverOptions, double[] initialGuess) {
        //deep copy the options
        this.solverOptions = new Options(solverOptions);
        this.equations = equations;
        //check if initial guess input is correct size
        assert solverOptions.getN() == initialGuess.length;
        //initiliaze algebra
        algebra = new Algebra();
        //calculate the maximum step size
        this.solverOptions.setMaxStep(new DenseDoubleMatrix1D(initialGuess));
        consecmax = 0;
        if (checkInitialGuess(new DenseDoubleMatrix1D(initialGuess))) {
            //
        } else {

        }
    }

    private DoubleMatrix2D getJacobian(DoubleMatrix1D x) {
        //check if analytical jacobian
        if (solverOptions.isAnalyticalHessian()) {
            return equations.getH(x);
        } else {
            //finite difference Jacobian
            DoubleMatrix2D finiteDifferenceJacobian = new DenseDoubleMatrix2D(solverOptions.getN(), solverOptions.getN());
            DoubleMatrix1D fx = equations.getF(x);
            DoubleMatrix1D dummyFx;
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
                x.setQuick(j, x.get(j) + stepSize);
                //to reduce finite difference precision errors
                stepSize = x.get(j) - temp;
                //get the value at x+stepsize
                dummyFx = equations.getF(x);
                for (int i = 0; i < solverOptions.getN(); i++) {
                    finiteDifferenceJacobian.setQuick(i, j, (dummyFx.get(i) - fx.get(i)) / stepSize);
                }
                //resetQuick the value
                x.setQuick(j, temp);
            }
            return finiteDifferenceJacobian;
        }
    }

    private boolean checkInitialGuess(DoubleMatrix1D initialGuess) {
        double maximumValue = Double.MIN_VALUE;
        DoubleMatrix1D functionValues = this.equations.getF(initialGuess);
        //calculate the devation of function at initial guess from 0
        for (int i = 0; i < solverOptions.getN(); i++) {
            maximumValue = Math.max(maximumValue, Math.abs(functionValues.get(i) * solverOptions.getTypicalX().get(i)));
        }
        return maximumValue <= 0.01 * solverOptions.getFunctionTolerance();
    }

    private double totalError(DoubleMatrix1D fx) {
        return 0.5 * fx.zDotProduct(fx) * solverOptions.getTypicalF().zDotProduct(solverOptions.getTypicalF());
    }

    private double norm(DoubleMatrix1D f) {
        return Math.sqrt(f.zDotProduct(f));
    }

    private DoubleMatrix2D oneD2twoD(DoubleMatrix1D oneD) {
        DoubleMatrix2D twoD = oneD.like2D(oneD.size(), 2);
        for (int i = 0; i < oneD.size(); i++) {
            twoD.setQuick(i, 0, oneD.get(i));
        }
        return twoD;
    }

    private DoubleMatrix1D twoD2oneD(DoubleMatrix2D twoD) {
        DoubleMatrix1D oneD = twoD.like1D(twoD.rows());
        for (int i = 0; i < oneD.size(); i++) {
            oneD.setQuick(i, twoD.get(i, 0));
        }
        return oneD;
    }

    private DoubleMatrix1D newtonStep(DoubleMatrix2D jacobian, DoubleMatrix1D fx, DoubleMatrix1D g) {
        /*Calculate QR decomposition of DfJ*/
        DoubleMatrix2D dummyJacobian = jacobian.like(jacobian.rows(), jacobian.columns());
        dummyJacobian.assign(jacobian);
        //Df*J
        for (int i = 0; i < dummyJacobian.rows(); i++) {
            for (int j = 0; j < dummyJacobian.columns(); j++) {
                dummyJacobian.setQuick(i, j, dummyJacobian.get(i, j) * solverOptions.getTypicalF().get(i));
            }
        }
        //qr decomposition
        QRDecomposition qrDecomposition = new QRDecomposition(dummyJacobian);
        /*get the condition number*/
        DoubleMatrix2D rMatrix = qrDecomposition.getR();
        //R*Dx^-1
        for (int j = 0; j < rMatrix.rows(); j++) {
            for (int i = 0; i < rMatrix.columns(); i++) {
                rMatrix.setQuick(i, j, rMatrix.get(i, j) / solverOptions.getTypicalX().get(j));
            }
        }
        DenseMatrix64F asd= new DenseMatrix64F(jacobian.toArray());
        double conditionNumber= NormOps.conditionP(asd, 1.0);
        //get the condition number
        //double conditionNumber = algebra.cond(rMatrix);
        //double conditionNumber = ;
        if (qrDecomposition.hasFullRank() && conditionNumber < 1e-3 / Math.sqrt(solverOptions.getMachineEpsilon())) {
            /*calculate the newton step -J*sn=fx*/
            //fix rMatrix
            for (int j = 0; j < rMatrix.rows(); j++) {
                for (int i = 0; i < rMatrix.columns(); i++) {
                    rMatrix.setQuick(i, j, rMatrix.get(i, j) * solverOptions.getTypicalX().get(j));
                }
            }
            DoubleMatrix1D newtonianStep = twoD2oneD(qrDecomposition.solve(oneD2twoD(fx)));
            return newtonianStep.assign(solverOptions.getTypicalF(), (double a, double b) -> -a * b);
        } else {
            //ill-conditioned or singular jacobian
            /* solve -H*sn=g where
                H=J^T*Sf^2*J+sqrt(n*machineEpsilon)*||J^T*Sf^2*J||*Sx^2
            */
            //H=J^T*Sf^2*J
            DoubleMatrix2D h = algebra.mult(algebra.transpose(jacobian), jacobian);
            for(int i=0;i<h.rows();i++){
                for(int j=0;j<h.columns();j++){
                    h.setQuick(i, j, h.get(i, j)*solverOptions.getTypicalF().get(i)*solverOptions.getTypicalF().get(i));
                }
            }
            double normH=algebra.norm1(h);
            for(int i=0;i<h.rows();i++){
                h.setQuick(i, i,h.get(i, i)+Math.sqrt(h.rows()*solverOptions.getMachineEpsilon())*normH*solverOptions.getTypicalX().get(i)*solverOptions.getTypicalX().get(i));
            }
            //cholesky decomposition
            CholeskyDecomposition cDecomposition= new CholeskyDecomposition(h);
            DoubleMatrix1D newtonianStep = twoD2oneD(cDecomposition.solve(oneD2twoD(g)));
            return newtonianStep.assign(solverOptions.getTypicalF(), (double a, double b) -> -a);
        }
    }

    //xp=x+lambda*sn such that f(xp) <= f(x) + alpha*lambda*g^T*p
    private DoubleMatrix1D lineSearch(DoubleMatrix1D g, DoubleMatrix1D sn, DoubleMatrix1D x) {
        //maximum step taken
        maxStepTaken = false;
        //solver status
        solverStatus = STATUS_INSIDE_FUNCTION;
        //alpha
        double alpha = 1e-4;
        //norm of sn*Sx
        double newtonLength = norm(sn.copy().assign(solverOptions.getTypicalX(), (double a, double b) -> a * b));
        if (newtonLength > solverOptions.getMaxStep()) {
            //newton step xp=x+sn is longer than maximum allowed
            for (int i = 0; i < sn.size(); i++) {
                sn.setQuick(i, sn.get(i) * (solverOptions.getMaxStep() / newtonLength));
            }
            newtonLength = solverOptions.getMaxStep();
        }
        //initial slope
        double initialSlope = algebra.mult(g, sn);
        //relative length of sn as calculated in the stopping routine
        double relativeLength = Double.MIN_VALUE;
        for (int i = 0; i < sn.size(); i++) {
            relativeLength = Math.max(relativeLength, Math.abs(sn.get(i)) / Math.max(Math.abs(x.get(i)), 1 / solverOptions.getTypicalX().get(i)));
        }
        //minimum allowable step length
        double minLambda = solverOptions.getStepTolerance() / relativeLength;
        /* find the new x values*/
        DoubleMatrix1D xPlus = new DenseDoubleMatrix1D(x.size());
        double lambda = 1.0;
        double initialFunctionNorm = norm(equations.getF(x));
        double lambdaPrevious = 1.0;
        double previousFunctionNorm = 0.0;
        //loop to check whether xp=x+lambda*sn is satisfactory (generate new lambda if required)
        while (solverStatus == STATUS_INSIDE_FUNCTION) {
            //new x values
            for (int i = 0; i < x.size(); i++) {
                xPlus.setQuick(i, x.get(i) + lambda * sn.get(i));
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
                    DoubleMatrix2D a = new DenseDoubleMatrix2D(new double[][]{{1 / (lambda * lambda), -1 / (lambdaPrevious * lambdaPrevious)}, {-lambdaPrevious / (lambda * lambda), lambda / (lambdaPrevious * lambdaPrevious)}});
                    DoubleMatrix1D b = new DenseDoubleMatrix1D(new double[]{newFunctionNorm - initialFunctionNorm - lambda * initialSlope, previousFunctionNorm - initialFunctionNorm - lambdaPrevious * initialSlope});
                    DoubleMatrix1D multiplication = algebra.mult(a, b);
                    multiplication.setQuick(0, multiplication.get(0) / (lambda - lambdaPrevious));
                    multiplication.setQuick(1, multiplication.get(1) / (lambda - lambdaPrevious));
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

    private void checkConvergence(int iteration, DoubleMatrix1D x, DoubleMatrix1D xPlus, DoubleMatrix1D fx, DoubleMatrix1D g) {
        double functionTolerance= Double.MIN_VALUE;
        for(int i=0;i<fx.size();i++){
            functionTolerance=Math.max(functionTolerance,solverOptions.getTypicalF().get(i)*Math.abs(fx.get(i)));
        }
        double lastStepMagnitude= Double.MIN_VALUE;
        for(int i=0;i<x.size();i++){
            lastStepMagnitude=Math.max(lastStepMagnitude,Math.abs(xPlus.get(i)-x.get(i))/Math.max(Math.abs(xPlus.get(i)),1/solverOptions.getTypicalX().get(i) ));
        }
        double localMinimum= Double.MIN_VALUE;
        double functionNorm= norm(fx);
        for(int i=0;i<g.size();i++){
            localMinimum=Math.max(localMinimum,Math.abs(g.get(i))*Math.max(Math.abs(xPlus.get(i)),1/solverOptions.getTypicalX().get(i) )/(Math.max(functionNorm,solverOptions.getN()/2)) );
        }
        
        if(solverStatus== STATUS_FAILED){
            terminationStatus = FAILED__CANNOT_DECREASE_F;
        } else if (functionTolerance<=solverOptions.getFunctionTolerance()){
            terminationStatus =CONVERGED__FUNCTION_TOLERANCE;
        } else if(lastStepMagnitude<= solverOptions.getStepTolerance()){
            terminationStatus=CONVERGED__STEP_TOLERANCE;
        } else if (iteration>= solverOptions.getMaxIterations()){
            terminationStatus=FAILED__ITERATION_LIMIT_REACHED;
        } else if (maxStepTaken){
            consecmax += 1;
            if (consecmax == 5){
                terminationStatus=FAILED__TOO_MANY_MAXSTEP;
            }
        } else {
            consecmax=0;
            if(localMinimum<=solverOptions.getMinTolerance()){
                terminationStatus=FAILED__ANOTHER_LOCAL_MINIMUM;
            }
        }
    }

    public void solve(DoubleMatrix1D initialGuess) {
        /*initialize iteration*/
        int iteration = 0;
        terminationStatus = SOLVER_RUNNING;
        //x
        DoubleMatrix1D x = initialGuess;
        //f(x)
        DoubleMatrix1D fx = equations.getF(x);
        //J(x)
        DoubleMatrix2D jacobian = getJacobian(x);
        //g(x)=J^T*F(x)
        DoubleMatrix1D g = algebra.mult(algebra.transpose(jacobian), fx);
        g.assign(solverOptions.getTypicalF(), (double a, double b) -> a * b * b);
        /*iterate until solver succeeds, fails or maximum iteration number is reached*/
        while (terminationStatus == SOLVER_RUNNING) {
            iteration += 1;
            //get the newton step
            DoubleMatrix1D newtonianStep = newtonStep(jacobian, fx,g);
            //get new x values
            DoubleMatrix1D xPlus;
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
            g = algebra.mult(algebra.transpose(jacobian), fx);
            g.assign(solverOptions.getTypicalF(), (double a, double b) -> a * b * b);
            //check for convergence
            checkConvergence(iteration, x, xPlus, fx,g);
            //update x
            x = xPlus;
            System.out.println(iteration);
        }
        System.out.println(x);
        //System.out.println(terminationStatus);
    }

}
