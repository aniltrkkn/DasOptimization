/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package solvers;

import optimization.functionImplementation.Options;
import optimization.functionImplementation.Results;
import org.ejml.data.DMatrixRMaj;
import org.ejml.dense.row.CommonOps_DDRM;
import descentAlgorithms.LineSearch;
import descentAlgorithms.NewtonStep;
import descentAlgorithms.StepAlgorithm;
import descentAlgorithms.TrussRegionDoubleDogleg;
import finiteDifferenceApproximations.FiniteDifference;
import optimization.functionImplementation.ObjectiveFunctionNonLinear;

/**
 *
 * @author O. Anil Turkkan <turkkan.1@osu.edu> - anilturkkan.com
 */
public class NonlinearEquationSolver implements Solver {

    private final Options solverOptions;
    private final ObjectiveFunctionNonLinear equations;
    private final Results results;
    //number of consecutive past steps with length maxStep
    private int consecmax;
    private DMatrixRMaj x;
    //f(x)
    private DMatrixRMaj fx;
    //J(x)
    private DMatrixRMaj jacobian;
    //termination status
    public static final int SOLVER_RUNNING = 0;
    public static final int CONVERGED__FUNCTION_TOLERANCE = 1;
    public static final int CONVERGED__STEP_TOLERANCE = 2;
    public static final int FAILED__CANNOT_DECREASE_F = 3;
    public static final int FAILED__ITERATION_LIMIT_REACHED = 4;
    public static final int FAILED__TOO_MANY_MAXSTEP = 5;
    public static final int FAILED__ANOTHER_LOCAL_MINIMUM = 6;
    private int terminationStatus;

    public NonlinearEquationSolver(ObjectiveFunctionNonLinear equations, Options solverOptions) {
        //deep copy the options
        this.solverOptions = solverOptions;
        this.equations = equations;
        this.consecmax = 0;
        this.results = new Results();
    }

    /**
     * Check the nonlinear equations at the initial point return true if the
     * maximum error is within the function tolerance
     *
     * @param initialGuess initial guess supplied by the user
     * @return true if function is zero at the initial guess
     */
    private boolean checkInitialGuess(DMatrixRMaj initialGuess) {
        double maximumValue = Double.MIN_VALUE;
        DMatrixRMaj functionValues = this.equations.getF(initialGuess);
        //calculate the devation of function at initial guess from 0
        for (int i = 0; i < solverOptions.getN(); i++) {
            maximumValue = Math.max(maximumValue, Math.abs(functionValues.get(i)));
        }
        return maximumValue <= 0.01 * solverOptions.getFunctionTolerance();
    }

    /**
     * The error norm (magnitude) for a given x
     *
     * @param x function variables where norm is calculated
     * @return the magnitude of the error
     */
    @Override
    public double functionNorm(DMatrixRMaj x) {
        if (this.solverOptions.isSaveIterationDetails()) {
            this.results.updateFunctionEvaluations();
        }
        DMatrixRMaj fxDummy = equations.getF(x);
        double functionValue = 0;
        for (int i = 0; i < fxDummy.numRows; i++) {
            functionValue += Math.pow(fxDummy.get(i, 0), 2);
        }
        return 0.5 * functionValue;
    }

    /**
     * Check the convergence of the current iteration Decide if the solver
     * fails, stalls or succeeds
     *
     * @param iteration current iteration number
     * @param x previous x
     * @param xPlus current x
     * @param fx`function values at the current x
     * @param g current function gradient
     * @param solverStatus status of the solver (linear search or truss-region)
     * @param maxStepTaken if last step was equal in magnitude to largest
     * allowed
     */
    private void checkConvergence(int iteration, DMatrixRMaj x, DMatrixRMaj xPlus, DMatrixRMaj fx, DMatrixRMaj g, boolean solverStatus, boolean maxStepTaken) {
        /*
        calculate the maximum component of the scaled function
         */
        double functionTolerance = Double.MIN_VALUE;
        for (int i = 0; i < fx.numRows; i++) {
            functionTolerance = Math.max(functionTolerance,  Math.abs(fx.get(i, 0)));
        }
        /*
        maximum value of the scaled step
         */
        double lastStepMagnitude = Double.MIN_VALUE;
        for (int i = 0; i < x.numRows; i++) {
            lastStepMagnitude = Math.max(lastStepMagnitude, Math.abs(xPlus.get(i, 0) - x.get(i, 0)) /  Math.max(Math.abs(xPlus.get(i)),1.0));
        }
        /*
        check if stuck in a local minimum
         */
        double localMinimum = Double.MIN_VALUE;
        double functionNorm = this.functionNorm(fx);
        for (int i = 0; i < g.numRows; i++) {
            localMinimum = Math.max(localMinimum, Math.abs(g.get(i, 0)) * Math.abs(xPlus.get(i, 0)) / (Math.max(functionNorm, solverOptions.getN() / 2)));
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
            /*if (localMinimum <= solverOptions.getMinTolerance()) {
                terminationStatus = FAILED__ANOTHER_LOCAL_MINIMUM;
            }*/
        }
    }

    /**
     * Run the solver with the options specified with the given options
     *
     * @param initialGuess initial guess for the solver
     */
    public void solve(DMatrixRMaj initialGuess) {
        StepAlgorithm descentAlgorithm;
        NewtonStep newtonStep = new NewtonStep();
        /* check initial guess */
        if (checkInitialGuess(new DMatrixRMaj(initialGuess))) {
            x = initialGuess;
            terminationStatus = NonlinearEquationSolver.CONVERGED__FUNCTION_TOLERANCE;
            return;
        } else {
            /* initiliaze solvers */
            switch (solverOptions.getAlgorithm()) {
                case Options.LINE_SEARCH:
                    descentAlgorithm = new LineSearch();
                    break;
                case Options.DOGLEG_TRUST_REGION:
                    descentAlgorithm = new TrussRegionDoubleDogleg();
                    break;
                default:
                    descentAlgorithm = new LineSearch();
            }
        }
        /*initialize iteration*/
        int iteration = 0;
        terminationStatus = SOLVER_RUNNING;
        //x
        x = initialGuess;
        //f(x)
        fx = equations.getF(x);
        //J(x)
        jacobian = FiniteDifference.getJacobian(x, equations, solverOptions);
        //g(x)=J^T*F(x)
        DMatrixRMaj g = new DMatrixRMaj(solverOptions.getN(), 1);
        CommonOps_DDRM.multTransA(jacobian, fx, g);
        /*iterate until solver succeeds, fails or maximum iteration number is reached*/
        DMatrixRMaj step;
        DMatrixRMaj lowerTriangle;
        while (terminationStatus == SOLVER_RUNNING) {
            iteration++;
            step = newtonStep.newtonStep(jacobian, fx, g, solverOptions);
            lowerTriangle = newtonStep.getLowerTriangleR();
            //get new x values
            DMatrixRMaj xPlus = descentAlgorithm.solve(g, step, x, lowerTriangle, solverOptions, this);
            //new function values
            DMatrixRMaj fxPlus = equations.getF(xPlus);
            DMatrixRMaj gPlus = new DMatrixRMaj(g.numRows, 1);
            //get new jacobian
            jacobian = FiniteDifference.getJacobian(xPlus, equations, solverOptions);
            //get new gradient g(x)=J^T*F(x)
            CommonOps_DDRM.multTransA(jacobian, fxPlus, gPlus);
            //check for convergence
            checkConvergence(iteration, x, xPlus, fxPlus, gPlus, descentAlgorithm.isSolverFailed(), descentAlgorithm.isMaxStepTaken());
            //update x and g
            x = xPlus;
            g = gPlus;
            fx = fxPlus;
            //update results
            if (this.solverOptions.isSaveIterationDetails()) {
                this.results.update(this.functionNorm(x), g.data, x.data, this.solverOptions.getTrussRegionRadius());
            }
            //System.out.println(iteration);
            //System.out.println(x);
            //System.out.println(fx);
        }
        /* update if solver is successful or not */
        if (this.solverOptions.isSaveIterationDetails()) {
            this.results.setSuccessful(this.terminationStatus < NonlinearEquationSolver.FAILED__ANOTHER_LOCAL_MINIMUM);
        }
    }

    public DMatrixRMaj getX() {
        return x;
    }

    public DMatrixRMaj getFx() {
        return fx;
    }

    public DMatrixRMaj getJacobian() {
        return jacobian;
    }

    public int getTerminationStatus() {
        return terminationStatus;
    }

    public Results getResults() {
        return results;
    }

}
