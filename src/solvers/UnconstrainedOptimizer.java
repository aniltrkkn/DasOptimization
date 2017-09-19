/*
 * Copyright (c) 2016-2017, O. Anil Turkkan. All Rights Reserved.
 *
 * This file is part of DASOptimization.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *   http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
package solvers;

import descentAlgorithms.LineSearch;
import descentAlgorithms.StepAlgorithm;
import descentAlgorithms.TrussRegionDoubleDogleg;
import finiteDifferenceApproximations.FiniteDifference;
import optimization.functionImplementation.Options;
import optimization.functionImplementation.Results;
import org.ejml.data.DMatrixRMaj;
import optimization.functionImplementation.ObjectiveFunctionUnconstrained;


public class UnconstrainedOptimizer implements Solver {

    private final Options solverOptions;
    private final ObjectiveFunctionUnconstrained equations;
    private final Results results;
    //number of consecutive past steps with length maxStep
    private int consecmax;
    private DMatrixRMaj x;
    private double fx; //f(x)
    private DMatrixRMaj gx; //g(x)
    private DMatrixRMaj hx; //H(x)
    //termination status
    public static final int SOLVER_RUNNING = 0;
    public static final int CONVERGED__FUNCTION_TOLERANCE = 1;
    public static final int CONVERGED__GRADIENT_TOLERANCE = 2;
    public static final int CONVERGED__STEP_TOLERANCE = 3;
    public static final int FAILED__CANNOT_DECREASE_F = 4;
    public static final int FAILED__ITERATION_LIMIT_REACHED = 5;
    public static final int FAILED__TOO_MANY_MAXSTEP = 6;
    public static final int FAILED__ANOTHER_LOCAL_MINIMUM = 7;
    private int terminationStatus;

    public UnconstrainedOptimizer(ObjectiveFunctionUnconstrained equations, Options solverOptions) {
        //deep copy the options
        this.solverOptions = solverOptions;
        this.equations = equations;
        this.consecmax = 0;
        this.results = new Results();
    }

    /**
     * solve (LL^T)step= -g(x) for step
     * @param step descent step
     * @param lowerTriangle Cholesky decomposition of the Hessian
    */
    public void solveCholesky(DMatrixRMaj step,DMatrixRMaj lowerTriangle){
        int n=step.numRows;
        /* Solve Ly=step */
        DMatrixRMaj y= new DMatrixRMaj(lowerTriangle.numCols, 1);
        y.set(0, 0, gx.get(0)/lowerTriangle.get(0,0));
        for(int i=1;i<n;i++){
            y.set(i,0,gx.get(i));
            for(int j=0;j<i;j++){
                y.set(i,0,y.get(i,0)-lowerTriangle.get(i,j)*y.get(j,0));
            }
            y.set(i,0,y.get(i,0)/lowerTriangle.get(i,i));
        }
        /* solve L^Ts = y */
        step.set(n-1, 0, y.get(n-1)/lowerTriangle.get(n-1,n-1));
        for(int i=n-2;i>=0;i--){
            step.set(i,0,y.get(i));
            for(int j=i+1;j<n;j++){
                step.set(i,0,step.get(i,0)-lowerTriangle.get(j,i)*step.get(j,0));
            }
            step.set(i,0,step.get(i,0)/lowerTriangle.get(i,i));
        }
        /*inverse the step sign */
        for(int i=0;i<n;i++){
            step.set(i,0,-step.get(i,0));
        }
    }
    
    /**
     * Perform a Cholesky decomposition and store it in the lower triangle
     * @param lowerTriangle lower triangle Cholesky decomposition
     * @param maximumMember maximum member of the Hessian
     * @return 
     */
    private double perturbedCholeskyDecomposition(DMatrixRMaj lowerTriangle,double maximumMember) {
        double minimumValue= Math.pow(this.solverOptions.getMachineEpsilon(),1.0/4.0)* maximumMember;
        if(maximumMember <1e-10){
            //this means H is positive definite
            for(int i=0;i<hx.numCols;i++){
                maximumMember=Math.max(maximumMember, Math.abs(hx.get(i,i)));
            }
            maximumMember = Math.sqrt(maximumMember);
        }
        double minimumValue2 = Math.sqrt(this.solverOptions.getMachineEpsilon())*maximumMember;
        double addition=0.0; //this value contains the maximum amount that will be added to H
        /* main loop to form column of lower triangle */
        for(int j=0;j<hx.numCols;j++){
            lowerTriangle.set(j,j,hx.get(j,j));
            for(int i=0;i<j;i++){
                lowerTriangle.set(j,j,lowerTriangle.get(j,j)-lowerTriangle.get(j,i)*lowerTriangle.get(j,i));
            }
            double minimumJJ =0.0;
            for(int i=j+1;i<hx.numCols;i++){
                lowerTriangle.set(i,j,hx.get(j,i));
                for(int k=0;k<j;k++){
                    lowerTriangle.set(i,j,lowerTriangle.get(i,j)-lowerTriangle.get(i,k)*lowerTriangle.get(j,k));
                }
                minimumJJ=Math.max(minimumJJ,Math.abs(lowerTriangle.get(i,j)));
            }
            minimumJJ = Math.max(minimumJJ/maximumMember, minimumValue);
            if(lowerTriangle.get(j,j) > minimumJJ*minimumJJ){
                //normal Cholesky iteration
                lowerTriangle.set(j,j,Math.sqrt(lowerTriangle.get(j,j)));
            } else {
                //augment H(j,j)
                if(minimumJJ <minimumValue2 ){
                    minimumJJ=minimumValue2;
                }
                addition=Math.max(addition,minimumJJ*minimumJJ - lowerTriangle.get(j,j));
                lowerTriangle.set(j,j,minimumJJ);
            }
            for(int i=j+1;i<hx.numCols;i++){
                lowerTriangle.set(i,j,lowerTriangle.get(i,j)/lowerTriangle.get(j,j));
            }
        }
        return addition;
    }

    /**
     * convert the Hessian to a positive definite matrix and solve using
     * Cholesky decomposition
     */
    private void perturbedCholeskyDecompositionDriver(DMatrixRMaj lowerTriangle) {
        double epsilon = Math.sqrt(this.solverOptions.getMachineEpsilon());
        /* find the maximum and minimum diagonal */
        double minimumDiagonal = Double.MAX_VALUE;
        double maximumDiagonal = -Double.MAX_VALUE;
        for (int i = 0; i < hx.numCols; i++) {
            minimumDiagonal = Math.min(minimumDiagonal, hx.get(i, i));
            maximumDiagonal = Math.max(maximumDiagonal, hx.get(i, i));
        }
        double maxPosDiagonal = Math.max(0.0, maximumDiagonal);
        /* find the amount to add to diagonal */
        double addition = 0;
        if (minimumDiagonal <= epsilon * maxPosDiagonal) {
            addition = 2 * (maxPosDiagonal - minimumDiagonal) * epsilon - minimumDiagonal;
            maximumDiagonal += addition;
        }
        //calculate the maximum offdiagonal
        double maxOffDiagonal = 0.0;
        for (int i = 0; i < hx.numCols; i++) {
            for (int j = i + 1; j < hx.numCols; j++) {
                maxOffDiagonal = Math.max(maxOffDiagonal, Math.abs(hx.get(i, j)));
            }
        }
        if (maxOffDiagonal * (1 + 2 * epsilon) > maximumDiagonal) {
            addition += (maxOffDiagonal - maximumDiagonal) + 2 * epsilon * maxOffDiagonal;
            maximumDiagonal = maxOffDiagonal * (1 + 2 * epsilon);
        }
        if (maximumDiagonal == 0.0) { //H=0 
            addition = 1.0;
            maximumDiagonal = 1.0;
        }
        /* add the additional value */
        for (int i = 0; i < hx.numCols; i++) {
            hx.set(i, i, hx.get(i, i) + addition);
        }
        /* call Cholesky decomposition */
        addition = this.perturbedCholeskyDecomposition(lowerTriangle,Math.sqrt(Math.max(maximumDiagonal, maxOffDiagonal / hx.numCols)));
        /* check if it is successful */
        if (addition > 1e-7) {
            /* still not positive definite */
            double minValue = this.hx.get(0, 0);
            double maxValue = this.hx.get(0, 0);
            for (int i = 0; i < hx.numRows; i++) {
                double offRow = 0.0;
                for (int j = 0; j < i ;j++) {
                    offRow += Math.abs(hx.get(j, i));
                }
                for (int j = i + 1; j < hx.numRows; j++) {
                    offRow += Math.abs(hx.get(i, j));
                }
                minValue = Math.min(minValue, hx.get(i, i) - offRow);
                maxValue = Math.max(maxValue, hx.get(i, i) + offRow);
            }
            double addition2 = (maxValue - minValue) * epsilon - minValue;
            addition=Math.min(Math.max(0,addition2),addition);
            //add the value to hessian and make it positive definite
            for (int i = 0; i < hx.numRows; i++) {
                hx.set(i, i, hx.get(i, i) + addition);
            }
            
        }
        this.perturbedCholeskyDecomposition(lowerTriangle,0.0);
    }

    /**
     * The function value for given x
     *
     * @param x function variables where norm is calculated
     * @return
     */
    @Override
    public double functionNorm(DMatrixRMaj x) {
        if (this.solverOptions.isSaveIterationDetails()) {
            this.results.updateFunctionEvaluations();
        }
        return equations.getF(x);
    }

    /**
     * Check the unconstrained equation at the initial point return true if the
     * gradient is within the gradient tolerance
     *
     * @param initialGuess initial guess supplied by the user
     * @return true if gradient is zero at the initial guess
     */
    private boolean checkInitialGuess() {
        double maximumValue = Double.MIN_VALUE;
        //calculate the devation of function at initial guess from 0
        for (int i = 0; i < solverOptions.getN(); i++) {
            maximumValue = Math.max(maximumValue, Math.abs(gx.get(i)) * Math.max(Math.abs(x.get(i)), 1.0) / Math.max(Math.abs(fx), 1.0));
        }
        return maximumValue <= 0.001 * solverOptions.getGradientTolerance();
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
    private void checkConvergence(int iteration, DMatrixRMaj x, DMatrixRMaj xPlus, double fx, DMatrixRMaj g, boolean solverStatus, boolean maxStepTaken) {
        /*
        calculate the maximum component of the scaled function
         */
        double functionTolerance = Double.MIN_VALUE;
        for (int i = 0; i < x.numRows; i++) {
            functionTolerance = Math.max(functionTolerance, Math.abs(fx));
        }
        /*
        maximum value of the scaled step
         */
        double lastStepMagnitude = Double.MIN_VALUE;
        for (int i = 0; i < x.numRows; i++) {
            lastStepMagnitude = Math.max(lastStepMagnitude, Math.abs(xPlus.get(i, 0) - x.get(i, 0)) / Math.max(Math.abs(xPlus.get(i)), 1.0));
        }
        /*
        check if gradient is zero
         */
        double gradientMagnitude = Double.MIN_VALUE;
        //calculate the devation of function at initial guess from 0
        for (int i = 0; i < solverOptions.getN(); i++) {
            gradientMagnitude = Math.max(gradientMagnitude, Math.abs(gx.get(i)) * Math.abs(Math.max(Math.abs(x.get(i)), 1.0) / Math.max(Math.abs(fx), 1.0)));
        }

        if (!solverStatus) {
            terminationStatus = FAILED__CANNOT_DECREASE_F;
        } else if (functionTolerance <= solverOptions.getFunctionTolerance()) {
            terminationStatus = CONVERGED__FUNCTION_TOLERANCE;
        } else if (lastStepMagnitude <= solverOptions.getStepTolerance()) {
            terminationStatus = CONVERGED__STEP_TOLERANCE;
        } else if (gradientMagnitude <= solverOptions.getMinTolerance()) {
            terminationStatus = CONVERGED__GRADIENT_TOLERANCE;
        } else if (iteration >= this.solverOptions.getMaxIterations()) {
            terminationStatus = FAILED__ITERATION_LIMIT_REACHED;
        } else if (maxStepTaken) {
            consecmax += 1;
            if (consecmax == 5) {
                terminationStatus = FAILED__TOO_MANY_MAXSTEP;
            }
        } else {
            consecmax = 0;
        }
    }

    /**
     * Run the solver with the options specified with the given options
     *
     * @param initialGuess initial guess for the solver
     */
    public void solve(DMatrixRMaj initialGuess) {
        StepAlgorithm descentAlgorithm;
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
        /*initialize iteration*/
        int iteration = 0;
        terminationStatus = SOLVER_RUNNING;
        //x
        x = initialGuess;
        //calculate objective function
        fx = equations.getF(x);
        //calculate gradient
        gx = FiniteDifference.getGradient(x, equations, solverOptions);
        //check initial guess
        if (checkInitialGuess()) {
            return;
        }
        //calculate hessian
        hx = FiniteDifference.getHessian(x, equations, solverOptions);
        /*iterate until solver succeeds, fails or maximum iteration number is reached*/
        DMatrixRMaj lowerTriangle = new DMatrixRMaj(hx.numRows, hx.numRows);
        DMatrixRMaj step = new DMatrixRMaj(hx.numRows, 1);
        //LinearSolver<DMatrixRMaj> cDecomposition = LinearSolverFactory_DDRM.chol(x.numRows);
        while (terminationStatus == SOLVER_RUNNING) {
            iteration++;
            /*Cholesky Decomposition*/
            this.perturbedCholeskyDecompositionDriver(lowerTriangle);
            this.solveCholesky(step, lowerTriangle);
            /*test */
            /*cDecomposition.setA(hx);
            cDecomposition.solve(gx, step);
            CommonOps_DDRM.changeSign(step);
            ((CholeskyDecomposition) cDecomposition.getDecomposition()).getT(lowerTriangle);
            if (!((CholeskyDecomposition) cDecomposition.getDecomposition()).isLower()) {
                CommonOps_DDRM.transpose(lowerTriangle);
            }*/
            /*get new x values */
            DMatrixRMaj xNew = descentAlgorithm.solve(gx, step, x, lowerTriangle, solverOptions, this);
            /*get new gradient(x)*/
            DMatrixRMaj gxNew = FiniteDifference.getGradient(xNew, equations, solverOptions);
            /* get new function magnitude */
            double fxNew = functionNorm(xNew);
            /*check convergence */
            checkConvergence(iteration, x, xNew, fxNew, gxNew, descentAlgorithm.isSolverFailed(), descentAlgorithm.isMaxStepTaken());
            /*update everythin */
            x = xNew;
            gx = gxNew;
            fx = fxNew;
            hx = FiniteDifference.getHessian(x, equations, solverOptions);
            //update results
            if (this.solverOptions.isSaveIterationDetails()) {
                this.results.update(this.functionNorm(x), gx.data, x.data, this.solverOptions.getTrussRegionRadius());
            }
        }
        /* update if solver is successful or not */
        if (this.solverOptions.isSaveIterationDetails()) {
            this.results.setSuccessful(this.terminationStatus <= UnconstrainedOptimizer.CONVERGED__STEP_TOLERANCE);
        }
    }

    public DMatrixRMaj getX() {
        return x;
    }

    public double getFx() {
        return fx;
    }

    public DMatrixRMaj getGx() {
        return gx;
    }

    public DMatrixRMaj getHx() {
        return hx;
    }

    public int getTerminationStatus() {
        return terminationStatus;
    }

    public Results getResults() {
        return results;
    }
}
