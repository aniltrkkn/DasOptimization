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
package descentAlgorithms;

import optimization.functionImplementation.Options;
import org.ejml.data.DMatrixRMaj;
import org.ejml.dense.row.CommonOps_DDRM;
import org.ejml.dense.row.NormOps_DDRM;
import solvers.Solver;

public class LineSearch implements StepAlgorithm {

    //solver status
    private static final int STATUS_ONGOING = 0;
    private static final int STATUS_FAILED = 1;
    private static final int STATUS_INSIDE_FUNCTION = 2;
    private int solverStatus;
    private boolean maxStepTaken;

    @Override
    public boolean isSolverFailed() {
        return solverStatus != STATUS_FAILED;
    }

    @Override
    public boolean isMaxStepTaken() {
        return maxStepTaken;
    }

    //xp=x+lambda*sn such that f(xp) <= f(x) + alpha*lambda*g^T*p
    @Override
    public DMatrixRMaj solve(DMatrixRMaj g, DMatrixRMaj sn, DMatrixRMaj x, DMatrixRMaj lowerTriangleR, Options solverOptions, Solver solver) {
        //maximum step taken
        maxStepTaken = false;
        //solver status
        solverStatus = STATUS_INSIDE_FUNCTION;
        //alpha
        double alpha = 1e-4;
        //norm of sn*Sx
        DMatrixRMaj dummySn = new DMatrixRMaj(sn);
        double newtonLength = NormOps_DDRM.normP2(dummySn);
        if (newtonLength > solverOptions.getMaxStep()) {
            //newton step xp=x+sn is longer than maximum allowed
            for (int i = 0; i < sn.numRows; i++) {
                sn.set(i, sn.get(i) * (solverOptions.getMaxStep() / newtonLength));
            }
            newtonLength = solverOptions.getMaxStep();
        }
        //initial slope
        double initialSlope = CommonOps_DDRM.dot(g, sn);
        //relative length of sn as calculated in the stopping routine
        double relativeLength = Double.MIN_VALUE;
        for (int i = 0; i < sn.numRows; i++) {
            relativeLength = Math.max(relativeLength, Math.abs(sn.get(i)) / Math.max(Math.abs(x.get(i)), 1.0));
        }
        //minimum allowable step length
        double minLambda = solverOptions.getStepTolerance() / relativeLength;
        /* find the new x values*/
        DMatrixRMaj xPlus = new DMatrixRMaj(x.numRows, x.numCols);
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
                    DMatrixRMaj a = new DMatrixRMaj(new double[][]{{1 / (lambda * lambda), -1 / (lambdaPrevious * lambdaPrevious)}, {-lambdaPrevious / (lambda * lambda), lambda / (lambdaPrevious * lambdaPrevious)}});
                    DMatrixRMaj b = new DMatrixRMaj(2, 1);
                    b.set(0, 0, newFunctionNorm - initialFunctionNorm - lambda * initialSlope);
                    b.set(1, 0, previousFunctionNorm - initialFunctionNorm - lambdaPrevious * initialSlope);
                    DMatrixRMaj multiplication = new DMatrixRMaj(2, 1);
                    CommonOps_DDRM.mult(a, b, multiplication);
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
