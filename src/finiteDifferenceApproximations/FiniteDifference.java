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
package finiteDifferenceApproximations;

import optimization.functionImplementation.Options;
import org.ejml.data.DMatrixRMaj;
import optimization.functionImplementation.ObjectiveFunctionNonLinear;
import optimization.functionImplementation.ObjectiveFunctionUnconstrained;


public class FiniteDifference {

    /**
     * Return gradient vector
     *
     * @param x The vector where gradient will be calculated
     * @param equations User defined objective function
     * @param solverOptions User preferences
     * @return analytical gradient if set by options. Else, return central or
     * forward finite difference approximation
     */
    public static DMatrixRMaj getGradient(DMatrixRMaj x, ObjectiveFunctionUnconstrained equations, Options solverOptions) {
        //check if analytical gradients are given
        if (solverOptions.isAnalyticalGradient()) {
            return equations.getG(x);
        }
        if (solverOptions.isCentralDifferenceGradient()) {
            return getCentralDifferenceGradient(x, equations, solverOptions);
        } else {
            return getForwardDifferenceGradient(x, equations, solverOptions);
        }
    }

    /**
     * Approximate the gradient vector by central difference method
     *
     * @param x The vector where gradient will be calculated
     * @param equations User defined objective function
     * @param solverOptions User preferences
     * @return
     */
    public static DMatrixRMaj getCentralDifferenceGradient(DMatrixRMaj x, ObjectiveFunctionUnconstrained equations, Options solverOptions) {
        DMatrixRMaj gradient = new DMatrixRMaj(solverOptions.getN(), 1);
        double sqrDelta = Math.sqrt(solverOptions.getMachineEpsilon());
        double dummyFx; // storing intermadiate function values
        double dummyFx2; // storing intermadiate function values
        //start calculating gradient 
        for (int j = 0; j < solverOptions.getN(); j++) {
            //get the sign of x
            double sign = Math.signum(x.get(j));
            if (sign == 0.0) {
                sign = 1.0;
            }
            //calculate the step size
            double stepSize = sqrDelta *  Math.max(Math.abs(x.get(j)),1.0) * sign;
            //save the initial value
            double temp = x.get(j);
            //calculate x+stepSize
            x.set(j, 0, x.get(j) + stepSize);
            //to reduce finite difference precision errors
            stepSize = x.get(j, 0) - temp;
            //get the value at x+stepsize
            dummyFx = equations.getF(x);
            //calculate x-stepSize
            x.set(j, 0, temp - stepSize);
            //get the value at x-stepsize
            dummyFx2 = equations.getF(x);
            //set the gradient
            gradient.set(j, 0, (dummyFx - dummyFx2) / (2 * stepSize));
            //reset the value
            x.set(j, 0, temp);
        }
        return gradient;
    }

    /**
     * Approximate the gradient vector by forward difference method
     *
     * @param x The vector where gradient will be calculated
     * @param equations User defined objective function
     * @param solverOptions User preferences
     * @return
     */
    public static DMatrixRMaj getForwardDifferenceGradient(DMatrixRMaj x, ObjectiveFunctionUnconstrained equations, Options solverOptions) {
        DMatrixRMaj gradient = new DMatrixRMaj(solverOptions.getN(), 1);
        double sqrDelta = Math.sqrt(solverOptions.getMachineEpsilon());
        double dummyFx; // storing intermadiate function values
        double fx = equations.getF(x);
        //start calculating gradient 
        for (int j = 0; j < solverOptions.getN(); j++) {
            //get the sign of x
            double sign = Math.signum(x.get(j));
            if (sign == 0.0) {
                sign = 1.0;
            }
            //calculate the step size
            double stepSize = sqrDelta *  Math.max(Math.abs(x.get(j)),1.0) * sign;
            //save the initial value
            double temp = x.get(j);
            //calculate x+stepSize
            x.set(j, 0, x.get(j) + stepSize);
            //to reduce finite difference precision errors
            stepSize = x.get(j, 0) - temp;
            //get the value at x+stepsize
            dummyFx = equations.getF(x);
            //set the gradient
            gradient.set(j, 0, (dummyFx - fx) / stepSize);
            //reset the value
            x.set(j, 0, temp);
        }
        return gradient;
    }

    /**
     * Return Hessian matrix
     *
     * @param x The vector where Hessian will be calculated
     * @param equations User defined objective function
     * @param solverOptions User preferences
     * @return analytical Hessian if set by options. Else, return forward finite
     * difference approximation
     */
    public static DMatrixRMaj getHessian(DMatrixRMaj x, ObjectiveFunctionUnconstrained equations, Options solverOptions) {
        //check if analytical gradients are given
        if (solverOptions.isAnalyticalHessian()) {
            return equations.getH(x);
        } else {
            return getForwardDifferenceHessian(x, equations, solverOptions);
        }
    }

    /**
     * Calculate the Hessian matrix using forward difference approximation
     *
     * @param x The vector where Hessian will be calculated
     * @param equations User defined objective function
     * @param solverOptions User preferences
     * @return
     */
    public static DMatrixRMaj getForwardDifferenceHessian(DMatrixRMaj x, ObjectiveFunctionUnconstrained equations, Options solverOptions) {
        DMatrixRMaj hessian = new DMatrixRMaj(solverOptions.getN(), solverOptions.getN());
        double cubeDelta = Math.pow(solverOptions.getMachineEpsilon(), 1.0 / 3.0);
        DMatrixRMaj stepSize = new DMatrixRMaj(solverOptions.getN(), 1);
        double fx = equations.getF(x);
        DMatrixRMaj fxNew = new DMatrixRMaj(solverOptions.getN(), 1);
        //calculate step sizes
        for (int i = 0; i < solverOptions.getN(); i++) {
            //get the sign of x
            double sign = Math.signum(x.get(i));
            if (sign == 0.0) {
                sign = 1.0;
            }
            //calculate the step size
            stepSize.set(i, 0, cubeDelta *  Math.max(Math.abs(x.get(i)),1.0) * sign);
            //save the initial value
            double temp = x.get(i);
            //calculate x+stepSize
            x.set(i, 0, x.get(i) + stepSize.get(i));
            //to reduce finite difference precision errors
            stepSize.set(i, 0, x.get(i, 0) - temp);
            //get the value at x+stepsize
            fxNew.set(i, 0, equations.getF(x));
            //reset the value
            x.set(i, 0, temp);
        }
        //calculate hessian
        for (int i = 0; i < solverOptions.getN(); i++) {
            /* calculate row i of Hessian */
            //save the initial value
            double temp = x.get(i);
            //calculate x+2*stepSize
            x.set(i, 0, x.get(i) + 2 * stepSize.get(i));
            //get the value at x+2*stepSize
            double dummyFx = equations.getF(x);
            //Hii
            hessian.set(i, i, ((fx - fxNew.get(i, 0)) + (dummyFx - fxNew.get(i, 0))) / (stepSize.get(i, 0) * stepSize.get(i, 0)));
            //calculate x+stepSize
            x.set(i, 0, temp + stepSize.get(i));
            for (int j = i + 1; j < solverOptions.getN(); j++) {
                //calculate H(i,j)
                double temp2 = x.get(j);
                //calculate x+stepSize
                x.set(j, 0, x.get(j) + stepSize.get(j));
                //get the value at x+2*stepSize
                double dummyFx2 = equations.getF(x);
                //Hij - Hji
                hessian.set(i, j, ((fx - fxNew.get(i, 0)) + (dummyFx2 - fxNew.get(j, 0))) / (stepSize.get(i, 0) * stepSize.get(j, 0)));
                hessian.set(j, i, hessian.get(i, j));
                //calculate x+stepSize
                x.set(j, 0, temp2);
            }
            //reset the value
            x.set(i, 0, temp);
        }
        return hessian;
    }

    /**
     * Returns the Jacobian at a given x If analytical Jacobian flag is true,
     * user calculated Jacobian will be returned Otherwise, forward finite
     * difference Jacobian will be calculated
     *
     * @param x The point where Jacobian will be calculated
     * @param equations User defined objective function
     * @param solverOptions User preferences
     * @return
     */
    public static DMatrixRMaj getJacobian(DMatrixRMaj x, ObjectiveFunctionNonLinear equations, Options solverOptions) {
        //check if analytical jacobian
        if (solverOptions.isAnalyticalJacobian()) {
            return equations.getJ(x);
        } else {
            //finite difference Jacobian
            DMatrixRMaj finiteDifferenceJacobian = new DMatrixRMaj(solverOptions.getN(), solverOptions.getN());
            DMatrixRMaj fxx = equations.getF(x);
            DMatrixRMaj dummyFx;
            double sqrDelta = Math.sqrt(solverOptions.getMachineEpsilon());
            for (int j = 0; j < solverOptions.getN(); j++) {
                //calculate column j of J
                //get the sign of x
                double sign = Math.signum(x.get(j));
                if (sign == 0.0) {
                    sign = 1.0;
                }
                //calculate the step size
                double stepSize = sqrDelta * Math.max(Math.abs(x.get(j)),1.0) * sign;
                //save the initial value
                double temp = x.get(j);
                //calculate x+stepSize
                x.set(j, 0, x.get(j) + stepSize);
                //to reduce finite difference precision errors
                stepSize = x.get(j, 0) - temp;
                //get the value at x+stepsize
                dummyFx = equations.getF(x);
                for (int i = 0; i < solverOptions.getN(); i++) {
                    finiteDifferenceJacobian.set(i, j, (dummyFx.get(i) - fxx.get(i)) / stepSize);
                }
                //reset the value
                x.set(j, 0, temp);
            }
            return finiteDifferenceJacobian;
        }
    }

}
