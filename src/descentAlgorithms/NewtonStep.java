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
import org.ejml.dense.row.factory.LinearSolverFactory_DDRM;
import org.ejml.interfaces.decomposition.CholeskyDecomposition;
import org.ejml.interfaces.decomposition.QRDecomposition;
import org.ejml.interfaces.linsol.LinearSolver;
import org.ejml.data.DMatrixRMaj;
import org.ejml.dense.row.CommonOps_DDRM;
import org.ejml.dense.row.NormOps_DDRM;


public class NewtonStep {

    private  DMatrixRMaj lowerTriangleR;
    public  boolean restart;

    /**
     * Get lower triangle decomposition of the Jacobian
     * @return lower triangle Jacobian matrix
     */
    public  DMatrixRMaj getLowerTriangleR() {
        return lowerTriangleR;
    }

    /**
     * Factor a model Jacobian and calculate Newton step
     * @param jacobian Jacobian matrix
     * @param fx function values at current x
     * @param g current x
     * @param solverOptions solver options
     * @return new x vector 
     */
    public  DMatrixRMaj newtonStep(DMatrixRMaj jacobian, DMatrixRMaj fx, DMatrixRMaj g, Options solverOptions) {
        /*Calculate QR decomposition of J*/
        DMatrixRMaj dummyJacobian = jacobian.copy();
        //qr decomposition
        LinearSolver<DMatrixRMaj> qrSolver = LinearSolverFactory_DDRM.qr(dummyJacobian.numRows, dummyJacobian.numCols);
        qrSolver.setA(dummyJacobian);
        /*get the condition number*/
        double conditionNumber = qrSolver.quality();
        //get the condition number
        if (conditionNumber > 1e-12) {
            /*calculate the newton step -J*sn=fx*/
            DMatrixRMaj newtonianStep = new DMatrixRMaj(dummyJacobian.numRows, 1);
            qrSolver.solve(fx, newtonianStep);
            CommonOps_DDRM.changeSign(newtonianStep);
            lowerTriangleR = new DMatrixRMaj(jacobian.numRows, jacobian.numCols);
            ((QRDecomposition) qrSolver.getDecomposition()).getR(lowerTriangleR, true);
            CommonOps_DDRM.transpose(lowerTriangleR);
            return newtonianStep;
        } else {
            //ill-conditioned or singular jacobian
            /* solve -H*sn=g where
                H=J^T*Sf^2*J+sqrt(n*machineEpsilon)*||J^T*Sf^2*J||*Sx^2
             */
            //H=J^T*Sf^2*J
            DMatrixRMaj h = new DMatrixRMaj(dummyJacobian.numRows, dummyJacobian.numCols);
            CommonOps_DDRM.multInner(jacobian, h);
            double normH = NormOps_DDRM.normP1(h);
            for (int i = 0; i < h.numRows; i++) {
                h.set(i, i, h.get(i, i) + Math.sqrt(h.numRows * solverOptions.getMachineEpsilon()) * normH);
            }
            //cholesky decomposition
            LinearSolver<DMatrixRMaj> cDecomposition = LinearSolverFactory_DDRM.chol(dummyJacobian.numRows);
            cDecomposition.setA(h);
            DMatrixRMaj newtonianStep = new DMatrixRMaj(dummyJacobian.numRows, 1);
            cDecomposition.solve(g, newtonianStep);
            CommonOps_DDRM.changeSign(newtonianStep);
            lowerTriangleR = new DMatrixRMaj(jacobian.numRows, jacobian.numCols);
            ((CholeskyDecomposition) cDecomposition.getDecomposition()).getT(lowerTriangleR);
            if (!((CholeskyDecomposition) cDecomposition.getDecomposition()).isLower()) {
                CommonOps_DDRM.transpose(lowerTriangleR);
            }
            return newtonianStep;
        }
    }
}
