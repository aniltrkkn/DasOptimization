/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
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
/**
 *
 * @author O. Anil Turkkan <turkkan.1@osu.edu> - anilturkkan.com
 */
public class NewtonStep {

    private  DMatrixRMaj lowerTriangleR;
    public  boolean restart;

    public  DMatrixRMaj getLowerTriangleR() {
        return lowerTriangleR;
    }

    public  DMatrixRMaj newtonStep(DMatrixRMaj jacobian, DMatrixRMaj fx, DMatrixRMaj g, Options solverOptions) {
        /*Calculate QR decomposition of J*/
        DMatrixRMaj dummyJacobian = jacobian.copy();
        //qr decomposition
        LinearSolver<DMatrixRMaj> qrSolver = LinearSolverFactory_DDRM.qr(dummyJacobian.numRows, dummyJacobian.numCols);
        qrSolver.setA(dummyJacobian);
        /*get the condition number*/
        double conditionNumber = qrSolver.quality();
        //get the condition number
        //double conditionNumber = algebra.cond(rMatrix);
        //double conditionNumber = ;
        //if (MatrixFeatures.rank(dummyJacobian) == solverOptions.getN()&& conditionNumber > 1e3*Math.sqrt(solverOptions.getMachineEpsilon())) {
        // if (conditionNumber > 1e3 * Math.sqrt(solverOptions.getMachineEpsilon())) {
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
