/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package solverImplementation;

import optimization.functionImplementation.Options;
import org.ejml.data.DenseMatrix64F;
import org.ejml.factory.LinearSolverFactory;
import org.ejml.interfaces.decomposition.CholeskyDecomposition;
import org.ejml.interfaces.decomposition.QRDecomposition;
import org.ejml.interfaces.linsol.LinearSolver;
import org.ejml.ops.CommonOps;
import org.ejml.ops.NormOps;

/**
 *
 * @author O. Anil Turkkan <turkkan.1@osu.edu> - anilturkkan.com
 */
public class NewtonStep {

    private static DenseMatrix64F lowerTriangleR;
    public static boolean restart;

    public static DenseMatrix64F getLowerTriangleR() {
        return lowerTriangleR;
    }

    public static DenseMatrix64F newtonStep(DenseMatrix64F jacobian, DenseMatrix64F fx, DenseMatrix64F g, Options solverOptions) {
        /*Calculate QR decomposition of DfJ*/
        DenseMatrix64F dummyJacobian = jacobian.copy();
        //Df*J
        for (int i = 0; i < dummyJacobian.numRows; i++) {
            for (int j = 0; j < dummyJacobian.numCols; j++) {
                dummyJacobian.set(i, j, dummyJacobian.get(i, j) * solverOptions.getTypicalF().get(i));
            }
        }
        //qr decomposition
        LinearSolver<DenseMatrix64F> qrSolver = LinearSolverFactory.qr(dummyJacobian.numRows, dummyJacobian.numCols);
        qrSolver.setA(dummyJacobian);
        /*get the condition number*/
        double conditionNumber = qrSolver.quality();
        //get the condition number
        //double conditionNumber = algebra.cond(rMatrix);
        //double conditionNumber = ;
        //if (MatrixFeatures.rank(dummyJacobian) == solverOptions.getN()&& conditionNumber > 1e3*Math.sqrt(solverOptions.getMachineEpsilon())) {
        // if (conditionNumber > 1e3 * Math.sqrt(solverOptions.getMachineEpsilon())) {
        if (conditionNumber > 1e-8) {
            /*calculate the newton step -J*sn=fx*/
            DenseMatrix64F newtonianStep = new DenseMatrix64F(dummyJacobian.numRows, 1);
            qrSolver.solve(fx, newtonianStep);
            CommonOps.elementMult(newtonianStep, solverOptions.getTypicalF());
            CommonOps.changeSign(newtonianStep);
            lowerTriangleR = new DenseMatrix64F(jacobian.numRows, jacobian.numCols);
            ((QRDecomposition) qrSolver.getDecomposition()).getR(lowerTriangleR, true);
            CommonOps.transpose(lowerTriangleR);
            return newtonianStep;
        } else {
            //ill-conditioned or singular jacobian
            /* solve -H*sn=g where
                H=J^T*Sf^2*J+sqrt(n*machineEpsilon)*||J^T*Sf^2*J||*Sx^2
             */
            //H=J^T*Sf^2*J
            DenseMatrix64F h = new DenseMatrix64F(dummyJacobian.numRows, dummyJacobian.numCols);
            CommonOps.multInner(jacobian, h);
            for (int i = 0; i < h.numRows; i++) {
                for (int j = 0; j < h.numCols; j++) {
                    h.set(i, j, h.get(i, j) * solverOptions.getTypicalF().get(i) * solverOptions.getTypicalF().get(i));
                }
            }
            double normH = NormOps.normP1(h);
            for (int i = 0; i < h.numRows; i++) {
                h.set(i, i, h.get(i, i) + Math.sqrt(h.numRows * solverOptions.getMachineEpsilon()) * normH * solverOptions.getTypicalX().get(i, 0) * solverOptions.getTypicalX().get(i, 0));
            }
            //cholesky decomposition
            LinearSolver<DenseMatrix64F> cDecomposition = LinearSolverFactory.chol(dummyJacobian.numRows);
            cDecomposition.setA(h);
            DenseMatrix64F newtonianStep = new DenseMatrix64F(dummyJacobian.numRows, 1);
            cDecomposition.solve(g, newtonianStep);
            CommonOps.changeSign(newtonianStep);
            lowerTriangleR = new DenseMatrix64F(jacobian.numRows, jacobian.numCols);
            ((CholeskyDecomposition) cDecomposition.getDecomposition()).getT(lowerTriangleR);
            if (!((CholeskyDecomposition) cDecomposition.getDecomposition()).isLower()) {
                CommonOps.transpose(lowerTriangleR);
            }
            return newtonianStep;
        }
    }

    public static DenseMatrix64F newtonStepBFGS(DenseMatrix64F jacobian, DenseMatrix64F fx, DenseMatrix64F g, Options solverOptions) {
        if (restart) {
            NewtonStep.restart=false;
            /*Calculate QR decomposition of DfJ*/
            DenseMatrix64F dummyJacobian = jacobian.copy();
            //Df*J
            for (int i = 0; i < dummyJacobian.numRows; i++) {
                for (int j = 0; j < dummyJacobian.numCols; j++) {
                    dummyJacobian.set(i, j, dummyJacobian.get(i, j) * solverOptions.getTypicalF().get(i));
                }
            }
            //qr decomposition
            LinearSolver<DenseMatrix64F> qrSolver = LinearSolverFactory.qr(dummyJacobian.numRows, dummyJacobian.numCols);
            qrSolver.setA(dummyJacobian);
            /*get the condition number*/
            double conditionNumber = qrSolver.quality();
            //get the condition number
            if (conditionNumber > 1e-8) {
                /*calculate the newton step -J*sn=fx*/
                BFGS.l = new DenseMatrix64F(jacobian.numRows, jacobian.numCols);
                ((QRDecomposition) qrSolver.getDecomposition()).getR(BFGS.l, true);
                BFGS.h=  new DenseMatrix64F(jacobian.numRows, jacobian.numCols);
                ((QRDecomposition) qrSolver.getDecomposition()).getQ(BFGS.h, true);
                CommonOps.transpose(BFGS.h);
            } else {
                //ill-conditioned or singular jacobian
                /* solve -H*sn=g where
                H=J^T*Sf^2*J+sqrt(n*machineEpsilon)*||J^T*Sf^2*J||*Sx^2
                 */
                //H=J^T*Sf^2*J
                DenseMatrix64F h = new DenseMatrix64F(dummyJacobian.numRows, dummyJacobian.numCols);
                CommonOps.multInner(jacobian, h);
                for (int i = 0; i < h.numRows; i++) {
                    for (int j = 0; j < h.numCols; j++) {
                        h.set(i, j, h.get(i, j) * solverOptions.getTypicalF().get(i) * solverOptions.getTypicalF().get(i));
                    }
                }
                double normH = NormOps.normP1(h);
                for (int i = 0; i < h.numRows; i++) {
                    h.set(i, i, h.get(i, i) + Math.sqrt(h.numRows * solverOptions.getMachineEpsilon()) * normH * solverOptions.getTypicalX().get(i, 0) * solverOptions.getTypicalX().get(i, 0));
                }
                //cholesky decomposition
                LinearSolver<DenseMatrix64F> cDecomposition = LinearSolverFactory.chol(dummyJacobian.numRows);
                cDecomposition.setA(h);
                DenseMatrix64F newtonianStep = new DenseMatrix64F(dummyJacobian.numRows, 1);
                cDecomposition.solve(g, newtonianStep);
                CommonOps.changeSign(newtonianStep);
                BFGS.l = new DenseMatrix64F(jacobian.numRows, jacobian.numCols);
                ((CholeskyDecomposition) cDecomposition.getDecomposition()).getT(BFGS.l);
                if (!((CholeskyDecomposition) cDecomposition.getDecomposition()).isLower()) {
                    CommonOps.transpose(BFGS.l);
                }
                return newtonianStep;
            }
        } 
        /* sn= -J*Df*F */
        DenseMatrix64F newtonianStep = new DenseMatrix64F(g.numRows, 1);
        for(int i=0;i<g.numRows;i++){
            newtonianStep.set(i,0,0.0);
            for(int j=0;j<g.numRows;j++){
                newtonianStep.set(i,0,newtonianStep.get(i,0)-BFGS.h.get(i,j)*solverOptions.getTypicalF().get(j,0)*fx.get(j, 0));
            }
        }
        /* calculate l^-1*sn */
        newtonianStep.set(g.numRows-1, 0,newtonianStep.get(g.numRows-1,0)/BFGS.l.get(g.numRows-1,g.numRows-1));
        for(int i=g.numRows-2;i>=0;i--){
            double temp=0.0;
            for(int j=i+1;j<g.numRows;j++){
                temp += BFGS.l.get(i,j)*newtonianStep.get(j,0);
            }
            newtonianStep.set(i, 0,(newtonianStep.get(i,0)-temp)/BFGS.l.get(i,i));
        }
        return newtonianStep;
    }

}
