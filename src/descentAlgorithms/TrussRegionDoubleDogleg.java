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


public class TrussRegionDoubleDogleg implements StepAlgorithm {

    private boolean firstDog;
    private double cauchyLen;
    private double etta;
    private DMatrixRMaj s;
    private DMatrixRMaj ssd;
    private DMatrixRMaj v;
    private DMatrixRMaj xPlus;
    private DMatrixRMaj xPrev;
    private double fPrev;
    private static final int STATUS_SUCCESSFULL = 0;
    private static final int STATUS_FAILED = 1;
    private static final int STATUS_REDUCE_DELTA = 2;
    private static final int STATUS_DOUBLE_DELTA = 3;
    private static final int STATUS_FIRST_STEP = 4;
    private int solverStatus = STATUS_FIRST_STEP;
    private boolean maxStepTaken;

    @Override
    public boolean isSolverFailed() {
        return solverStatus != STATUS_FAILED;
    }

    @Override
    public boolean isMaxStepTaken() {
        return maxStepTaken;
    }


    /**
     * Main function for trust region double dogleg solver Find an xp such that
     * f(xp) < f(x)+alpha*g^T*(xp-x) @param g gradie
     *
     * nt vector
     * @param g
     * @param sn newton step vector
     * @param x current x vector
     * @param lowerTriangleR lower half of the decomposition
     * @param solverOptions options for the solver
     * @param solver class calling the function
     * @return
     */
    @Override
    public DMatrixRMaj solve(DMatrixRMaj g, DMatrixRMaj sn, DMatrixRMaj x, DMatrixRMaj lowerTriangleR, Options solverOptions, Solver solver) {
        //find x on the double dogleg curve
        //f(x)<= f(x)+alpha*gt*(x+-xc)
        solverStatus = STATUS_FIRST_STEP;
        firstDog = true;
        s = new DMatrixRMaj(g.numRows, 1);
        ssd = new DMatrixRMaj(g.numRows, 1);
        v = new DMatrixRMaj(g.numRows, 1);
        xPlus = x.copy();
        fPrev = solver.functionNorm(x);
        DMatrixRMaj dummySn = new DMatrixRMaj(sn);
        double newtonLength = NormOps_DDRM.normP2(dummySn);

        while (solverStatus >= STATUS_REDUCE_DELTA) {
            boolean newtonTaken = dogStep(newtonLength, g, sn, lowerTriangleR, solverOptions);
            trustRegup(newtonTaken, g, x, lowerTriangleR, solverOptions, solver);
        }
        return xPlus;
    }

    /**
     * Find an approximate solution to minimize g^T*s+1/2*s^TLL^Ts subject to
     * ||Dxs|| <= delta
     * @
     *
     * param newtonLength magnitude of newton step
     * @param g gradient vector
     * @param sn newton step vector
     * @param lowerTriangleR lower half of the decomposition
     * @param solverOptions options for the solver
     * @return return if newton step is preferred
     */
    private boolean dogStep(double newtonLength, DMatrixRMaj g, DMatrixRMaj sn, DMatrixRMaj lowerTriangleR, Options solverOptions) {
        //find approximate solution to
        //minimize gT*s+1/2sTLLTs
        boolean newTaken;
        if (newtonLength <= solverOptions.getTrussRegionRadius()) {
            //s is newton step
            newTaken = true;
            s = sn.copy();
            solverOptions.setTrussRegionRadius(newtonLength);
        } else {
            //newton step too long ,s on double dogleg curve
            newTaken = false;
            if (firstDog) {
                //calculate double dogle curve
                firstDog = false;
                //alpha=||Dx^-1 g||^2
                DMatrixRMaj dummyG = new DMatrixRMaj(g);
                double alpha = Math.pow(NormOps_DDRM.normP2(dummyG), 2);
                //implement b = || L^tDx^-2g||^2
                double beta = 0.0;
                for (int i = 0; i < g.numRows; i++) {
                    double temp = 0.0;
                    for (int j = i; j < g.numRows; j++) {
                        temp += lowerTriangleR.get(j, i) * g.get(j, 0);
                    }
                    beta += temp * temp;
                }
                //ssd is Cauchy step in scaled metric
                for (int i = 0; i < g.numRows; i++) {
                    ssd.set(i, 0, -(alpha / beta)  * g.get(i, 0));
                }
                cauchyLen = alpha * Math.sqrt(alpha) / beta;
                double dummy = 0.0;
                for (int i = 0; i < g.numRows; i++) {
                    dummy += g.get(i, 0) * sn.get(i, 0);
                }
                etta = 0.2 + (0.8 * alpha * alpha / (beta * Math.abs(dummy)));
                for (int i = 0; i < g.numRows; i++) {
                    v.set(i, 0, etta  * sn.get(i, 0) - ssd.get(i, 0));
                }
                if (solverOptions.getTrussRegionRadius() == -1.0) {
                    //no user delta
                    solverOptions.setTrussRegionRadius(Math.min(cauchyLen, solverOptions.getMaxStep()));
                }
            }
            //calculate double dogleg curve
            if (etta * newtonLength <= solverOptions.getTrussRegionRadius()) {
                //take partial step in newton direction
                for (int i = 0; i < sn.numRows; i++) {
                    s.set(i, 0, solverOptions.getTrussRegionRadius() / newtonLength * sn.get(i, 0));
                }
            } else if (cauchyLen >= solverOptions.getTrussRegionRadius()) {
                //take steepest direction
                for (int i = 0; i < sn.numRows; i++) {
                    s.set(i, 0, solverOptions.getTrussRegionRadius() / cauchyLen * ssd.get(i, 0) );
                }
            } else {
                //calculate  convex combination2
                double temp = CommonOps_DDRM.dot(v, ssd);
                double tempV = CommonOps_DDRM.dot(v, v);
                double lambda = (-temp + Math.sqrt(temp * temp - tempV * (cauchyLen * cauchyLen - solverOptions.getTrussRegionRadius() * solverOptions.getTrussRegionRadius()))) / tempV;
                for (int i = 0; i < sn.numRows; i++) {
                    s.set(i, 0, (ssd.get(i, 0) + lambda * v.get(i, 0)) );
                }
            }
        }
        return newTaken;
    }

    /**
     * Given a step s produced by the double dogstep decide whether xp=x+s
     * should be accepted or not if accepted, choose the initial trust region
     * for next iteration if not, decrease or increase trust region for current
     * iteration
     *
     * @param newtonTaken if newton step is taken or not
     * @param g gradient vector
     * @param x x vector
     * @param lowerTriangleR lower half of the decomposed matrix
     * @param solverOptions options for the solver
     * @param solver class calling the solver
     */
    private void trustRegup(boolean newtonTaken, DMatrixRMaj g, DMatrixRMaj x, DMatrixRMaj lowerTriangleR, Options solverOptions, Solver solver) {
        //given a step s produced by dogleg algorithm decide whether x+ is
        //accepted or not
        maxStepTaken = false;
        double alpha = 1e-4;
        //stepLength=||S||
        double stepLength = NormOps_DDRM.normP2(s);
        //xPlus=x+s
        CommonOps_DDRM.add(x, s, xPlus);
        //deltaF=F(xPlus)-F(x)
        double fPlus = solver.functionNorm(xPlus);
        double deltaF = fPlus - solver.functionNorm(x);
        //initSlope=g*s
        double initSlope = CommonOps_DDRM.dot(g, s);
        if (solverStatus != STATUS_DOUBLE_DELTA) {
            fPrev = 0.0;
        }
        if (solverStatus == STATUS_DOUBLE_DELTA && (fPlus >= fPrev || deltaF > alpha * initSlope)) {
            //reset xp to xpre and terminate global step
            solverStatus = STATUS_SUCCESSFULL;
            xPlus = xPrev.copy();
            solverOptions.setTrussRegionRadius(solverOptions.getTrussRegionRadius()/2.0);
        } else if (deltaF >= alpha * initSlope) {
            //f(xPlus) too large
            double relLength = Double.MIN_VALUE;
            for (int i = 0; i < g.numRows; i++) {
                relLength = Math.max(relLength, Math.abs(s.get(i, 0)) /  Math.max(Math.abs(xPlus.get(i)),1.0) );
            }
            if (relLength < solverOptions.getStepTolerance()) {
                //xp-xc too small, terminate
                solverStatus = STATUS_FAILED;
                xPlus = x.copy();
            } else {
                //reduce delta, continue
                solverStatus = STATUS_REDUCE_DELTA;
                double deltaTemp = -(initSlope * stepLength) / (2 * (deltaF - initSlope));
                if (deltaTemp < 0.1 * solverOptions.getTrussRegionRadius()) {
                    solverOptions.setTrussRegionRadius(solverOptions.getTrussRegionRadius()*0.1);
                } else if (deltaTemp > 0.5 * solverOptions.getTrussRegionRadius()) {
                    solverOptions.setTrussRegionRadius(solverOptions.getTrussRegionRadius()*0.5);
                } else {
                    solverOptions.setTrussRegionRadius(deltaTemp);
                }
            }
        } else {
            //f(xP) sufficiently small
            double deltaFPred = initSlope;
            for (int i = 0; i < s.numRows; i++) {
                double temp = 0.0;
                for (int j = i; j < s.numRows; j++) {
                    temp += lowerTriangleR.get(j, i) * s.get(j, 0);
                }
                deltaFPred += (temp * temp / 2);
            }
            if (solverStatus != STATUS_REDUCE_DELTA && (Math.abs(deltaFPred - deltaF) <= Math.abs(deltaF) * 0.1 || (deltaF) <= (initSlope)) && !newtonTaken && solverOptions.getTrussRegionRadius() <= 0.99 * solverOptions.getMaxStep()) {
                //double delta and continue
                solverStatus = STATUS_DOUBLE_DELTA;
                xPrev = xPlus.copy();
                fPrev = fPlus;
                solverOptions.setTrussRegionRadius(Math.min(2*solverOptions.getTrussRegionRadius(), solverOptions.getMaxStep()));
            } else {
                //accept xP as the new iterate
                solverStatus = STATUS_SUCCESSFULL;
                if (stepLength > 0.99 * solverOptions.getMaxStep()) {
                    maxStepTaken = true;
                }
                if (deltaF >= 0.1 * deltaFPred) {
                    solverOptions.setTrussRegionRadius(solverOptions.getTrussRegionRadius()/2.0);
                } else if (deltaF <= 0.75 * deltaFPred) {
                    solverOptions.setTrussRegionRadius(Math.min(2*solverOptions.getTrussRegionRadius(), solverOptions.getMaxStep()));
                }
            }
        }
    }

}
