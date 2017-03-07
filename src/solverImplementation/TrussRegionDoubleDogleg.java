/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package solverImplementation;

import optimization.functionImplementation.Options;
import org.ejml.data.DenseMatrix64F;
import org.ejml.ops.CommonOps;
import org.ejml.ops.NormOps;
import solvers.Solver;

/**
 *
 * @author O. Anil Turkkan <turkkan.1@osu.edu> - anilturkkan.com
 */
public class TrussRegionDoubleDogleg {

    private static boolean firstDog;
    private static double delta = -1;
    private static double cauchyLen;
    private static double etta;
    private static DenseMatrix64F s;
    private static DenseMatrix64F ssd;
    private static DenseMatrix64F v;
    private static DenseMatrix64F xPlus;
    private static DenseMatrix64F xPrev;
    private static double fPrev;
    private static final int STATUS_SUCCESSFULL = 0;
    private static final int STATUS_FAILED = 1;
    private static final int STATUS_REDUCE_DELTA = 2;
    private static final int STATUS_DOUBLE_DELTA = 3;
    private static final int STATUS_FIRST_STEP = 4;
    private static int solverStatus = STATUS_FIRST_STEP;
    private static boolean maxStepTaken;

    public static boolean getSolverStatus() {
        return solverStatus != STATUS_FAILED;
    }

    public static boolean isMaxTaken() {
        return maxStepTaken;
    }

    public static DenseMatrix64F dogDriver(DenseMatrix64F g, DenseMatrix64F sn, DenseMatrix64F x, DenseMatrix64F lowerTriangleR, Options solverOptions, Solver solver) {
        //find x on the double dogleg curve
        //f(x)<= f(x)+alpha*gt*(x+-xc)
        solverStatus = STATUS_FIRST_STEP;
        firstDog = true;
        s = new DenseMatrix64F(g.numRows, 1);
        ssd = new DenseMatrix64F(g.numRows, 1);
        v = new DenseMatrix64F(g.numRows, 1);
        xPlus = x.copy();
        fPrev = solver.functionNorm(x);
        DenseMatrix64F dummySn = new DenseMatrix64F(sn);
        CommonOps.elementMult(dummySn, solverOptions.getTypicalX());
        double newtonLength = NormOps.normP2(dummySn);

        while (solverStatus >= STATUS_REDUCE_DELTA) {
            boolean newtonTaken = dogStep(newtonLength, g, sn, lowerTriangleR, solverOptions);
            trustRegup(newtonTaken, g, x, lowerTriangleR, solverOptions, solver);
        }
        return xPlus;
    }

    private static boolean dogStep(double newtonLength, DenseMatrix64F g, DenseMatrix64F sn, DenseMatrix64F lowerTriangleR, Options solverOptions) {
        //find approximate solution to
        //minimize gT*s+1/2sTLLTs
        boolean newTaken;
        if (newtonLength <= delta) {
            //s is newton step
            newTaken = true;
            s = sn.copy();
            delta = newtonLength;
        } else {
            //newton step too long ,s on double dogleg curve
            newTaken = false;
            if (firstDog) {
                //calculate double dogle curve
                firstDog = false;
                //alpha=||Dx^-1 g||^2
                DenseMatrix64F dummyG = new DenseMatrix64F(g);
                CommonOps.elementDiv(dummyG, solverOptions.getTypicalX());
                double alpha = Math.pow(NormOps.normP2(dummyG), 2);
                //implement b = || L^tDx^-2g||^2
                double beta = 0.0;
                for (int i = 0; i < g.numRows; i++) {
                    double temp = 0.0;
                    for (int j  = i; j < g.numRows; j++) {
                        temp += lowerTriangleR.get(j, i) * g.get(j, 0) / (solverOptions.getTypicalX().get(j, 0) * solverOptions.getTypicalX().get(j, 0));
                    }
                    beta += temp * temp;
                }
                //ssd is Cauchy step in scaled metric
                for (int i = 0; i < g.numRows; i++) {
                    ssd.set(i, 0, -(alpha / beta) / solverOptions.getTypicalX().get(i, 0) * g.get(i, 0));
                }
                cauchyLen = alpha * Math.sqrt(alpha) / beta;
                double dummy = 0.0;
                for (int i = 0; i < g.numRows; i++) {
                    dummy += Math.abs(g.get(i, 0) * sn.get(i, 0));
                }
                etta = 0.2 + (0.8 * alpha * alpha / (beta * dummy));
                for (int i = 0; i < g.numRows; i++) {
                    v.set(i, 0, etta * solverOptions.getTypicalX().get(i, 0) * sn.get(i, 0) - ssd.get(i, 0));
                }
                if (delta == -1.0) {
                    //no user delta
                    delta = Math.min(cauchyLen, solverOptions.getMaxStep());
                }
            }
            //calculate double dogleg curve
            if (etta * newtonLength <= delta) {
                //take partial step in newton direction
                for (int i = 0; i < sn.numRows; i++) {
                    s.set(i, 0, delta / newtonLength * sn.get(i, 0));
                }
            } else if (cauchyLen >= delta) {
                //take steepest direction
                for (int i = 0; i < sn.numRows; i++) {
                    s.set(i, 0, delta / cauchyLen * ssd.get(i, 0) / solverOptions.getTypicalX().get(i, 0));
                }
            } else {
                //calculate  convex combination2
                double temp = CommonOps.dot(v, ssd);
                double tempV = CommonOps.dot(v, v);
                double lambda = (-temp + Math.sqrt(temp * temp - tempV * (cauchyLen * cauchyLen - delta * delta))) / tempV;
                for (int i = 0; i < sn.numRows; i++) {
                    s.set(i, 0, (ssd.get(i, 0) + lambda * v.get(i, 0)) / solverOptions.getTypicalX().get(i, 0));
                }
            }
        }
        return newTaken;
    }

    private static void trustRegup(boolean newtonTaken, DenseMatrix64F g, DenseMatrix64F x, DenseMatrix64F lowerTriangleR, Options solverOptions, Solver solver) {
        //given a step s produced by dogleg algorithm decide whether x+ is
        //accepted or not
        maxStepTaken = false;
        double alpha = 1e-4;
        //stepLength=||DxS||
        DenseMatrix64F dummyS = new DenseMatrix64F(s);
        CommonOps.elementDiv(dummyS, solverOptions.getTypicalX());
        double stepLength = NormOps.normP2(s);
        //xPlus=x+s
        CommonOps.add(x, s, xPlus);
        //deltaF=F(xPlus)-F(x)
        double fPlus = solver.functionNorm(xPlus);
        double deltaF = fPlus - solver.functionNorm(x);
        //initSlope=g*s
        double initSlope = CommonOps.dot(g, s);
        if (solverStatus != STATUS_DOUBLE_DELTA) {
            fPrev = 0.0;
        }
        if (solverStatus == STATUS_DOUBLE_DELTA && (fPlus >= fPrev || deltaF > alpha * initSlope)) {
            //reset xp to xpre and terminate global step
            solverStatus = STATUS_SUCCESSFULL;
            xPlus = xPrev.copy();
            delta /= 2.0;
        } else if (deltaF >= alpha * initSlope) {
            //f(xPlus) too large
            double relLength = Double.MIN_VALUE;
            for (int i = 0; i < g.numRows; i++) {
                relLength = Math.max(relLength, Math.abs(s.get(i, 0)) / Math.max(Math.abs(xPlus.get(i, 0)), 1 / solverOptions.getTypicalX().get(i, 0)));
            }
            if (relLength < solverOptions.getStepTolerance()) {
                //xp-xc too small, terminate
                solverStatus = STATUS_FAILED;
                xPlus = x.copy();
            } else {
                //reduce delta, continue
                solverStatus = STATUS_REDUCE_DELTA;
                double deltaTemp = -(initSlope * stepLength) / (2 * (deltaF - initSlope));
                if (deltaTemp < 0.1 * delta) {
                    delta = 0.1 * delta;
                } else if (deltaTemp > 0.5 * delta) {
                    delta = 0.5 * delta;
                } else {
                    delta = deltaTemp;
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
            if (solverStatus != STATUS_REDUCE_DELTA && (Math.abs(deltaFPred - deltaF) <= Math.abs(deltaF) * 0.1 ||  Math.abs(deltaF) <= Math.abs(initSlope)) && newtonTaken == false && delta <= 0.99 * solverOptions.getMaxStep()) {
                //double delta and continue
                solverStatus = STATUS_DOUBLE_DELTA;
                xPrev = xPlus.copy();
                fPrev = fPlus;
                delta = Math.min(2 * delta, solverOptions.getMaxStep());
            } else {
                //accept xP as the new iterate
                solverStatus = STATUS_SUCCESSFULL;
                if (stepLength > 0.99 * solverOptions.getMaxStep()) {
                    maxStepTaken = true;
                }
                if (deltaF >= 0.1 * deltaFPred) {
                    delta /= 2;
                } else if (deltaF <= 0.75 * deltaFPred) {
                    delta = Math.min(2 * delta, solverOptions.getMaxStep());
                }
            }
        }
    }

}
