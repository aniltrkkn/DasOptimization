/*
 * To change this license header, choose license Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package solverImplementation;

import optimization.functionImplementation.Options;
import org.ejml.data.DenseMatrix64F;
import org.ejml.ops.CommonOps;
import org.ejml.ops.NormOps;

/**
 *
 * @author O. Anil Turkkan <turkkan.1@osu.edu> - anilturkkan.com
 */
public class BFGS {
    
    private static DenseMatrix64F s;
    private static DenseMatrix64F y;
    public static DenseMatrix64F l;
    public static DenseMatrix64F h;
    
    public static void CholeskyDecomposition(DenseMatrix64F x, DenseMatrix64F xPlus,
            DenseMatrix64F g, DenseMatrix64F gPlus, Options solverOptions) {
        s = new DenseMatrix64F(x.numRows, 1);
        y = new DenseMatrix64F(x.numRows, 1);
        l = new DenseMatrix64F(x.numRows, x.numCols);
        DenseMatrix64F t = new DenseMatrix64F(x.numRows, 1);
        DenseMatrix64F u = new DenseMatrix64F(x.numRows, 1);
        /* assign s and y*/
        CommonOps.subtract(xPlus, x, s);
        CommonOps.subtract(gPlus, g, y);
        /*check if we should skip the update */
        double temp1 = CommonOps.dot(y, s);
        if (temp1 >= Math.sqrt(solverOptions.getMachineEpsilon()) * NormOps.fastNormP2(s) * NormOps.fastNormP2(y)) {
            /* t=l^Ts */
            for (int i = 0; i < x.numRows; i++) {
                t.set(i, 0, 0);
                for (int j = 0; j < x.numRows; j++) {
                    t.set(i, 0, t.get(i, 0) + l.get(j, i) * s.get(j, 0));
                }
            }
            /* s^Tll^Ts*/
            double temp2 = CommonOps.dot(t, t);
            /* sqrt(y^Ts/s^Tll^Ts */
            double alpha = Math.sqrt(temp1 / temp2);
            /* tolerances */
            double tolerance;
            if (solverOptions.isAnalyticalGradient()) {
                tolerance = solverOptions.getMachineEpsilon();
            } else {
                tolerance = Math.sqrt(solverOptions.getMachineEpsilon());
            }
            boolean skipUpdate = true;
            for (int i = 0; i < x.numRows; i++) {
                /* (ll^Ts)[i] */
                double temp3 = 0;
                for (int j = 0; j <= i; j++) {
                    temp3 += l.get(i, j) * t.get(j, 0);
                }
                if (Math.abs(y.get(i, 0) - temp3) >= tolerance * Math.max(Math.abs(g.get(i, 0)), Math.abs(gPlus.get(i, 0)))) {
                    skipUpdate = false;
                }
                u.set(i, 0, y.get(i, 0) - alpha * temp3);
            }
            /* skip the update or update the l */
            if (!skipUpdate) {
                /* 1/sqrt(y^Ts*s^Tll^Ts) */
                double temp3 = 1.0 / Math.sqrt(temp1 * temp2);
                for (int i = 0; i < x.numRows; i++) {
                    t.set(i, 0, t.get(i, 0) * temp3);
                }
                /* copy l^t into upper triangle of l */
                for (int i = 1; i < x.numRows; i++) {
                    for (int j = 0; j < i; j++) {
                        l.set(j, i, l.get(i, j));
                    }
                }
                /* calculate QR factorization of L^T+ut^T */
                BFGS.qrUpdate(t, u,false);
                /* copy transpose of upper triangle of l into l */
                for (int i = 1; i < x.numRows; i++) {
                    for (int j = 0; j < i - 1; j++) {
                        l.set(i, j, l.get(j, i));
                    }
                }
            }
        }
    }
    
    private static void qrUpdate(DenseMatrix64F u, DenseMatrix64F v, boolean nonlinear) {
        for (int i = 1; i < u.numRows; i++) {
            BFGS.l.set(i, i - 1, 0.0);
        }
        /* find the largest k such that u[k] != 0 */
        int k = u.numRows - 1;
        while (u.get(k, 0) == 0 && k > 0) {
            k--;
        }
        /* transform R+uv^T to upper Hessenberg */
        for (int i = k - 1; i >= 0; i--) {
            BFGS.jacobiRotation(u.numRows, i, u.get(i, 0), -u.get(i + 1, 0),nonlinear);
            if (u.get(i, 0) == 0) {
                u.set(i, 0, Math.abs(u.get(i + 1, 0)));
            } else {
                u.set(i, 0, Math.sqrt(u.get(i, 0) * u.get(i, 0) + u.get(i + 1, 0) * u.get(i + 1, 0)));
            }
        }
        for (int j = 0; j < u.numRows; j++) {
            BFGS.l.set(0, j, BFGS.l.get(0, j) + u.get(0, 0) * v.get(j, 0));
        }
        /*transform upper Hessenberg matrix to upper triangular */
        for (int i = 0; i <= k - 1; i++) {
            BFGS.jacobiRotation(u.numRows, i, BFGS.l.get(i, i), -BFGS.l.get(i + 1, i),nonlinear);
        }
    }
    
    private static void jacobiRotation(int n, int i, double a, double b, boolean nonlinear) {
        double c;
        double ss;
        if (a == 0) {
            c = 0;
            ss = Math.signum(b);
        } else {
            double den = Math.sqrt(a * a + b * b);
            c = a / den;
            ss = b / den;
        }
        /* premultiple l by Jacobi rotation */
        for (int j = 0; j < n; j++) {
            double yy = BFGS.l.get(i, j);
            double w = BFGS.l.get(i + 1, j);
            BFGS.l.set(i, j, c * yy - ss * w);
            BFGS.l.set(i + 1, j, ss * yy + c * w);
        }
        if(nonlinear){
            /*premultiply H by jacobi rotation */
            for (int j = 0; j < n; j++) {
                double yy = BFGS.h.get(i, j);
                double w = BFGS.h.get(i + 1, j);
                BFGS.h.set(i, j, c * yy - ss * w);
                BFGS.h.set(i + 1, j, ss * yy + c * w);
            }
        }
    }
    
    public static void updateJacobian(DenseMatrix64F x, DenseMatrix64F xPlus,
            DenseMatrix64F g, DenseMatrix64F gPlus, Options solverOptions) {
        s = new DenseMatrix64F(x.numRows, 1);
        y = new DenseMatrix64F(x.numRows, 1);
        h = new DenseMatrix64F(x.numRows, x.numCols);
        DenseMatrix64F t = new DenseMatrix64F(x.numRows, 1);
        /* assign s and y*/
        CommonOps.subtract(xPlus, x, s);
        CommonOps.subtract(gPlus, g, y);
        /*check if we should skip the update */
        double temp1 = CommonOps.dot(y, s);
        if (temp1 >= Math.sqrt(solverOptions.getMachineEpsilon()) * NormOps.fastNormP2(s) * NormOps.fastNormP2(y)) {
            /* tolerances */
            double tolerance;
            if (solverOptions.isAnalyticalGradient()) {
                tolerance = solverOptions.getMachineEpsilon();
            } else {
                tolerance = Math.sqrt(solverOptions.getMachineEpsilon());
            }
            boolean skipUpdate = true;
            /* t[i]=Hs[i] */
            for (int i = 0; i < x.numRows; i++) {
                t.set(i, 0, 0.0);
                for (int j = 0; j <= i; j++) {
                    t.set(i, 0, t.get(i, 0) + h.get(j, i) * s.get(j, 0));
                }
                for (int j = i + 1; j < x.numRows; j++) {
                    t.set(i, 0, t.get(i, 0) + h.get(i, j) * s.get(j, 0));
                }
                if (Math.abs(y.get(i, 0) - t.get(i, 0)) >= tolerance * Math.max(Math.abs(g.get(i, 0)), Math.abs(gPlus.get(i, 0)))) {
                    skipUpdate = false;
                }
            }
            /* skip the update or update the h */
            if (!skipUpdate) {
                /* s^THs */
                double temp2 = CommonOps.dot(s, t);
                for (int i = 0; i < x.numRows; i++) {
                    for (int j = i; j < x.numRows; j++) {
                        h.set(i, j, h.get(i, j) + (y.get(i, 0) * y.get(j, 0)) / temp1 - (t.get(i, 0) * t.get(j, 0)) / temp2);
                    }
                }
            }
        }
    }
    
    public static void updateJacobianNonlinear(DenseMatrix64F x, DenseMatrix64F xPlus,
            DenseMatrix64F f, DenseMatrix64F fPlus, Options solverOptions) {
        s = new DenseMatrix64F(x.numRows, 1);
        y = new DenseMatrix64F(x.numRows, 1);
        /* assign s and y*/
        CommonOps.subtract(xPlus, x, s);
        CommonOps.subtract(fPlus, f, y);
        /*check if we should skip the update */
        boolean skipUpdate = true;
        /*t= Rs*/
        DenseMatrix64F t = new DenseMatrix64F(x.numRows, 1);
        for (int i = 0; i < x.numRows; i++) {
            t.set(i, 0, 0.0);
            for (int j = i; j < x.numRows; j++) {
                t.set(i, 0, t.get(i, 0) + l.get(i, j) * s.get(j, 0));
            }
        }
        /*w = S*y-H^T t */
        DenseMatrix64F w = new DenseMatrix64F(x.numRows, 1);
        for (int i = 0; i < x.numRows; i++) {
            w.set(i, 0, solverOptions.getTypicalF().get(i,0) * (fPlus.get(i) - f.get(i)));
            for (int j = 0; j < x.numRows; j++) {
                w.set(i, 0, w.get(i, 0) - h.get(j, i) * t.get(j, 0));
            }
            double tolerance;
            if (solverOptions.isAnalyticalGradient()) {
                tolerance = solverOptions.getMachineEpsilon();
            } else {
                tolerance = Math.sqrt(solverOptions.getMachineEpsilon());
            }
            if (Math.abs(w.get(i, 0)) >= tolerance * solverOptions.getTypicalF().get(i,0) * (Math.abs(fPlus.get(i)) + Math.abs(f.get(i)))) {
                skipUpdate = false;
            } else {
                w.set(i, 0, 0.0);
            }
        }
        /* skip the update or update the h */
        if (!skipUpdate) {
            /*t= Hw*/
            for (int i = 0; i < x.numRows; i++) {
                t.set(i, 0, 0.0);
                for (int j = 0; j < x.numRows; j++) {
                    t.set(i, 0, t.get(i, 0) + BFGS.h.get(i, j) * w.get(j, 0));
                }
            }
            double denom = 0.0;
            for (int i = 0; i < x.numRows; i++) {
                denom += Math.pow(s.get(i) * solverOptions.getTypicalX().get(i,0), 2);
            }
            for (int i = 0; i < x.numRows; i++) {
                s.set(i, 0, s.get(i, 0) * solverOptions.getTypicalX().get(i,0) * solverOptions.getTypicalX().get(i,0) / denom);
            }
            BFGS.qrUpdate(t, s, true);
        }
    }
    
    public static DenseMatrix64F newtonStep(DenseMatrix64F g) {
        DenseMatrix64F ss = new DenseMatrix64F(g.numRows, 1);
        /* Solve Ly=g */
        ss.set(0, 0, -g.get(0, 0) / l.get(0, 0));
        for (int i = 1; i < g.numRows; i++) {
            double temp = 0;
            for (int j = 0; j <= i - 1; j++) {
                temp += l.get(i, j) * ss.get(j, 0);
            }
            ss.set(i, 0, (-g.get(i, 0) - temp) / l.get(i, i));
        }
        /* solve L^Ts=y */
        ss.set(g.numRows - 1, 0, ss.get(g.numRows - 1, 0) / l.get(g.numRows - 1, g.numRows - 1));
        for (int i = g.numRows - 2; i >= 0; i--) {
            double temp = 0;
            for (int j = i + 1; j < g.numRows; j++) {
                temp += l.get(j, i) * ss.get(j, 0);
            }
            ss.set(i, 0, (ss.get(i, 0) - temp) / l.get(i, i));
        }
        return ss;
    }
}
