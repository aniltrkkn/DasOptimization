/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package optimization;

import optimization.functionImplementation.Options;
import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.Math.pow;
import static java.lang.Math.sqrt;
import org.ejml.data.DMatrixRMaj;
import optimization.functionImplementation.ObjectiveFunctionNonLinear;

/**
 *
 * @author anill
 */
public class DoglegTrustSolver {

    //Options
    private final Options options;
    //jacobian matrix
    private final double[][] jac;
    /*
    *run time arguments
    */
    //machine epsilon
    private double mechineEpsilon;
    
    
    
    private double[] xC;
    private double[] xP;
    private double[] fvC;
    private double[] fvP;
    private double[] fvPrev;
    private double fC;
    private double fP;
    private double fPrev;
    private double iterLimit;
    
    private int n;
    private double delta = 0.0;
    private int iteration;
    private double[] g;
    private int termCode;
    
    private double minTol;
    private double stepTol;
    private ObjectiveFunctionNonLinear function;
    
    private double[] sx;
    private double[] sf;
    private double etta;
    private double functionTolerance;
    private int consecmax;
    private double[][] m;
    private double[] m1;
    private double[] m2;
    private double[][] h;
    private double[][] l;
    private double[] sn;
    private boolean firstDog;
    private double cauchyLen;
    private double[] s;
    private double[] ssd;
    private double[] v;
    private double[] xPrev;
    private double maxStep;
    private boolean maxTaken;
    
    /**
     * Constructor for the nonlinear solver
     * @param function a class that contains objective function (gradient, hessian optional)
     * @param options options class
     * @param initialGuess initial guess
     */
    public DoglegTrustSolver(ObjectiveFunctionNonLinear function, Options options, double[] initialGuess) {
        this.function = function;
        this.options=options;
        this.n = initialGuess.length;
        this.xC = new double[n];
        xP = new double[n];
        System.arraycopy(initialGuess, 0, this.xC, 0, n);
        fvC = new double[n];
        fvP = new double[n];
        jac = new double[n][n];
        g = new double[n];
        m = new double[n][n];
        m1 = new double[n];
        m2 = new double[n];
        h = new double[n][n];
        l = new double[n][n];
        sn = new double[n];
        s = new double[n];
        ssd = new double[n];
        v = new double[n];
        xPrev = new double[n];
        fvPrev = new double[n];
    }
    /**
     * set the machine epsilon
     */
    private void setMechineEpsilon() {
          mechineEpsilon = 1.0;
        while ((mechineEpsilon / 2.0 + 1.0) != 1.0) {
            mechineEpsilon /= 2.0;
        }
        mechineEpsilon = mechineEpsilon * 2.0;
    }

    private int validateInput() {
        //validate the input
        if (delta == 0.0) {
            delta = -1.0;
        }
        //typical x
        sx = new double[n];
        for (int i = 0; i < n; i++) {
                sx[i] = 1.0;
            }

        //typical f
        sf = new double[n];
        for (int i = 0; i < n; i++) {
                sf[i] = 1.0;
            }
        //function tolerance
        functionTolerance = pow(mechineEpsilon, 2.0 / 3.0);
        //step tolerance
        stepTol = pow(mechineEpsilon, 2.0 / 3.0);
        //max step
        double dummy=0.0;
        for(int i=0;i<n;i++){
            dummy+=sx[i]*xC[i]*sx[i]*xC[i];
        }
        dummy = sqrt(dummy);
        double dummy1=0.0;
        for(int i=0;i<n;i++){
            dummy1+=sx[i]*sx[i];
        }
        dummy1 = sqrt(dummy1);
        maxStep=1000*max(dummy,dummy1);
        iterLimit=10000;
        //min tolerance
        minTol = pow(mechineEpsilon, 2.0 / 3.0);
        if (n < 1) {
            return -1;
        } else {
            return 0;
        }

    }

    private double[] fillFx(double[] x) {
        //find the function values
        DMatrixRMaj xx= new DMatrixRMaj(n,1);
        for(int i=0;i<n;i++){
            xx.set(i,0,x[i]);
        }
        return function.getF(xx).getData();
    }

    private double sumOfSquares() {
        //calculate the sum of squares 
        fvP = fillFx(xC);
        double sum = 0.0;
        for (int i = 0; i < n; i++) {
            sum += 0.5 * ((sf[i] * fvP[i]) * (sf[i] * fvP[i]));
        }
        return sum;
    }

    private boolean stop0() {
        //stopping cnditions for iteration 0
        consecmax = 0;
        double dummy=Double.MIN_VALUE;
        for (int i = 0; i < n; i++) {
            dummy = max(dummy,abs(sf[i] * fvP[i]));
        }
        return false;
        //return dummy <= 0.01 * functionTolerance;
    }

    private void jacobian(double [] x,double[] f) {
        //finite difference jacobian approximation
        double[] dummyFx;
        double sqrDelta = sqrt(mechineEpsilon);
        for (int j = 0; j < n; j++) {
            //calculate column j of J
            double sign = getSign(x[j]);
            double stepSize = sqrDelta * max(abs(x[j]), 1 / sx[j]) * sign;
            double temp = x[j];
            x[j] += stepSize;
            //to reduce finite difference precision errors
            stepSize = x[j] - temp;
            dummyFx = fillFx(x);
            for (int i = 0; i < n; i++) {
                jac[i][j] = (dummyFx[i] - f[i]) / stepSize;
            }
            //reset the value
            x[j] = temp;
        }
    }

    private double getSign(double x) {
        //return sign of a number
        double sign = Math.signum(x);
        if (sign == 0.0) {
            sign = 1.0;
        }
        return sign;
    }

    private boolean qr() {
        boolean singular = false;
        for (int k = 0; k < n - 1; k++) {
            double dummy = Double.MIN_VALUE;
            for (int i = k; i < n; i++) {
                dummy = max(dummy,abs(m[i][k]));
            }
            if (dummy == 0.0) {
                //matrix is singular
                m1[k] = 0.0;
                m2[k] = 0.0;
                singular = true;
            } else {
                //form Q and premultiply M by it
                for (int i = k; i < n; i++) {
                    m[i][k] = m[i][k] / dummy;
                }
                double sigma = 0.0;
                for (int i = k; i < n; i++) {
                    sigma += m[i][k] * m[i][k];
                }
                sigma = sqrt(sigma) * getSign(m[k][k]);
                m[k][k] = m[k][k] + sigma;
                m1[k] = sigma * m[k][k];
                m2[k] = -dummy * sigma;
                for (int j = k + 1; j < n; j++) {
                    double dummy1 = 0.0;
                    for (int i = k; i < n; i++) {
                        dummy1 += m[i][k] * m[i][j];
                    }
                    dummy1 /= m1[k];
                    for (int i = k; i < n; i++) {
                        m[i][j] = m[i][j] - dummy1 * m[i][k];
                    }
                }
            }
        }
        if (m[n - 1][n - 1] == 0.0) {
            singular = true;
        }
        m2[n - 1] = m[n - 1][n - 1];
        return singular;
    }

    private double[] rSolve(double[] b) {
        //solve Rx=b
        b[n - 1] = b[n - 1] / m2[n - 1];
        for (int i = n - 2; i >= 0; i--) {
            double dummy = 0.0;
            for (int j = i + 1; j < n; j++) {
                dummy += m[i][j] * b[j];
            }
            b[i] = (b[i] - dummy) / m2[i];
        }
        return b;
    }

    private double conditionNumber() {
        //estimate the condition number of a upper triangular matrix
        double[] p = new double[n];
        double[] pm = new double[n];
        double[] x = new double[n];
        double estimate = abs(m2[0]);
        for (int j = 1; j < n; j++) {
            double temp = 0.0;
            for (int i = 0; i <= j - 1; i++) {
                temp += abs(m[i][j]);
            }
            temp += abs(m2[j]);
            estimate = max(estimate, temp);
        }
        //solve Rx=e
        x[0] = 1.0 / m2[0];
        for (int i = 1; i < n; i++) {
            p[i] = m[0][i] * x[0];
        }
        for (int j = 1; j < n; j++) {
            //select e and calculate x[j]
            double xp = (1 - p[j]) / m2[j];
            double xm = (-1 - p[j]) / m2[j];
            double temp = abs(xp);
            double tempM = abs(xm);
            for (int i = j + 1; i< n; i++) {
                pm[i] = p[i] + m[j][i] * xm;
                tempM += abs(pm[i]) / abs(m2[i]);
                p[i] += m[j][i] * xp;
                temp += abs(p[i]) / abs(m2[i]);
            }
            if (temp >= tempM) {
                x[j] = xp;
            } else {
                x[j] = xm;
                System.arraycopy(pm, j + 1, p, j + 1, n - (j + 1));
            }
        }
        double xNorm = 0.0;
        for (int j = 0; j < n; j++) {
            xNorm += abs(x[j]);
        }
        estimate /= xNorm;
        //calculate R-1x
        x = rSolve(x);
        xNorm = 0.0;
        for (int j = 0; j < n; j++) {
            xNorm += abs(x[j]);
        }
        estimate *= xNorm;
        return estimate;
    }

    private double cholDecomp(double maxOfFl) {
        //perturbed cholesky decomposition
        double minL = pow(mechineEpsilon, 0.25) * maxOfFl;
        if (maxOfFl == 0.0) {
            //h is known to be positive definite
            maxOfFl = Double.MIN_VALUE;
            for (int i = 0; i < n; i++) {
                maxOfFl=max(maxOfFl,h[i][i]);
            }
            maxOfFl = sqrt(maxOfFl);
        }
        double minL2 = sqrt(mechineEpsilon) * maxOfFl;
        double maxAdd = 0.0;
        for (int j = 0; j < n; j++) {
            //form column j of L
            double dummy = 0.0;
            for (int i = 0; i <= j - 1; i++) {
                dummy += l[j][i] * l[j][i];
            }
            l[j][j] = h[j][j] - dummy;
            double minljj = 0.0;
            for (int i = j + 1; i < n; i++) {
                double dummy2 = 0.0;
                for (int k = 0; k <= j - 1; k++) {
                    dummy2 += l[i][k] * l[j][k];
                }
                l[i][j] = h[j][i] - dummy2;
                minljj = max(minljj, abs(l[i][j]));
            }
            minljj = max(minljj / maxOfFl, minL);
            if (l[j][j] > (minljj * minljj)) {
                //normal cholesky iteration
                l[j][j] = sqrt(l[j][j]);
            } else {
                //augment hjj
                if (minljj < minL2) {
                    //only possible when input =0
                    minljj = minL2;
                }
                maxAdd = max(maxAdd, minljj * minljj - l[j][j]);
                l[j][j] = minljj;
            }
            for (int i = j + 1; i < n; i++) {
                l[i][j] = l[i][j] / l[j][j];
            }
        }
        return maxAdd;
    }

    private double[] lSolve(double[] b) {
        //solve Ly=b for y
        double[] y = new double[n];
        y[0] = b[0] / l[0][0];
        for (int i = 1; i < n; i++) {
            double dummy = 0.0;
            for (int j = 0; j <= i - 1; j++) {
                dummy += l[i][j] * y[j];
            }
            y[i] = (b[i] - dummy) / l[i][i];
        }
        return y;
    }

    private double[] lTSolve(double[] y) {
        //solve Ly=b for y
        double[] x = new double[n];
        x[n - 1] = y[n - 1] / l[n - 1][n - 1];
        for (int i = n - 2; i >= 0; i--) {
            double dummy = 0.0;
            for (int j = i + 1; j < n; j++) {
                dummy += l[j][i] * x[j];
            }
            x[i] = (y[i] - dummy) / l[i][i];
        }
        return x;
    }

    private double[] cholSolve() {
        //cholesky solve driver
        //solve ly=g
        double[] s = lSolve(g);
        //solve Lt.s=y
        s = lTSolve(s);
        for (int i = 0; i < n; i++) {
            s[i] = -s[i];
        }
        return s;
    }

    private double[] qrSolve(double[] b) {
        //solve (QR)x = b for x
        //b=QTb
        for (int j = 0; j < n - 1; j++) {
            // b=Qjb
            double dummy = 0.0;
            for (int i = j; i < n; i++) {
                dummy += m[i][j] * b[i];
            }
            double sigma = dummy / m1[j];
            for (int i = j; i < n; i++) {
                b[i] -= sigma * m[i][j];
            }
        }
        //b=R-1b
        return rSolve(b);
    }

    private void neModel() {
        //formulation of the affine model
        //condition the jacobian
        //calculate the newton step

        //m=DfJ
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                m[i][j] = sf[i] * jac[i][j];
            }
        }
        //qr factorization
        boolean singular = qr();
        double estimate;
        if (!singular) {
            //estimate K1(RD^-1)
            for (int j = 0; j < n; j++) {
                for (int i = 0; i <= j - 1; i++) {
                    m[i][j] = m[i][j] / sx[j];
                }
                m2[j] = m2[j] / sx[j];
            }
            estimate = conditionNumber();
        } else {
            estimate = 0.0;
        }
        if (singular || estimate > 1 / sqrt(mechineEpsilon)) {
            //perturb jacobian 
            //H =JtDf^2J

            for (int i = 0; i < n; i++) {
                for (int j = i; j < n; j++) {
                    h[i][j] = 0.0;
                    for (int k = 0; k < n; k++) {
                        h[i][j] += jac[k][i] * jac[k][j] * sf[k] * sf[k];
                    }
                }
            }
            //Hnorm || Dx-1HDx-1||
            double hNorm = 0.0;
            for (int j = 0; j < n; j++) {
                hNorm = abs(h[0][j]) / sx[j];
            }
            hNorm = (1 / sx[0]) * hNorm;
            for (int i = 1; i < n; i++) {
                double temp1 = 0.0;
                for (int j = 0; j <= i; j++) {
                    temp1 = abs(h[j][i]) / sx[j];
                }
                double temp2 = 0.0;
                for (int j = i + 1; j < n; j++) {
                    temp2 = abs(h[i][j]) / sx[j];
                }
                double temp = (1 / sx[i]) * (temp1 + temp2);
                hNorm = max(temp, hNorm);
            }
            // H = H+ sqrt(n*eps)*Hnorm*Dx^2
            for (int i = 0; i < n; i++) {
                h[i][i] += sqrt(n * mechineEpsilon) * hNorm * sx[i] * sx[i];
            }
            cholDecomp(0.0);
            sn = cholSolve();
        } else {
            //calculate normal Newton step
            //reset R=RDx
            for (int j = 0; j < n; j++) {
                for (int i = 0; i <= j - 1; i++) {
                    m[i][j] *= sx[j];
                }
                m2[j] *= sx[j];
            }
            for (int i = 0; i < n; i++) {
                sn[i] = -sf[i] * fvC[i];
            }
            sn = qrSolve(sn);
            //lower triangle of m=RT
            for (int i = 0; i < n; i++) {
                m[i][i] = m2[i];
                for (int j = 0; j <= i - 1; j++) {
                    m[i][j] = m[j][i];
                }
            }
        }

    }

    public boolean dogStep(double newtonLength) {
        //find approximate solution to
        //minimize gT*s+1/2sTLLTs
        boolean newTaken;
        if (newtonLength <= delta) {
            //s is newton step
            newTaken = true;
            System.arraycopy(sn, 0, s, 0, n);
            delta = newtonLength;
        } else {
            //newton step too long ,s on double dogleg curve
            newTaken = false;
            if (firstDog) {
                //calculate double dogle curve
                firstDog = false;

                double alpha = 0.0;
                for (int i = 0; i < n; i++) {
                    alpha += g[i] / sx[i] * g[i] / sx[i];
                }
                //implement b = || LtDx-2g||
                double bmechineEpsilon = 0.0;
                for (int i = 0; i < n; i++) {
                    double temp = 0.0;
                    for (int j = i; j < n; j++) {
                        temp += m[j][i] * g[j] / (sx[j] * sx[j]);
                    }
                    bmechineEpsilon += temp * temp;
                }
                for (int i = 0; i < n; i++) {
                    ssd[i] = -(alpha / bmechineEpsilon) / sx[i] * g[i];
                }
                cauchyLen = alpha * sqrt(alpha) / bmechineEpsilon;
                double dummy = 0.0;
                for (int i = 0; i < n; i++) {
                    dummy += abs(g[i] * sn[i]);
                }
                etta = 0.2 + (0.8 * alpha * alpha / (bmechineEpsilon * dummy));
                for (int i = 0; i < n; i++) {
                    v[i] = etta * sx[i] * sn[i] - ssd[i];
                }
                if(delta == -1.0){
                    //no user delta
                    delta=min(cauchyLen,maxStep);
                }  
            }
            //calculate double dogleg curve
            if(etta*newtonLength<=delta){
                //take partial step in newton direction
                for (int i = 0; i < n; i++) {
                    s[i] = delta/newtonLength*sn[i];
                }
            } else if(cauchyLen >= delta){
                //take steepest direction
                for (int i = 0; i < n; i++) {
                    s[i] = delta/cauchyLen*ssd[i]/sx[i];
                }
            } else {
                //calculate  convex combination2
                double temp=0.0;
                for (int i = 0; i < n; i++) {
                    temp += ssd[i]*ssd[i];
                }
                double tempV=0.0;
                for (int i = 0; i < n; i++) {
                    tempV += v[i]*v[i];
                }
                double lambda=(-temp+sqrt(temp*temp-tempV*(cauchyLen*cauchyLen-delta*delta)))/tempV;
                for (int i = 0; i < n; i++) {
                    s[i] = (ssd[i]+lambda*v[i])/sx[i];
                } 
               /* double a=0.0;
                for (int i = 0; i < n; i++) {
                    a += ssd[i]*ssd[i];
                }
                double b=0.0;
                for (int i = 0; i < n; i++) {
                   b += sn[i]*sn[i];
                }
                double c=0.0;
                for (int i = 0; i < n; i++) {
                   c += (sn[i]-ssd[i])*(sn[i]-ssd[i]);
                }
                double d=(a+b-c)/2.0;
                double alpha=(b-delta*delta)/(b-d+sqrt(d*d-a*b+delta*delta*c));
                for (int i = 0; i < n; i++) {
                    s[i] = (alpha*ssd[i]+(1-alpha)*sn[i])/sx[i];
                }*/
            }
        }
        return newTaken;
    }
    
    public double functionNorm(double[] f){
        //the minimazation function for nonlinear functions
        double norm=0.0;
        for(int i=0;i<n;i++){
            norm += 0.5 * (( f[i] * sf[i]) * ( f[i] * sf[i]));
        }
        return norm;
    }
    public int trustRegup(boolean newtonTaken,int code){
        //given a step s produced by dogleg algorithm decide whether x+ is
        //accapted or not
        int retCode=0;
        maxTaken=false;
        double alpha= 1e-4;
        double stepLength=0.0;
        for(int i=0;i<n;i++){
            stepLength += sx[i]*sx[i]*s[i]*s[i];
        }
        stepLength=sqrt(stepLength);
        for(int i=0;i<n;i++){
            xP[i]=xC[i]+s[i];
        }
        fvP=fillFx(xP);
        fP=functionNorm(fvP);
        fC=functionNorm(fvC);
        double deltaF=fP-fC;
        double initSlope=0.0;
        for(int i=0;i<n;i++){
            initSlope += g[i]*s[i];
        }
        if(code !=3){
            fPrev=0.0;
        }
        if(code ==3 && (fP >= fPrev  || deltaF > alpha*initSlope )){
            //reset xp to xpre and terminate globa step
            retCode=0;
            System.arraycopy(xPrev, 0, xP, 0, n);
            fP=fPrev;
            System.arraycopy(fvPrev, 0, fvP, 0, n);
            delta /= 2.0;
        } else if (deltaF >= alpha*initSlope){
            //f(xP) too large
            double relLength=Double.MIN_VALUE;
            for(int i=0;i<n;i++){
                relLength=max(relLength,abs(s[i])/max(abs(xP[i]),1/sx[i]) );
            }
            if (relLength <stepTol){
                //xp-xc too small, terminate
                retCode=1;
                System.arraycopy(xC, 0, xP, 0, n);
            } else {
                //reduce delta, continue
                retCode=2;
                double deltaTemp=-(initSlope*stepLength)/(2*(deltaF-initSlope));
                if (deltaTemp<0.1*delta){
                    delta=0.1*delta;
                } else if (deltaTemp>0.5*delta){
                    delta=0.5*delta;
                } else {
                    delta= deltaTemp;
                }
            }
        } else {
            //f(xP) sufficiently small
            double deltaFPred = initSlope;
            for(int i=0;i<n;i++){
                double temp=0.0;
                for(int j=i;j<n;j++){
                    temp += m[j][i]*s[j];
                }
                deltaFPred += (temp*temp/2);
            }
            if(retCode != 2 && (abs(deltaFPred-deltaF)<= abs(deltaF)*0.1 || deltaF <= initSlope) && newtonTaken == false && delta<= 0.99*maxStep ){
                //double delta and continue
                retCode=3;
                System.arraycopy(xP, 0, xPrev, 0, n);
                fPrev=fP;
                System.arraycopy(fvP, 0, fvPrev, 0, n);
                delta= min(2*delta,maxStep);
            } else {
                //accept xP as the new iterate
                retCode=0;
                if (stepLength >0.99*maxStep){
                    maxTaken=true;
                } 
                if (deltaF >= 0.1*deltaFPred){
                    delta /=2;
                } else if (deltaF <= 0.75*deltaFPred){
                    delta = min(2*delta,maxStep);
                }
            }
        }
        return retCode;
    }
    
    public int dogDriver() {
        //find x on the double dogleg curve
        //f(x)<= f(x)+alpha*gt*(x+-xc)
        int retCode = 4;
        firstDog = true;
        double newtonLength = 0.0;
        for (int i = 0; i < n; i++) {
            newtonLength += sx[i] * sn[i] * sx[i] * sn[i];
        }
        newtonLength = sqrt(newtonLength);
        while(retCode >=2){
            boolean newtonTaken=dogStep(newtonLength);
            retCode= trustRegup(newtonTaken, retCode);
        }
        return retCode;
    }
    
    private void neStop(int retCode){
        //termination criteria
        termCode=0;
        double dummy= Double.MIN_VALUE;
        for(int i=0;i<n;i++){
            dummy=max(dummy,sf[i]*abs(fvP[i]));
        }
        double dummy1= Double.MIN_VALUE;
        for(int i=0;i<n;i++){
            dummy1=max(dummy1,abs(xP[i]-xC[i])/max(abs(xP[i]),1/sx[i] ));
        }
        double dummy2= Double.MIN_VALUE;
        for(int i=0;i<n;i++){
            dummy2=max(dummy2,abs(g[i])*max(abs(xP[i]),1/sx[i] )/(max(fP,n/2)) );
        }
        
        if(retCode == 1){
            termCode = 3;
        } else if (dummy<=functionTolerance){
            termCode =1;
        } else if(dummy1<= stepTol){
            termCode=2;
        } else if (iteration>= iterLimit){
            termCode=4;
        } else if (maxTaken){
            consecmax += 1;
            if (consecmax == 5){
                termCode=5;
            }
        } else {
            consecmax=0;
            if(dummy2<=minTol){
                termCode=6;
            }
        }
    }
    
    public void solve() {
        //solve given nonlinear equations   
        setMechineEpsilon();
        if (validateInput() < 0) {
            //errors found
            return;
        }
        //iteration 0
        iteration = 0;
        fP=sumOfSquares();
        boolean terminate = stop0();
        if (terminate) {
            //initial point is the solution
            return;
        } else {
            jacobian(xC,fvP);
            //Jt*Df^2*F
            for (int i = 0; i < n; i++) {
                g[i] = 0.0;
                for (int j = 0; j < n; j++) {
                    g[i] += jac[j][i] * fvP[j] * sf[j] * sf[j];
                }
            }
            System.arraycopy(fvP,0,fvC,0,n);
        }
        boolean restart = true;
        //iteration section
        while (!terminate) {
            iteration += 1;
            neModel();
            int retCode =dogDriver();
            jacobian(xP,fvP);
            //Jt*Df^2*F
            for (int i = 0; i < n; i++) {
                g[i] = 0.0;
                for (int j = 0; j < n; j++) {
                    g[i] += jac[j][i] * fvP[j] * sf[j] * sf[j];
                }
            }
            neStop(retCode);
            if (termCode >0){
                terminate=true;
            } else {
                restart =false;
            }
            System.arraycopy(xP, 0, xC, 0, n);
            fC=fP;
            System.arraycopy(fvP, 0, fvC, 0, n);
            System.out.println("iteration");
            System.out.println(iteration);
            System.out.println("x");
            for(int i=0;i<xC.length;i++){
                System.out.println(xC[i]);
            }
            System.out.println("g");
            for(int i=0;i<xC.length;i++){
                System.out.println(g[i]);
            }
            System.out.println("f");
            System.out.println(fC);
            System.out.println("delta");
            System.out.println(delta);
            System.out.println("code");
            System.out.println(termCode);
        }
    }
}
