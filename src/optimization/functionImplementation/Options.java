/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package optimization.functionImplementation;

import cern.colt.matrix.DoubleMatrix1D;
import cern.colt.matrix.impl.DenseDoubleMatrix1D;

/**
 *
 * @author anill
 */
public class Options {

    private final int n;
    /* algorithm options */
    public static final int LINE_SEARCH = 0;
    public static final int DOGLEG_TRUST_REGION = 1;
    private int algorithm;
    /* analytical or finite difference gradients */
    private boolean analyticalGradient;
    private boolean analyticalHessian;
    /* machine epsilon */
    private double machineEpsilon;
    /* tolerances */
    private double gradientTolerance;
    private double stepTolerance;
    private double functionTolerance;
    private double minTolerance;
    private double trussRegionRadius;
    private boolean trussRegionRadiusSet;
    /* typical values */
    private DoubleMatrix1D typicalX;
    private DoubleMatrix1D typicalF;
    /* maximum iteration*/
    private int maxIterations;
    /*maximum step*/
    private double maxStep;

    public Options(int n) {
        this.n = n;
        //make sure that n is greater than 0 
        assert n > 0;
        //default algorithm is DOGLEG_TRUST_REGION
        algorithm = Options.DOGLEG_TRUST_REGION;
        //default gradients are calculated with finite difference
        analyticalGradient = false;
        analyticalHessian = false;
        //calculate tolerances
        machineEpsilon = 1.0;
        while ((machineEpsilon / 2.0 + 1.0) != 1.0) {
            machineEpsilon /= 2.0;
        }
        machineEpsilon = machineEpsilon * 2.0;
        gradientTolerance = Math.pow(machineEpsilon, 1.0 / 3.0);
        stepTolerance = Math.pow(machineEpsilon, 2.0 / 3.0);
        functionTolerance = Math.pow(machineEpsilon, 1.0 / 3.0);
        minTolerance = Math.pow(machineEpsilon, 2.0 / 3.0);
        trussRegionRadiusSet = false;
        trussRegionRadius = -1.0;
        //typical values
        typicalX = new DenseDoubleMatrix1D(n);
        typicalF =  new DenseDoubleMatrix1D(n);
        typicalX.assign(1.0);
        typicalF.assign(1.0);
        //maximum iterations
        maxIterations = 100;
    }

    public Options(Options solverOptions) {
        n = solverOptions.getN();
        algorithm = solverOptions.getAlgorithm();
        analyticalGradient = solverOptions.isAnalyticalGradient();
        analyticalHessian = solverOptions.isAnalyticalHessian();
        machineEpsilon = solverOptions.getMachineEpsilon();
        gradientTolerance= solverOptions.getGradientTolerance();
        stepTolerance=solverOptions.getStepTolerance();
        trussRegionRadius = solverOptions.getTrussRegionRadius();
        trussRegionRadiusSet=solverOptions.isTrussRegionRadiusSet();
        typicalX=solverOptions.getTypicalX().copy();
        typicalF=solverOptions.getTypicalF().copy();
        maxIterations=solverOptions.getMaxIterations();
        maxStep=solverOptions.getMaxStep();
    }

    public double getMachineEpsilon() {
        return machineEpsilon;
    }

    public int getN() {
        return n;
    }

    public int getAlgorithm() {
        return algorithm;
    }

    public boolean isAnalyticalGradient() {
        return analyticalGradient;
    }

    public boolean isAnalyticalHessian() {
        return analyticalHessian;
    }

    public double getGradientTolerance() {
        return gradientTolerance;
    }

    public double getStepTolerance() {
        return stepTolerance;
    }

    public double getFunctionTolerance() {
        return functionTolerance;
    }

    public double getMinTolerance() {
        return minTolerance;
    }
    
    public double getTrussRegionRadius() {
        return trussRegionRadius;
    }

    public DoubleMatrix1D getTypicalX() {
        return typicalX;
    }

    public DoubleMatrix1D getTypicalF() {
        return typicalF;
    }

    public boolean isTrussRegionRadiusSet() {
        return trussRegionRadiusSet;
    }

    public int getMaxIterations() {
        return maxIterations;
    }

    public double getMaxStep() {
        return maxStep;
    }

    public void setAlgorithm(int algorithm) {
        assert algorithm >= LINE_SEARCH;
        assert algorithm <= DOGLEG_TRUST_REGION;
        this.algorithm = algorithm;
    }

    public void setAnalyticalGradient(boolean analyticalGradient) {
        this.analyticalGradient = analyticalGradient;
    }

    public void setAnalyticalHessian(boolean analyticalHessian) {
        this.analyticalHessian = analyticalHessian;
    }

    public void setGradientTolerance(double gradientTolerance) {
        this.gradientTolerance = gradientTolerance;
    }

    public void setStepTolerance(double stepTolerance) {
        this.stepTolerance = stepTolerance;
    }

    public void setFunctionTolerance(double functionTolerance) {
        this.functionTolerance = functionTolerance;
    }

    public void setMinTolerance(double minTolerance) {
        this.minTolerance = minTolerance;
    }
    
    public void setTrussRegionRadius(double trussRegionRadius) {
        trussRegionRadiusSet = true;
        this.trussRegionRadius = trussRegionRadius;
    }

    public void setTypicalX(DoubleMatrix1D typicalX) {
        assert typicalX.size() == n;
        this.typicalX = typicalX.copy();
        for(int i=0;i<typicalX.size();i++){
            this.typicalX.set(i, 1.0/this.typicalX.get(i));
        }
    }

    public void setTypicalF(DoubleMatrix1D typicalF) {
        assert typicalF.size() == n;
        this.typicalF = typicalF.copy();
        for(int i=0;i<typicalF.size();i++){
            this.typicalF.set(i, 1.0/this.typicalF.get(i));
        }
    }

    public void setMaxIterations(int maxIterations) {
        assert maxIterations > 0;
        this.maxIterations = maxIterations;
    }

    public void setMaxStep(DoubleMatrix1D initialGuess) {
        maxStep = Math.max(1.0,Math.max(Math.sqrt(initialGuess.zDotProduct(initialGuess) * typicalX.zDotProduct(typicalX)),
                    Math.sqrt(initialGuess.zDotProduct(initialGuess)))) * 1000;
    }

}
