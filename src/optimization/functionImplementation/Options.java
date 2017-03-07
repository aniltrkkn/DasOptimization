/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package optimization.functionImplementation;

import org.ejml.data.DenseMatrix64F;
import org.ejml.ops.CommonOps;
import org.ejml.ops.NormOps;



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
    private boolean BFGSHessian;
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
    private final DenseMatrix64F typicalX;
    private final DenseMatrix64F typicalF;
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
        typicalX = new DenseMatrix64F(n,1);
        typicalF =  new DenseMatrix64F(n,1);
        CommonOps.fill(typicalX,1.0);
        CommonOps.fill(typicalF,1.0);
        //maximum iterations
        maxIterations = 100;
        //max step
        maxStep=NormOps.fastNormP2(typicalX)*1000;
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

    public boolean isBFGSHessian() {
        return BFGSHessian;
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

    public DenseMatrix64F getTypicalX() {
        return typicalX;
    }

    public DenseMatrix64F getTypicalF() {
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
        this.BFGSHessian=!analyticalHessian;
    }

    public void setBFGSHessian(boolean BFGSHessian) {
        this.BFGSHessian = BFGSHessian;
        this.analyticalHessian=!BFGSHessian;
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

    public void setTypicalX(DenseMatrix64F typicalX) {
        assert typicalX.getNumRows()== n;
        assert typicalX.getNumCols()== 1;
        this.typicalX.set(typicalX.copy());
        maxStep=Math.max(NormOps.fastNormP2(typicalX)*1000,maxStep);
    }

    public void setTypicalF(DenseMatrix64F typicalF) {
        assert typicalF.getNumRows()== n;
        assert typicalF.getNumCols()== 1;
        this.typicalF.set(typicalF.copy());
    }

    public void setMaxIterations(int maxIterations) {
        assert maxIterations > 0;
        this.maxIterations = maxIterations;
    }

}
