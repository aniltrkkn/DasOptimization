/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package optimization.functionImplementation;


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
    private boolean centralDifferenceGradient;
    private boolean analyticalJacobian;
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
    /* maximum iteration*/
    private int maxIterations;
    /*maximum step*/
    private double maxStep;
    /*save iteration details*/
    private boolean saveIterationDetails;

    public Options(int n) {
        this.n = n;
        //make sure that n is greater than 0 
        assert n > 0;
        //default algorithm is DOGLEG_TRUST_REGION
        algorithm = Options.LINE_SEARCH;
        //default gradients are calculated with finite difference
        analyticalGradient = false;
        centralDifferenceGradient = false;
        analyticalHessian = false;
        analyticalJacobian = false;
        //calculate tolerances
        machineEpsilon = 1.0;
        while ((machineEpsilon / 2.0 + 1.0) != 1.0) {
            machineEpsilon /= 2.0;
        }
        machineEpsilon = machineEpsilon * 2.0;
        gradientTolerance = 1e-8;
        stepTolerance = 1e-8;
        functionTolerance = 1e-8;
        minTolerance = 1e-8;
        trussRegionRadiusSet = false;
        trussRegionRadius = -1.0;
        //maximum iterations
        maxIterations = 100;
        //max step
        maxStep=1000.0;
        //details
        this.saveIterationDetails=false;
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
        maxIterations=solverOptions.getMaxIterations();
        maxStep=solverOptions.getMaxStep();
        this.saveIterationDetails=solverOptions.isSaveIterationDetails();
    }

    public void setAllTolerances(double tolerance){
        gradientTolerance =tolerance;
        stepTolerance = tolerance;
        functionTolerance = tolerance;
        minTolerance = tolerance;
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

    public boolean isTrussRegionRadiusSet() {
        return trussRegionRadiusSet;
    }

    public int getMaxIterations() {
        return maxIterations;
    }

    public void setMaxStep(double maxStep) {
        this.maxStep = maxStep;
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

    public boolean isCentralDifferenceGradient() {
        return centralDifferenceGradient;
    }

    public void setCentralDifferenceGradient(boolean centralDifferenceGradient) {
        this.centralDifferenceGradient = centralDifferenceGradient;
    }  

    public boolean isAnalyticalJacobian() {
        return analyticalJacobian;
    }

    public void setAnalyticalJacobian(boolean analyticalJacobian) {
        this.analyticalJacobian = analyticalJacobian;
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

    public void setMaxIterations(int maxIterations) {
        assert maxIterations > 0;
        this.maxIterations = maxIterations;
    }

    public boolean isSaveIterationDetails() {
        return saveIterationDetails;
    }

    public void setSaveIterationDetails(boolean saveIterationDetails) {
        this.saveIterationDetails = saveIterationDetails;
    }
}
