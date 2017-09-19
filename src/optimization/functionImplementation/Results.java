/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package optimization.functionImplementation;

import java.util.ArrayList;
import java.util.List;

/**
 *
 * @author O. Anil Turkkan <turkkan.1@osu.edu> - anilturkkan.com
 */
public class Results {

    private final List<Double> functionNorm;
    private final List<List<Double>> gradient;
    private final List<List<Double>> x;
    private final List<Double> trustRadius;
    private int functionEvaluations;
    private boolean successful;

    public Results() {
        this.functionNorm = new ArrayList<>();
        this.gradient = new ArrayList<>();
        this.x = new ArrayList<>();
        this.trustRadius = new ArrayList<>();
        this.functionEvaluations=0;
    }

    public void update(double functionNorm, double[] gradient, double[] x, double trustRadius) {
        this.functionNorm.add(functionNorm);
        this.gradient.add(new ArrayList<>());
        for (double g : gradient) {
            this.gradient.get(this.gradient.size() - 1).add(g);
        }
        this.x.add(new ArrayList<>());
        for (double xx : x) {
            this.x.get(this.x.size() - 1).add(xx);
        }
        if (trustRadius != -1.0 ){
            this.trustRadius.add(trustRadius);
        }
    }
    
    @Override
    public String toString(){
        if (this.successful){
            return "Solver successful after " + String.valueOf(x.size()) + " iterations" + System.lineSeparator()+ "Function Evaluations:" + String.valueOf(this.functionEvaluations)+ System.lineSeparator();
        } else {
            return "Solver successful after " + String.valueOf(x.size()) + " iterations" + System.lineSeparator()+ "Function Evaluations:" + String.valueOf(this.functionEvaluations)+ System.lineSeparator();
        }
    }
    
    public void updateFunctionEvaluations(){
        this.functionEvaluations++;
    }

    public List<Double> getFunctionNorm() {
        return functionNorm;
    }

    public List<List<Double>> getGradient() {
        return gradient;
    }

    public List<List<Double>> getX() {
        return x;
    }

    public List<Double> getTrustRadius() {
        return trustRadius;
    }

    public boolean isSuccessful() {
        return successful;
    }

    public void setSuccessful(boolean successful) {
        this.successful = successful;
    }
    
}
