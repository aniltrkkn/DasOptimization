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
package optimization.functionImplementation;

import java.util.ArrayList;
import java.util.List;

public class Results {

    private final List<Double> functionNorm;
    private final List<List<Double>> gradient;
    private final List<List<Double>> x;
    private final List<Double> trustRadius;
    private int functionEvaluations;
    private boolean successful;
    private String stopReason;

    public Results() {
        this.functionNorm = new ArrayList<>();
        this.gradient = new ArrayList<>();
        this.x = new ArrayList<>();
        this.trustRadius = new ArrayList<>();
        this.functionEvaluations = 0;
    }

    /**
     * Update the results object with new iteration values
     *
     * @param functionNorm function value after current iteration
     * @param gradient gradient vector
     * @param x x vector
     * @param trustRadius truss radius if applicable
     */
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
        if (trustRadius != -1.0) {
            this.trustRadius.add(trustRadius);
        }
    }

    @Override
    public String toString() {
        if (this.successful) {
            return "Solver successful after " + String.valueOf(x.size()) + " iterations" + System.lineSeparator() + "Function Evaluations:" + String.valueOf(this.functionEvaluations) + System.lineSeparator() + "Stop Reason:" + this.stopReason + System.lineSeparator();
        } else {
            return "Solver failed after " + String.valueOf(x.size()) + " iterations" + System.lineSeparator() + "Function Evaluations:" + String.valueOf(this.functionEvaluations) + System.lineSeparator() + "Stop Reason:" + this.stopReason + System.lineSeparator();
        }
    }

    public void updateFunctionEvaluations() {
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

    public String getStopReason() {
        return stopReason;
    }

    public void setStopReason(String stopReason) {
        this.stopReason = stopReason;
    }

    public int getFunctionEvaluations() {
        return functionEvaluations;
    }
}
