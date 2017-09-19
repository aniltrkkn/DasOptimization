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
package test;

import optimization.functionImplementation.ObjectiveFunctionUnconstrained;
import optimization.functionImplementation.Options;
import org.ejml.data.DMatrixRMaj;
import solvers.UnconstrainedOptimizer;


public class UnconstrainedTest {

    public static void bealeFunction(int algorithm, boolean analyticalGradient, boolean analyticalHessian, double initialMultiplier) {
        //https://www.sfu.ca/~ssurjano/beale.html
        ObjectiveFunctionUnconstrained function = new ObjectiveFunctionUnconstrained() {
            @Override
            public double getF(DMatrixRMaj x) {
                double x1 = x.get(0);
                double x2 = x.get(1);
                return (1.5 - x1 + x1 * x2) * (1.5 - x1 + x1 * x2) + (2.25 - x1 + x1 * x2 * x2) * (2.25 - x1 + x1 * x2 * x2) + (2.625 - x1 + x1 * x2 * x2 * x2) * (2.625 - x1 + x1 * x2 * x2 * x2);
            }

            @Override
            public DMatrixRMaj getG(DMatrixRMaj x) {
                double x1 = x.get(0);
                double x2 = x.get(1);
                DMatrixRMaj g = new DMatrixRMaj(2, 1);
                g.set(0, 0, 2.0 * (x2 * x2 - 1.0) * (x1 * x2 * x2 - x1 + 9.0 / 4.0) + 2.0 * (x2 * x2 * x2 - 1) * (x1 * x2 * x2 * x2 - x1 + 21.0 / 8.0) + 2 * (x2 - 1) * (x1 * x2 - x1 + 3.0 / 2.0));
                g.set(1, 0, 2.0 * x1 * (x1 * x2 - x1 + 3.0 / 2.0) + 4.0 * x1 * x2 * (x1 * x2 * x2 - x1 + 9.0 / 4.0) + 6 * x1 * x2 * x2 * (x1 * x2 * x2 * x2 - x1 + 21.0 / 8.0));
                return g;
            }

            @Override
            public DMatrixRMaj getH(DMatrixRMaj x) {
                double x1 = x.get(0);
                double x2 = x.get(1);
                DMatrixRMaj h = new DMatrixRMaj(2, 2);
                h.set(0, 0, 2 * (x2 - 1.0) * (x2 - 1.0) + 2.0 * (x2 * x2 - 1.0) * (x2 * x2 - 1.0) + 2.0 * (x2 * x2 * x2 - 1.0) * (x2 * x2 * x2 - 1.0));
                h.set(0, 1, 9.0 * x2 - 4.0 * x1 - 4.0 * x1 * x2 - 12.0 * x1 * x2 * x2 + 8.0 * x1 * x2 * x2 * x2 + 12 * x1 * x2 * x2 * x2 * x2 * x2 + (63.0 * x2 * x2) / 4.0 + 3.0);
                h.set(1, 0, 9.0 * x2 - 4.0 * x1 - 4.0 * x1 * x2 - 12.0 * x1 * x2 * x2 + 8.0 * x1 * x2 * x2 * x2 + 12.0 * x1 * x2 * x2 * x2 * x2 * x2 + (63.0 * x2 * x2) / 4.0 + 3.0);
                h.set(1, 1, (x1 * (63.0 * x2 - 4.0 * x1 - 24.0 * x1 * x2 + 24.0 * x1 * x2 * x2 + 60.0 * x1 * x2 * x2 * x2 * x2 + 18.0)) / 2.0);
                return h;
            }
        };
        /* set up initial guess */
        DMatrixRMaj initialGuess = new DMatrixRMaj(2, 1);
        initialGuess.set(0, 0, -4.5 * initialMultiplier);
        initialGuess.set(1, 0, -4.5 * initialMultiplier);
        Options options = new Options(2);
        /* options */
        options.setAlgorithm(algorithm);
        options.setAnalyticalGradient(analyticalGradient);
        options.setAnalyticalHessian(analyticalHessian);
        options.setSaveIterationDetails(true);
        /* start optimizer */
        UnconstrainedOptimizer solver = new UnconstrainedOptimizer(function, options);
        solver.solve(new DMatrixRMaj(initialGuess));
        System.out.println(solver.getResults());
        System.out.println("F: " + solver.getFx());
        System.out.println("x: " + solver.getX());
    }

    public static void helicalValleyFunction(int algorithm, boolean analyticalGradient, boolean analyticalHessian, double initialMultiplier) {
        ObjectiveFunctionUnconstrained function = new ObjectiveFunctionUnconstrained() {
            @Override
            public double getF(DMatrixRMaj x) {
                double x1 = x.get(0);
                double x2 = x.get(1);
                double x3 = x.get(2);
                double theta;
                if (x1 >= 0) {
                    theta = 1.0 / (2 * Math.PI) * Math.atan(x2 / x1);
                } else {
                    theta = 1.0 / (2 * Math.PI) * Math.atan(x2 / x1) + 0.5;
                }
                return 100.0 * ((x3 - 10.0 * theta) * (x3 - 10.0 * theta) + (Math.sqrt(x1 * x1 + x2 * x2) - 1) * (Math.sqrt(x1 * x1 + x2 * x2) - 1)) + x3 * x3;
            }

            @Override
            public DMatrixRMaj getG(DMatrixRMaj x) {
                return null;
            }

            @Override
            public DMatrixRMaj getH(DMatrixRMaj x) {
                return null;
            }
        };
        /* set up initial guess */
        DMatrixRMaj initialGuess = new DMatrixRMaj(3, 1);
        initialGuess.set(0, 0, -1.0 * initialMultiplier);
        initialGuess.set(1, 0, 0 * initialMultiplier);
        initialGuess.set(2, 0, 0 * initialMultiplier);
        Options options = new Options(3);
        /* options */
        options.setAlgorithm(algorithm);
        options.setAnalyticalGradient(false);
        options.setAnalyticalHessian(false);
        options.setSaveIterationDetails(true);
        /* start optimizer */
        UnconstrainedOptimizer solver = new UnconstrainedOptimizer(function, options);
        solver.solve(new DMatrixRMaj(initialGuess));
        System.out.println(solver.getResults());
        System.out.println("F: " + solver.getFx());
        System.out.println("x: " + solver.getX());
    }

    public static void woodFunction(int algorithm, boolean analyticalGradient, boolean analyticalHessian, double initialMultiplier) {
        ObjectiveFunctionUnconstrained function = new ObjectiveFunctionUnconstrained() {
            @Override
            public double getF(DMatrixRMaj x) {
                double x1 = x.get(0);
                double x2 = x.get(1);
                double x3 = x.get(2);
                double x4 = x.get(3);
                return 100 * (x1 * x1 - x2) * (x1 * x1 - x2) + (1 - x1) * (1 - x1) + 90 * (x3 * x3 - x4) * (x3 * x3 - x4) + (1 - x3) * (1 - x3)
                        + 10.1 * ((1 - x2) * (1 - x2) + (1 - x4) * (1 - x4)) + 19.8 * (1 - x2) * (1 - x4);
            }

            @Override
            public DMatrixRMaj getG(DMatrixRMaj x) {
                double x1 = x.get(0);
                double x2 = x.get(1);
                double x3 = x.get(2);
                double x4 = x.get(3);
                DMatrixRMaj g = new DMatrixRMaj(4, 1);
                g.set(0, 0, 2 * x1 - 400.0 * x1 * (-x1 * x1 + x2) - 2);
                g.set(1, 0, -200.0 * x1 * x1 + (1101.0 * x2) / 5.0 + (99.0 * x4) / 5.0 - 40.0);
                g.set(2, 0, 2.0 * x3 - 360.0 * x3 * (-x3 * x3 + x4) - 2.0);
                g.set(3, 0, -180.0 * x3 * x3 + (99.0 * x2) / 5.0 + (1001.0 * x4) / 5.0 - 40.0);
                return g;
            }

            @Override
            public DMatrixRMaj getH(DMatrixRMaj x) {
                double x1 = x.get(0);
                double x2 = x.get(1);
                double x3 = x.get(2);
                double x4 = x.get(3);
                DMatrixRMaj h = new DMatrixRMaj(4, 4);
                h.set(0, 0, 1200.0 * x1 * x1 - 400.0 * x2 + 2);
                h.set(0, 1, -400 * x1);
                h.set(1, 0, -400 * x1);
                h.set(1, 1, 1101.0 / 5.0);
                h.set(1, 3, 99.0 / 5.0);
                h.set(2, 2, 1080.0 * x3 * x3 - 360.0 * x4 + 2.0);
                h.set(2, 3, -360.0 * x3);
                h.set(3, 1, 99.0 / 5.0);
                h.set(3, 2, -360.0 * x3);
                h.set(3, 3, 1001.0 / 5.0);
                return h;
            }
        };
        DMatrixRMaj initialGuess = new DMatrixRMaj(4, 1);
        initialGuess.set(0, 0, -3.0 * initialMultiplier);
        initialGuess.set(1, 0, -1.0 * initialMultiplier);
        initialGuess.set(2, 0, -3.0 * initialMultiplier);
        initialGuess.set(3, 0, -1.0 * initialMultiplier);
        Options options = new Options(4);
        options.setAnalyticalGradient(analyticalGradient);
        options.setAnalyticalHessian(analyticalHessian);
        options.setAlgorithm(algorithm);
        options.setSaveIterationDetails(true);
        UnconstrainedOptimizer solver = new UnconstrainedOptimizer(function, options);
        for (int i = 0; i < 1; i++) {
            solver.solve(new DMatrixRMaj(initialGuess));
        }
        solver.solve(new DMatrixRMaj(initialGuess));
        System.out.println(solver.getResults());
        System.out.println("F: " + solver.getFx());
        System.out.println("x: " + solver.getX());
    }

    public static void rosenbrockFunction(int algorithm, boolean analyticalGradient, boolean analyticalHessian, double initialMultiplier) {
        ObjectiveFunctionUnconstrained function = new ObjectiveFunctionUnconstrained() {
            @Override
            public double getF(DMatrixRMaj x) {
                double x1 = x.get(0);
                double x2 = x.get(1);
                double x3 = x.get(2);
                double x4 = x.get(3);
                return (100.0 * (x2 - x1 * x1) * (x2 - x1 * x1) + (x1 - 1.0) * (x1 - 1.0)) + (100.0 * (x3 - x2 * x2) * (x3 - x2 * x2) + (x2 - 1.0) * (x2 - 1.0)) + (100.0 * (x4 - x3 * x3) * (x4 - x3 * x3) + (x3 - 1.0) * (x3 - 1.0));
            }

            @Override
            public DMatrixRMaj getG(DMatrixRMaj x) {
                double x1 = x.get(0);
                double x2 = x.get(1);
                double x3 = x.get(2);
                double x4 = x.get(3);
                DMatrixRMaj g = new DMatrixRMaj(4, 1);
                g.set(0, 0, 2.0*x1 - 2.0*x1*(- 100.0*x1*x1 + 100*x2) - 200*x1*(- x1*x1 + x2) - 2);
                g.set(1, 0,- 200*x1*x1 + 202*x2 - 2*x2*(100*x3 - 100*x2*x2) - 200*x2*(x3 - x2*x2) - 2);
                g.set(2, 0, - 200*x2*x2 + 202*x3 - 2*x3*(100*x4 - 100*x3*x3) - 200*x3*(x4 - x3*x3) - 2);
                g.set(3, 0, - 200*x3*x3 + 200*x4);
                return g;
            }

            @Override
            public DMatrixRMaj getH(DMatrixRMaj x) {
                double x1 = x.get(0);
                double x2 = x.get(1);
                double x3 = x.get(2);
                double x4 = x.get(3);
                DMatrixRMaj h = new DMatrixRMaj(4, 4);
                h.set(0, 0, 1200*x1*x1 - 400*x2 + 2);
                h.set(0, 1, -400 * x1);
                h.set(1, 0, -400 * x1);
                h.set(1, 1, 1200.0*x2*x2-400*x3+202);
                h.set(1, 3, -400*x2);
                h.set(2, 1,  -400*x2);
                h.set(2, 2, 1200*x3*x3 - 400*x4 + 202);
                h.set(2, 3, -400*x3);
                h.set(3, 2, -400*x3);
                h.set(3, 3, 200);
                return h;
            }
        };
        /* set up initial guess */
        DMatrixRMaj initialGuess = new DMatrixRMaj(4, 1);
        initialGuess.set(0, 0, -2.0 * initialMultiplier);
        initialGuess.set(1, 0, 2.0 * initialMultiplier);
        initialGuess.set(2, 0, -2.0 * initialMultiplier);
        initialGuess.set(3, 0, 2.0 * initialMultiplier);
        Options options = new Options(4);
        /* options */
        options.setAllTolerances(1e-12);
        options.setMaxIterations(10000);
        options.setAlgorithm(algorithm);
        options.setAnalyticalGradient(analyticalGradient);
        options.setAnalyticalHessian(analyticalHessian);
        options.setSaveIterationDetails(true);
        /* start optimizer */
        UnconstrainedOptimizer solver = new UnconstrainedOptimizer(function, options);
        solver.solve(new DMatrixRMaj(initialGuess));
        System.out.println(solver.getResults());
        System.out.println("F: " + solver.getFx());
        System.out.println("x: " + solver.getX());
    }
    
     public static void powellSingularFunction(int algorithm, boolean analyticalGradient, boolean analyticalHessian, double initialMultiplier) {
        ObjectiveFunctionUnconstrained function = new ObjectiveFunctionUnconstrained() {
            @Override
            public double getF(DMatrixRMaj x) {
                double x1 = x.get(0);
                double x2 = x.get(1);
                double x3 = x.get(2);
                double x4 = x.get(3);
                return Math.pow((x1+10*x2),2)+5*Math.pow((x3-x4),2)+Math.pow((x2-2*x3),4)+10*Math.pow((x1-x4),4);
            }

            @Override
            public DMatrixRMaj getG(DMatrixRMaj x) {
                double x1 = x.get(0);
                double x2 = x.get(1);
                double x3 = x.get(2);
                double x4 = x.get(3);
                DMatrixRMaj g = new DMatrixRMaj(4, 1);
                g.set(0, 0, 2*x1 + 20*x2 + 40*Math.pow((x1 - x4),3));
                g.set(1, 0,20*x1 + 200*x2 + 4*Math.pow((x2 - 2*x3),3));
                g.set(2, 0, 10*x3 - 10*x4 - 8*Math.pow((x2 - 2*x3),3));
                g.set(3, 0,  10*x4 - 10*x3 - 40*Math.pow((x1 - x4),3));
                return g;
            }

            @Override
            public DMatrixRMaj getH(DMatrixRMaj x) {
                double x1 = x.get(0);
                double x2 = x.get(1);
                double x3 = x.get(2);
                double x4 = x.get(3);
                DMatrixRMaj h = new DMatrixRMaj(4, 4);
                h.set(0, 0, 120*(x1 - x4)*(x1 - x4) + 2);
                h.set(0, 1, 20);
                h.set(0, 3, -120*(x1 - x4)*(x1 - x4));
                h.set(1, 0 ,20);
                h.set(1, 1, 12*(x2 - 2*x3)*(x2 - 2*x3) + 200);
                h.set(1, 2,   -24*(x2 - 2*x3)*(x2 - 2*x3));
                h.set(2, 1, -24*(x2 - 2*x3)*(x2 - 2*x3));
                h.set(2, 2, 48*(x2 - 2*x3)*(x2 - 2*x3) + 10);
                h.set(2, 3, -10);
                h.set(3, 0,  -120*(x1 - x4)*(x1 - x4));
                h.set(3, 2,  -10);
                h.set(3, 3,  120*(x1 - x4)*(x1 - x4) + 10);
                return h;
            }
        };
        /* set up initial guess */
        DMatrixRMaj initialGuess = new DMatrixRMaj(4, 1);
        initialGuess.set(0, 0, 3.0);
        initialGuess.set(1, 0,-1.0);
        initialGuess.set(2, 0, 0.0);
        initialGuess.set(3, 0, 1.0);
        Options options = new Options(4);
        /* options */
        options.setAllTolerances(1e-12);
        options.setAlgorithm(algorithm);
        options.setAnalyticalGradient(analyticalGradient);
        options.setAnalyticalHessian(analyticalHessian);
        options.setSaveIterationDetails(true);
        /* start optimizer */
        UnconstrainedOptimizer solver = new UnconstrainedOptimizer(function, options);
        solver.solve(new DMatrixRMaj(initialGuess));
        System.out.println(solver.getResults());
        System.out.println("F: " + solver.getFx());
        System.out.println("x: " + solver.getX());
    }



}
