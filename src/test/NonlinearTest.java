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

import optimization.functionImplementation.Options;
import org.ejml.data.DMatrixRMaj;
import solvers.NonlinearEquationSolver;
import optimization.functionImplementation.ObjectiveFunctionNonLinear;


public class NonlinearTest {

    public static void extendedRosenbrockFunction(int numberOfVariables, int solver, boolean analyticalJacobian) {
        //input class
        ObjectiveFunctionNonLinear f = new ObjectiveFunctionNonLinear() {
            @Override
            public DMatrixRMaj getF(DMatrixRMaj x) {
                DMatrixRMaj f = new DMatrixRMaj(numberOfVariables, 1);
                for (int i = 0; i < numberOfVariables / 2; i++) {
                    f.set(2 * i, 0, 10 * (x.get(2 * i + 1, 0) - x.get(2 * i, 0) * x.get(2 * i, 0)));
                    f.set(2 * i + 1, 0, 1 - x.get(2 * i, 0));
                }
                return f;
            }

            @Override
            public DMatrixRMaj getJ(DMatrixRMaj x) {
                DMatrixRMaj H = new DMatrixRMaj(numberOfVariables, numberOfVariables);
                for (int i = 0; i < numberOfVariables / 2; i++) {
                    H.set(2 * i, 2 * i, -20 * x.get(2 * i));
                    H.set(2 * i, 2 * i + 1, 10);
                    H.set(2 * i + 1, 2 * i, -1);
                }

                return H;
            }

        };
        //initial guess
        DMatrixRMaj initialGuess = new DMatrixRMaj(numberOfVariables, 1);
        for (int i = 0; i < numberOfVariables/ 2; i++) {
            initialGuess.set(2 * i, -1.2 );
            initialGuess.set(2 * i + 1, 1.0);
        }
        //options
        Options options = new Options(numberOfVariables);
        options.setAnalyticalJacobian(analyticalJacobian);
        options.setAlgorithm(solver);
        options.setSaveIterationDetails(true);
        options.setAllTolerances(1e-12);
        NonlinearEquationSolver nonlinearSolver = new NonlinearEquationSolver(f, options);
        //solve and print output
        nonlinearSolver.solve(new DMatrixRMaj(initialGuess));
        System.out.println(nonlinearSolver.getResults());
        System.out.println("F: " + nonlinearSolver.getFx());
        System.out.println("x: " + nonlinearSolver.getX());

    }

    public static void powellSingularFunction(int numberOfVariables, int solver, boolean analyticalJacobian) {
        //input class
        ObjectiveFunctionNonLinear f = new ObjectiveFunctionNonLinear() {
            @Override
            public DMatrixRMaj getF(DMatrixRMaj x) {
                DMatrixRMaj f = new DMatrixRMaj(numberOfVariables, 1);
                for (int i = 0; i < numberOfVariables/ 4; i++) {
                    f.set(4 * i, 0, x.get(4 * i, 0) + 10 * x.get(4 * i + 1, 0));
                    f.set(4 * i + 1, 0, Math.sqrt(5) * (x.get(4 * i + 2, 0) - x.get(4 * i + 3, 0)));
                    f.set(4 * i + 2, 0, Math.pow(x.get(4 * i + 1, 0) - 2 * x.get(4 * i + 2, 0), 2));
                    f.set(4 * i + 3, 0, Math.sqrt(10) * Math.pow(x.get(4 * i, 0) - x.get(4 * i + 3, 0), 2));
                }
                return f;
            }

            @Override
            public DMatrixRMaj getJ(DMatrixRMaj x) {
                return null;
            }
        };
        //initial guess
        DMatrixRMaj initialGuess = new DMatrixRMaj(numberOfVariables, 1);
        for (int i = 0; i < numberOfVariables / 4; i++) {
            initialGuess.set(4 * i, 3.0 );
            initialGuess.set(4 * i + 1, -1.0 );
            initialGuess.set(4 * i + 2, 0.0 );
            initialGuess.set(4 * i + 3, 1.0 );
        }
        //options
        Options options = new Options(numberOfVariables);
        options.setAnalyticalJacobian(false);
        options.setAlgorithm(solver);
        options.setSaveIterationDetails(true);
        options.setAllTolerances(1e-12);
        NonlinearEquationSolver nonlinearSolver = new NonlinearEquationSolver(f, options);
        //solve and print output
        nonlinearSolver.solve(new DMatrixRMaj(initialGuess));
        System.out.println(nonlinearSolver.getResults());
        System.out.println("F: " + nonlinearSolver.getFx());
        System.out.println("x: " + nonlinearSolver.getX());
    }

    public static void trigonometricFunction(int numberOfVariables, int solver, boolean analyticalJacobian) {
        //input class
        ObjectiveFunctionNonLinear f = new ObjectiveFunctionNonLinear() {
            @Override
            public DMatrixRMaj getF(DMatrixRMaj x) {
                DMatrixRMaj f = new DMatrixRMaj(numberOfVariables, 1);
                for (int i = 0; i < numberOfVariables; i++) {
                    double temp = 0.0;
                    for (int j = 0; j < numberOfVariables; j++) {
                        temp += Math.cos(x.get(j, 0)) + ((double) i) * (1.0 - Math.cos(x.get(i, 0))) - Math.sin(x.get(i, 0));
                    }
                    f.set(i, 0, numberOfVariables - temp);
                }
                return f;
            }

            @Override
            public DMatrixRMaj getJ(DMatrixRMaj x) {
                return null;
            }
        };
        //initial guess
        DMatrixRMaj initialGuess = new DMatrixRMaj(numberOfVariables, 1);
        for (int i = 0; i < numberOfVariables; i++) {
            initialGuess.set(i, 1.0/ numberOfVariables);
        }
        //options
        Options options = new Options(numberOfVariables);
        options.setAnalyticalJacobian(false);
        options.setAlgorithm(solver);
        options.setSaveIterationDetails(true);
        options.setAllTolerances(1e-12);
        NonlinearEquationSolver nonlinearSolver = new NonlinearEquationSolver(f, options);
        //solve and print output
        nonlinearSolver.solve(new DMatrixRMaj(initialGuess));
        System.out.println(nonlinearSolver.getResults());
        System.out.println("F: " + nonlinearSolver.getFx());
        System.out.println("x: " + nonlinearSolver.getX());
    }

    public static void helicalValleyFunction(int numberOfVariables, int solver, boolean analyticalJacobian) {
        //input class
        ObjectiveFunctionNonLinear f = new ObjectiveFunctionNonLinear() {
            @Override
            public DMatrixRMaj getF(DMatrixRMaj x) {
                DMatrixRMaj f = new DMatrixRMaj(numberOfVariables, 1);
                double fValue;
                if (x.get(0, 0) > 0.0) {
                    fValue = (1 / (2 * Math.PI)) * Math.atan(x.get(1, 0) / x.get(0, 0));
                } else {
                    fValue = (1 / (2 * Math.PI)) * Math.atan(x.get(1, 0) / x.get(0, 0)) + 0.5;
                }
                f.set(0, 0, 10.0 * (x.get(2, 0) - 10 * fValue));
                f.set(1, 0, 10.0 * (Math.sqrt(x.get(0, 0) * x.get(0, 0) + x.get(1, 0) * x.get(1, 0)) - 1));
                f.set(2, 0, x.get(2, 0));
                return f;
            }

            @Override
            public DMatrixRMaj getJ(DMatrixRMaj x) {
                return null;
            }
        };
        //initial guess
        DMatrixRMaj initialGuess = new DMatrixRMaj(numberOfVariables, 1);
        initialGuess.set(0, 0, -1.0 );
        initialGuess.set(1, 0, 0.0 );
        initialGuess.set(2, 0, 0.0 );
        //options
        Options options = new Options(numberOfVariables);
        options.setAnalyticalJacobian(false);
        options.setAlgorithm(solver);
        options.setSaveIterationDetails(true);
        options.setAllTolerances(1e-12);
        NonlinearEquationSolver nonlinearSolver = new NonlinearEquationSolver(f, options);
        //solve and print output
        nonlinearSolver.solve(new DMatrixRMaj(initialGuess));
        System.out.println(nonlinearSolver.getResults());
        System.out.println("F: " + nonlinearSolver.getFx());
        System.out.println("x: " + nonlinearSolver.getX());
    }

    
}
