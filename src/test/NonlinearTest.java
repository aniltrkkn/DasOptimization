/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package test;

import optimization.functionImplementation.Options;
import org.ejml.data.DMatrixRMaj;
import solvers.NonlinearEquationSolver;
import optimization.functionImplementation.ObjectiveFunctionNonLinear;
import optimization.functionImplementation.ObjectiveFunctionUnconstrained;
import solvers.UnconstrainedOptimizer;

/**
 *
 * @author turkkan.1
 */
public class NonlinearTest {

    public static void testExtendedRosenbrockFunction(int n, int solver, double initialGuessMultiplier) {
        System.out.println("**********************************************");
        System.out.println("Number of Variables=" + Integer.toString(n));
        ObjectiveFunctionNonLinear f = new ObjectiveFunctionNonLinear() {
            @Override
            public DMatrixRMaj getF(DMatrixRMaj x) {
                DMatrixRMaj f = new DMatrixRMaj(n, 1);
                for (int i = 0; i < n / 2; i++) {
                    f.set(2 * i, 0, 10 * (x.get(2 * i + 1, 0) - x.get(2 * i, 0) * x.get(2 * i, 0)));
                    f.set(2 * i + 1, 0, 1 - x.get(2 * i, 0));
                }
                return f;
            }

            @Override
            public DMatrixRMaj getJ(DMatrixRMaj x) {
                DMatrixRMaj H = new DMatrixRMaj(n, n);
                for (int i = 0; i < n / 2; i++) {
                    H.set(2 * i, 2 * i, -20 * x.get(2 * i));
                    H.set(2 * i, 2 * i + 1, 10);
                    H.set(2 * i + 1, 2 * i, -1);
                }

                return H;
            }

        };
        DMatrixRMaj initialGuess = new DMatrixRMaj(n, 1);
        for (int i = 0; i < n / 2; i++) {
            initialGuess.set(2 * i, -1.2 * initialGuessMultiplier);
            initialGuess.set(2 * i + 1, 1.0 * initialGuessMultiplier);
        }
        Options options = new Options(n);
        options.setAnalyticalHessian(false);
        options.setAlgorithm(solver);
        NonlinearEquationSolver nonlinearSolver = new NonlinearEquationSolver(f, options);
        long startTime = System.nanoTime();
        for (int i = 0; i < 1; i++) {
            nonlinearSolver.solve(new DMatrixRMaj(initialGuess));
        }
        long endTime = System.nanoTime();
        System.out.println(nonlinearSolver.getX());
        //System.out.println(solver.getFx());
        //System.out.println(solver.getJacobian());
        System.out.println((endTime - startTime) / 1e9);
    }

    public static void testPowellSingularFunction(int n, int solver, double initialGuessMultiplier) {
        System.out.println("**********************************************");
        System.out.println("Number of Variables=" + Integer.toString(n));
        ObjectiveFunctionNonLinear f = new ObjectiveFunctionNonLinear() {
            @Override
            public DMatrixRMaj getF(DMatrixRMaj x) {
                DMatrixRMaj f = new DMatrixRMaj(n, 1);
                for (int i = 0; i < n / 4; i++) {
                    f.set(4 * i, 0, x.get(4 * i, 0) + 10 * x.get(4 * i + 1, 0));
                    f.set(4 * i + 1, 0, Math.sqrt(5) * (x.get(4 * i + 2, 0) - x.get(4 * i + 3, 0)));
                    f.set(4 * i + 2, 0, Math.pow(x.get(4 * i + 1, 0) - 2 * x.get(4 * i + 2, 0), 2));
                    f.set(4 * i + 3, 0, Math.sqrt(10) * Math.pow(x.get(4 * i, 0) - x.get(4 * i + 3, 0), 2));
                }
                return f;
            }

            @Override
            public DMatrixRMaj getJ(DMatrixRMaj x) {
                DMatrixRMaj H = new DMatrixRMaj(n, n);
                for (int i = 0; i < n / 2; i++) {
                    H.set(2 * i, 2 * i, -20 * x.get(2 * i));
                    H.set(2 * i, 2 * i + 1, 10);
                    H.set(2 * i + 1, 2 * i, -1);
                }

                return H;
            }
        };
        DMatrixRMaj initialGuess = new DMatrixRMaj(n, 1);
        for (int i = 0; i < n / 4; i++) {
            initialGuess.set(4 * i, 3.0 * initialGuessMultiplier);
            initialGuess.set(4 * i + 1, -1.0 * initialGuessMultiplier);
            initialGuess.set(4 * i + 2, 0.0 * initialGuessMultiplier);
            initialGuess.set(4 * i + 3, 1.0 * initialGuessMultiplier);
        }
        Options options = new Options(n);
        options.setAnalyticalHessian(true);
        System.out.println("Algoritm: Line Search");
        options.setAlgorithm(solver);
        NonlinearEquationSolver nonlinearSolver = new NonlinearEquationSolver(f, options);
        long startTime = System.nanoTime();
        for (int i = 0; i < 1; i++) {
            nonlinearSolver.solve(new DMatrixRMaj(initialGuess));
        }
        long endTime = System.nanoTime();
        System.out.println(nonlinearSolver.getX());
        //System.out.println(solver.getFx());
        //System.out.println(solver.getJacobian());
        System.out.println((endTime - startTime) / 1e9);
    }

    public static void testTrigonometricFunction(int n, int solver, double initialGuessMultiplier) {
        System.out.println("**********************************************");
        System.out.println("Number of Variables=" + Integer.toString(n));
        ObjectiveFunctionNonLinear f = new ObjectiveFunctionNonLinear() {
            @Override
            public DMatrixRMaj getF(DMatrixRMaj x) {
                DMatrixRMaj f = new DMatrixRMaj(n, 1);
                for (int i = 0; i < n; i++) {
                    double temp = 0.0;
                    for (int j = 0; j < n; j++) {
                        temp += Math.cos(x.get(j, 0)) + ((double) i) * (1.0 - Math.cos(x.get(i, 0))) - Math.sin(x.get(i, 0));
                    }
                    f.set(i, 0, n - temp);
                }
                return f;
            }

            @Override
            public DMatrixRMaj getJ(DMatrixRMaj x) {
                DMatrixRMaj H = new DMatrixRMaj(n, n);
                return H;
            }
        };
        DMatrixRMaj initialGuess = new DMatrixRMaj(n, 1);
        for (int i = 0; i < n; i++) {
            initialGuess.set(i, 1 / n * initialGuessMultiplier);
        }
        Options options = new Options(n);
        options.setAnalyticalHessian(true);
        System.out.println("Algoritm: Line Search");
        options.setAlgorithm(solver);
        NonlinearEquationSolver nonlinearSolver = new NonlinearEquationSolver(f, options);
        long startTime = System.nanoTime();
        for (int i = 0; i < 1; i++) {
            nonlinearSolver.solve(new DMatrixRMaj(initialGuess));
        }
        long endTime = System.nanoTime();
        System.out.println(nonlinearSolver.getX());
        //System.out.println(solver.getFx());
        //System.out.println(solver.getJacobian());
        System.out.println((endTime - startTime) / 1e9);
    }

    public static void testHelicalValleyFunction(int n, int solver, double initialGuessMultiplier) {
        System.out.println("**********************************************");
        System.out.println("Number of Variables=" + Integer.toString(n));
        ObjectiveFunctionNonLinear f = new ObjectiveFunctionNonLinear() {
            @Override
            public DMatrixRMaj getF(DMatrixRMaj x) {
                DMatrixRMaj f = new DMatrixRMaj(n, 1);
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
                DMatrixRMaj H = new DMatrixRMaj(n, n);
                return H;
            }
        };
        DMatrixRMaj initialGuess = new DMatrixRMaj(n, 1);
        initialGuess.set(0, 0, -1.0 * initialGuessMultiplier);
        initialGuess.set(1, 0, 0.0 * initialGuessMultiplier);
        initialGuess.set(2, 0, 0.0 * initialGuessMultiplier);
        Options options = new Options(n);
        options.setAnalyticalHessian(false);
        System.out.println("Algoritm: Line Search");
        options.setAlgorithm(solver);
        NonlinearEquationSolver nonlinearSolver = new NonlinearEquationSolver(f, options);
        long startTime = System.nanoTime();
        for (int i = 0; i < 1; i++) {
            nonlinearSolver.solve(new DMatrixRMaj(initialGuess));
        }
        long endTime = System.nanoTime();
        System.out.println(nonlinearSolver.getX());
        //System.out.println(solver.getFx());
        //System.out.println(solver.getJacobian());
        System.out.println((endTime - startTime) / 1e9);
    }

    
}
