/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package descentAlgorithms;

import optimization.functionImplementation.Options;
import org.ejml.data.DMatrixRMaj;
import solvers.Solver;

/**
 *
 * @author turkkan.1
 */
public interface StepAlgorithm {
    public boolean isMaxStepTaken();
    public boolean isSolverFailed();
    public DMatrixRMaj solve(DMatrixRMaj g, DMatrixRMaj sn, DMatrixRMaj x, DMatrixRMaj lowerTriangleR, Options solverOptions, Solver solver);
}
