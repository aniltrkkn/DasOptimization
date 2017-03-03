/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package optimization.functionImplementation;

import cern.colt.matrix.DoubleMatrix1D;
import cern.colt.matrix.DoubleMatrix2D;

/**
 *
 * @author anill
 */
public interface ObjectiveFunction {
    /*objective function*/
    public DoubleMatrix1D getF(DoubleMatrix1D x);
    public DoubleMatrix1D getD(DoubleMatrix1D x);
    public DoubleMatrix2D getH(DoubleMatrix1D x);
}
