/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package optimization.functionImplementation;

import org.ejml.data.DenseMatrix64F;



/**
 *
 * @author anill
 */
public interface ObjectiveFunction {
    /*objective function*/
    public DenseMatrix64F getF(DenseMatrix64F x);
    public DenseMatrix64F getD(DenseMatrix64F x);
    public DenseMatrix64F getH(DenseMatrix64F x);
}
