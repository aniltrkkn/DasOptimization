/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package optimization.functionImplementation;

import org.ejml.data.DMatrixRMaj;



/**
 *
 * @author anill
 */
public interface ObjectiveFunctionUnconstrained {
    /*objective function*/
    public double getF(DMatrixRMaj x);
    public DMatrixRMaj getG(DMatrixRMaj x);
    public DMatrixRMaj getH(DMatrixRMaj x);
}
