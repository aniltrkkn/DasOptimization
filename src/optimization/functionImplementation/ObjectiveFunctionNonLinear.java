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
public interface ObjectiveFunctionNonLinear {
    /*objective function*/
    public DMatrixRMaj getF(DMatrixRMaj x);
    public DMatrixRMaj getJ(DMatrixRMaj x);
}
