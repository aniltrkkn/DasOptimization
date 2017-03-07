/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package solvers;

import org.ejml.data.DenseMatrix64F;

/**
 *
 * @author O. Anil Turkkan <turkkan.1@osu.edu> - anilturkkan.com
 */
public interface Solver {
     public double functionNorm(DenseMatrix64F x);
}
