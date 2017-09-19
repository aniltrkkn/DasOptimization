/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package optimization;

import optimization.functionImplementation.Options;
import test.NonlinearTest;
import test.UnconstrainedTest;

/**
 *
 * @author O. Anil Turkkan <turkkan.1@osu.edu> - anilturkkan.com
 */
public class Optimization {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
      // UnconstrainedTest.bealeFunction(Options.DOGLEG_TRUST_REGION, true, true, 1.0);
      // UnconstrainedTest.helicalValleyFunction(Options.DOGLEG_TRUST_REGION, true, true, 1.0);
      // UnconstrainedTest.woodFunction(Options.DOGLEG_TRUST_REGION, false, false, 1.0);
      UnconstrainedTest.powellSingularFunction(Options.DOGLEG_TRUST_REGION, true, true, 1.0);
    }

}
