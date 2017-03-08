/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package optimization;
import optimization.functionImplementation.Options;
import test.NonlinearTest;

/**
 *
 * @author O. Anil Turkkan <turkkan.1@osu.edu> - anilturkkan.com
 */
public class Optimization {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        //
       // NonlinearTest.testExtendedRosenbrockFunction(500,Options.DOGLEG_TRUST_REGION,100);
       // NonlinearTest.testPowellSingularFunction(400, Options.LINE_SEARCH, 1);
      // NonlinearTest.testPowellSingularFunction(10, Options.DOGLEG_TRUST_REGION, 1);
      NonlinearTest.testHelicalValleyFunction(3, Options.LINE_SEARCH, 10);
    }

}
