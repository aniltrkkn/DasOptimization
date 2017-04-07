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
        NonlinearTest.testExtendedRosenbrockFunction(500,Options.LINE_SEARCH,1);
        //NonlinearTest.testPowellSingularFunction(8, Options.testTrigonometricFunction, 1);
     // NonlinearTest.testTrigonometricFunction(2, Options.LINE_SEARCH, 1);
     //NonlinearTest.testHelicalValleyFunction(3, Options.LINE_SEARCH, 10);
    }

}
