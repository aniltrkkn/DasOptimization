/*
 * Copyright (c) 2016-2017, O. Anil Turkkan. All Rights Reserved.
 *
 * This file is part of DASOptimization.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *   http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
package optimization;

import optimization.functionImplementation.Options;
import test.NonlinearTest;
import test.UnconstrainedTest;


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
