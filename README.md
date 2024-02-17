# DasOptimization

## Differences
The difference between my version and the original is that I moved the library to Gradle and <a href="https://jitpack.io/">jitpack.io</a> so that it would be easier for people to use this library in their projects. Also here is the user guide for DasOptimization.

## Download

### Gradle:
- Add the JitPack repository to your `build.gradle` file:
```gradle
repositories {
    ...
    maven { url 'https://jitpack.io' }
}
```
- Add the dependency:
```gradle
dependencies {
    ...
    implementation 'com.github.Tamada4a:DasOptimization:1.0'
}
```

### Maven:
- Add the JitPack repository to your `pom.xml` file:
```xml
<repositories>
    ...
    <repository>
        <id>jitpack.io</id>
        <url>https://jitpack.io</url>
    </repository>
</repositories>
```
- Add the dependency:
```xml
<dependencies>
    ...
    <dependency>
        <groupId>com.github.Tamada4a</groupId>
        <artifactId>DasOptimization</artifactId>
         <version>1.0</version>
    </dependency>
</dependencies>
```

## User Guide for DasOptimization
### Solving System of Nonlinear Equations
You need to use <em>NonlinearEquationSolver</em> class to solve a system of nonlinear equations. Before creating a <em>NonlinearEquationSolver object</em>, you need to define your set of equations and the Jacobian (if available). You also need to create an <em>Options</em> object to specify user defined options for the solver.

#### Options Object
Create an <em>Options</em> object as:
```java
Options options = new Options(n);
```
where <em>n</em> is the number of variables.<br/>

This will be enough for most applications. However, you can fully customize solver preferences and here are some of the important parameters:
```java
options.setAnalyticalJacobian(true); //specify if you will supply the analytical Jacobian (default:false)
options.setAlgorithm(Options.TRUST_REGION); //set the algorithm; Options.TRUST_REGION or Options.LINE_SEARCH (default: Options.TRUST_REGION)
options.setSaveIterationDetails(true); //save iteration details to a Results object (default:false)
options.setAllTolerances(1e-12); //set convergence tolerances (default:1e-8)
options.setMaxIterations(1000); //set maximum number of iterations (default:100)
```

#### ObjectiveFunctionNonLinear Object
In order to define a set of nonlinear equations that you want to solve, you need to implement the <em>ObjectiveFunctionNonLinear</em> interface.  You also have to define the Jacobian matrix if you had set the analytical Jacobian flag at the <em>Options</em> object.
```java
    ObjectiveFunctionNonLinear f = new ObjectiveFunctionNonLinear() {
        @Override
        public DMatrixRMaj getF(DMatrixRMaj x) {
            DMatrixRMaj f = new DMatrixRMaj(2, 1);
            f.set(0, 0, 10 * (x.get(1, 0) - x.get(0, 0) * x.get(0, 0)));
            f.set(1, 0, 1 - x.get(0, 0));
            return f;
        }

        @Override
        public DMatrixRMaj getJ(DMatrixRMaj x) {
            DMatrixRMaj J = new DMatrixRMaj(numberOfVariables, numberOfVariables);
            J.set(0,0, -20 * x.get(0);
            J.set(0,1, 10);
            J.set(1, 0, -1);
            return J;
        }
    }
```
If you do not want to define the Jacobian and want to use finite difference Jacobian, just return <em>null</em> in <em>getJ</em> method.

#### Solving System of Nonlinear Equations
Last step before solving the set of nonlinear equations is defining an initial guess vector:
```java
DMatrixRMaj initialGuess = new DMatrixRMaj(2, 1);
initialGuess.set(0, 0, -1.2);
initialGuess.set(1, 0, 1.0);
```
Now, we are ready to create an <em>NonlinearEquationSolver</em> object using previously defined option and equation definitions:
```java
NonlinearEquationSolver nonlinearSolver = new NonlinearEquationSolver(f, options);
```
Finally, solve the system of nonlinear equations at the initial guess:
```java
nonlinearSolver.solve(new DMatrixRMaj(initialGuess));
```

### Nonlinear Unconstrained Optimization
You need to use <em>UnconstrainedOptimizer</em> class for solving a nonlinear unconstrained optimization problem. Before creating a <em>UnconstrainedOptimizer object</em>, you need to define the objective function, gradient vector (if available) and the Hessian matrix (if available). You also need to create an <em>Options</em> object to specify user defined options for the nonlinear unconstrained optimization solver.

#### Options Object
Create an <em>Options</em> object as:
```java
Options options = new Options(n);
```
where <em>n</em> is the number of variables.<br/>

Options for the nonlinear unconstrained optimization include specifying if analytical gradient and/or Hessian are implemented:
```java
options.setAnalyticalGradient(true); //specify if you will supply the analytical gradient (default:false)
options.setAnalyticalHessian(true); //specify if you will supply the analytical Hessian (default:false)
options.setAlgorithm(Options.TRUST_REGION); //set the algorithm; Options.TRUST_REGION or Options.LINE_SEARCH (default: Options.TRUST_REGION)
options.setSaveIterationDetails(true); //save iteration details to a Results object (default:false)
options.setAllTolerances(1e-12); //set convergence tolerances (default:1e-8)
options.setMaxIterations(1000); //set maximum number of iterations (default:100)
```

#### ObjectiveFunctionUnconstrained Object
In order to define a set of nonlinear equations that you want to solve, you need to implement the <em>ObjectiveFunctionNonLinear</em> interface.  You also have to define the gradient vector and/or the Hessian matrix if you set the appropriate flags in the <em>Options</em> object.
```java
ObjectiveFunctionUnconstrained function = new ObjectiveFunctionUnconstrained() {
        @Override
        public double getF(DMatrixRMaj x) {
            double x1 = x.get(0);
            double x2 = x.get(1);
            return (1.5 - x1 + x1 * x2) * (1.5 - x1 + x1 * x2) + (2.25 - x1 + x1 * x2 * x2) * (2.25 - x1 + x1 * x2 * x2) + (2.625 - x1 + x1 * x2 * x2 * x2) * (2.625 - x1 + x1 * x2 * x2 * x2);
        }

        @Override
        public DMatrixRMaj getG(DMatrixRMaj x) {
            double x1 = x.get(0);
            double x2 = x.get(1);
            DMatrixRMaj g = new DMatrixRMaj(2, 1);
            g.set(0, 0, 2.0 * (x2 * x2 - 1.0) * (x1 * x2 * x2 - x1 + 9.0 / 4.0) + 2.0 * (x2 * x2 * x2 - 1) * (x1 * x2 * x2 * x2 - x1 + 21.0 / 8.0) + 2 * (x2 - 1) * (x1 * x2 - x1 + 3.0 / 2.0));
            g.set(1, 0, 2.0 * x1 * (x1 * x2 - x1 + 3.0 / 2.0) + 4.0 * x1 * x2 * (x1 * x2 * x2 - x1 + 9.0 / 4.0) + 6 * x1 * x2 * x2 * (x1 * x2 * x2 * x2 - x1 + 21.0 / 8.0));
            return g;
        }

        @Override
        public DMatrixRMaj getH(DMatrixRMaj x) {
            double x1 = x.get(0);
            double x2 = x.get(1);
            DMatrixRMaj h = new DMatrixRMaj(2, 2);
            h.set(0, 0, 2 * (x2 - 1.0) * (x2 - 1.0) + 2.0 * (x2 * x2 - 1.0) * (x2 * x2 - 1.0) + 2.0 * (x2 * x2 * x2 - 1.0) * (x2 * x2 * x2 - 1.0));
            h.set(0, 1, 9.0 * x2 - 4.0 * x1 - 4.0 * x1 * x2 - 12.0 * x1 * x2 * x2 + 8.0 * x1 * x2 * x2 * x2 + 12 * x1 * x2 * x2 * x2 * x2 * x2 + (63.0 * x2 * x2) / 4.0 + 3.0);
            h.set(1, 0, 9.0 * x2 - 4.0 * x1 - 4.0 * x1 * x2 - 12.0 * x1 * x2 * x2 + 8.0 * x1 * x2 * x2 * x2 + 12.0 * x1 * x2 * x2 * x2 * x2 * x2 + (63.0 * x2 * x2) / 4.0 + 3.0);
            h.set(1, 1, (x1 * (63.0 * x2 - 4.0 * x1 - 24.0 * x1 * x2 + 24.0 * x1 * x2 * x2 + 60.0 * x1 * x2 * x2 * x2 * x2 + 18.0)) / 2.0);
            return h;
        }
    };
```
Return <em>null</em> in the gradient and/or the Hessian function if you want to use finite element difference approximations.

#### Unconstrained Nonlinear Optimization
Last step before solving the set of nonlinear equations is defining an initial guess vector:
```java
DMatrixRMaj initialGuess = new DMatrixRMaj(2, 1);
initialGuess.set(0, 0, -4.5);
initialGuess.set(1, 0, -4.5);
```
Now, we are ready to create an <em>UnconstrainedOptimizer</em> object using previously defined option and equation definitions:
```java
UnconstrainedOptimizer unconstrainedSolver= new UnconstrainedOptimizer(f, options);
```
Finally, solve the nonlinear unconstrained optimization problem:
```java
unconstrainedSolver.solve(new DMatrixRMaj(initialGuess));
```

### Post-Processing
After calling <em>unconstrainedSolver.solve(initialGuess)</em> or <em>nonlinearSolver.solve(initialGuess)</em>, you can access the solution details by:
```java
solver.getX(); //return final x values
solver.getFx(); //return final function values
solver.getJacobian(); //return final Jacobian matrix (only in nonlinear equations solver)
solver.getGx(); // return final gradient vector (only in nonlinear unconstrained optimization)
solver.getHx(); // return final Hessian matrix (only in nonlinear unconstrained optimization)
solver.getTerminationString(); // return the convergence (or failure) details
```

#### Results Object
You can set to save the iteration details (minimal effect on performance) if you set:
```java
options.setSaveIterationDetails(true);
```
Iteration details can be obtained by:
```java
Results results = solver.getResults(); //get results from the solver
solver.getX(); //get list of x values
solver.getFunctionNorm(); //get list of gradient vectors
solver.getFunctionEvaluations(); //get number of function evaluations
```

### System of Nonlinear Equations Solver Tests
<em>NonlinearTest</em> class under the <em>test</em> package contains a number of complicated test problems. The general structure of the test problems are:
```java
NonlinearTest.test(int numberOfVariables, int solver, boolean analyticalJacobian);
```
where <em>numberOfVariables</em> are number of equations, <em>solver</em> is either Options.TRUST_REGION or Options.LINE_SEARCH and <em>analyticalJacobian</em> should be set to true if analytical <em>Jacobian</em> is desired for the problem.<br/>

The list of problems are:
```java
NonlinearTest.extendedRosenbrockFunction(int numberOfVariables, int solver, boolean analyticalJacobian) //numberOfVariables must be multiple of two
NonlinearTest.powellSingularFunction(int numberOfVariables, int solver, boolean analyticalJacobian) //numberOfVariables must be multiple of four, analytical Jacobian not available
NonlinearTest.trigonometricFunction(int numberOfVariables, int solver, boolean analyticalJacobian) //analytical Jacobian not available
NonlinearTest.helicalValleyFunction(int numberOfVariables, int solver, boolean analyticalJacobian) //numberOfVariables does not have any affect, analytical Jacobian not available
```

### Nonlinear Unconstrained Optimization Tests
<em>UnconstrainedTest</em> class under the <em>test</em> package contains a number of complicated test problems. The general structure of the test problems are:
```java
NonlinearTest.test(int algorithm, boolean analyticalGradient, boolean analyticalHessian);
```
where <em>algorithm</em> is either Options.TRUST_REGION or Options.LINE_SEARCH, <em>analyticalGradient</em> and <em>analyticalHessian</em> should be set to true if analytical gradient and/or Hessian is desired for the problem.<br/>

The list of problems are:
```java
UnconstrainedTest.bealeFunction(int algorithm, boolean analyticalGradient, boolean analyticalHessian)
UnconstrainedTest.helicalValleyFunction(int algorithm, boolean analyticalGradient, boolean analyticalHessian) //analytical gradient and Hessian are not available
UnconstrainedTest.woodFunction(int algorithm, boolean analyticalGradient, boolean analyticalHessian)
UnconstrainedTest.rosenbrockFunction(int algorithm, boolean analyticalGradient, boolean analyticalHessian)
UnconstrainedTest.powellSingularFunction(int algorithm, boolean analyticalGradient, boolean analyticalHessian)
```

## Code Examples
Code example of solving system of nonlinear equations:
1. [Main class](https://github.com/Tamada4a/DasOptimization/blob/master/src/Main.java) of this repository.
2. [My study project on Kotlin](https://github.com/Tamada4a/SimovinARMA/blob/main/Kotlin/src/main/kotlin/Main.kt).

## Authors
- Original <a href="https://github.com/aniltrkkn/DasOptimization">repository</a>.
- <a href="https://github.com/aniltrkkn/DasOptimization">Author</a> of original repository.

## Original README
DasOptimization (V1.0) is a lightweight, robust and scalable library capable of solving a system of nonlinear equations and nonlinear unconstrained optimization. Some of the capabilities are:

- DasOptimization implemented two descent algorithms; Line Search and Trust-Region-Dogleg algorithm.
- Analytical gradient, the Jacobian or the Hessian can be supplied but if they are not available, they will be calculated with finite differences. 
- DasOptimization can recover from singular Hessian and Jacobian matrices.
- DasOptimization is coded only in Java and therefore, it can be used in Android projects.

For more information, please visit the [user guide](http://compliantanalysis.com/dasOptimization).

DasOptimization is created for and successfully employed in [DAS Mechanism Analysis](http://compliantanalysis.com) which is a fast kinematic analysis software for mechanisms.

EJML ([Efficient Java Matrix Library](http://ejml.org/wiki/index.php?title=Main_Page)) is used for matrix calculations. You need to include or build appropriate .jar files in your project library before running any DasOptimization code. 

DasOptimization (V2.0) will include a nonlinear constrained optimizer capable of handling nonlinear constraints. 

DasOptimization is developed by [Design, Innovation and Simulation Lab.](https://disl.osu.edu/) of the Ohio State University.
