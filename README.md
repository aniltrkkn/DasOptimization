# DasOptimization

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
