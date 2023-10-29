import optimization.functionImplementation.ObjectiveFunctionNonLinear;
import optimization.functionImplementation.Options;
import org.ejml.data.DMatrixRMaj;
import solvers.NonlinearEquationSolver;

public class Main {
    public static void main(String[] args) {
        int var = 4;
        Options options = new Options(var);
        options.setAnalyticalJacobian(false); //specify if you will supply the analytical Jacobian (default:false)
        options.setAlgorithm(Options.TRUST_REGION); //set the algorithm; Options.TRUST_REGION or Options.LINE_SEARCH (default: Options.TRUST_REGION)
        options.setSaveIterationDetails(true); //save iteration details to a Results object (default:false)
        options.setAllTolerances(1e-12); //set convergence tolerances (default:1e-8)
        options.setMaxIterations(1000); //set maximum number of iterations (default:100)

        ObjectiveFunctionNonLinear f = new ObjectiveFunctionNonLinear() {
            @Override
            public DMatrixRMaj getF(DMatrixRMaj x) {
                DMatrixRMaj f = new DMatrixRMaj(var, 1);
                f.set(0, 0, x.get(0, 0) * x.get(0, 0) + x.get(1, 0) * x.get(1, 0) + x.get(2, 0) * x.get(2, 0) + x.get(3, 0) * x.get(3, 0) - 339.78679);
                f.set(1, 0, x.get(1, 0) * x.get(0, 0) + x.get(1, 0) * x.get(2, 0) + x.get(3, 0) * x.get(2, 0) - 97.29910);
                f.set(2, 0, x.get(0, 0) * x.get(2, 0) + x.get(1, 0) * x.get(3, 0) + 85.07386);
                f.set(3, 0, x.get(0, 0) * x.get(3, 0) + 28.87565);
                return f;
            }

            @Override
            public DMatrixRMaj getJ(DMatrixRMaj x) {
                return null;
            }
        };

        DMatrixRMaj initialGuess = new DMatrixRMaj(var, 1);
        initialGuess.set(0, 0, 18.433307);
        initialGuess.set(1, 0, 0.0);
        initialGuess.set(2, 0, 0.0);
        initialGuess.set(3, 0, 0.0);

        NonlinearEquationSolver nonlinearSolver = new NonlinearEquationSolver(f, options);
        //solve and print output
        nonlinearSolver.solve(new DMatrixRMaj(initialGuess));
        System.out.println("Results: " + nonlinearSolver.getResults());
        System.out.println("x: " + nonlinearSolver.getX());
    }
}