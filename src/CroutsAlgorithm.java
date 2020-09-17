import java.util.ArrayList;

public class CroutsAlgorithm {
    public static void main(String[] args) {

        double[][] A = new double[][]{{11, 2, -5, 6, 48}, {1, 0, 17, 29, -21}, {-3, 4, 55, -61, 0}, {41, 97, -32, 47, 23}, {-6, 9, -4, -8, 50}};
        double[] b1 = new double[]{4, 0, -7, -2, -11};
        double[] b2 = new double[]{2, 77, -1003, -7, -10};
        Decomposer p = new Decomposer(A);
        p.LUdecomp();

        double[][] L = p.getL();
        double[][] U = p.getU();

        LinearEquationSolver l1 = new LinearEquationSolver(b1, L, U);
        double[] y1 = l1.getY();
        double[] x1 = l1.getX();
        LinearEquationSolver l2 = new LinearEquationSolver(b2, L, U);
        double[] y2 = l2.getY();
        double[] x2 = l2.getX();

        System.out.println();
        System.out.println("   x1" + " " + "    y1");
        for(int i = 0; i < y1.length; i++){
            System.out.println(Decomposer.formattedDouble(x1[i]) + " " + Decomposer.formattedDouble(y1[i]));
        }
        System.out.println();
        System.out.println("   x2" + " " + "    y2");
        for(int i = 0; i < y2.length; i++){
            System.out.println(Decomposer.formattedDouble(x2[i]) + " " + Decomposer.formattedDouble(y2[i]));
        }

    }
}

class Decomposer {
    double[][] A;
    double[][] L;
    double[][] U;
    int n;

    public Decomposer() {

    }

    public Decomposer(double[][] A) {
        this.A = A;
        n = A.length;
        L = new double[n][n];
        U = new double[n][n];
    }


    public double[][] getL() {
        return L;
    }

    public double[][] getU() {
        return U;
    }

    public void LUdecomp() {

        for (int j = 0; j < n; j++) {
            setUValues(j);
            int maximalRow = determineRow(j);
            switchRows(j, maximalRow);
            switchRowsOfL(j, maximalRow);
            setUValues(j);
            setLValues(j);
        }
        printLUDecomp();
    }

    public static String formattedDouble(double d) {
        return String.format("%.4f", d);
    }

    private void printLUDecomp() {
        for (int j = 0; j < n; j++) {
            System.out.print("| ");
            for (int i = 0; i < n; i++) {
                System.out.print(formattedDouble(L[j][i]) + " ");
            }
            System.out.print("| ");
            System.out.print("| ");
            for (int i = 0; i < n; i++) {
                System.out.print(formattedDouble(U[j][i]) + " ");
            }
            System.out.println("|");
        }
    }

    private void switchRows(int i, int j) {
        double temp;
        for (int k = 0; k < A.length; k++) {
            temp = A[i][k];
            A[i][k] = A[j][k];
            A[j][k] = temp;
        }
    }


    private void switchRowsOfL(int i, int j) {
        double temp;
        for (int k = 0; k < L.length; k++) {
            temp = L[i][k];
            L[i][k] = L[j][k];
            L[j][k] = temp;
        }
    }

    private int determineRow(int j) {
        Resulter results = new Resulter();

        for (int i = j; i < n; i++) {
            results.add(A[i][j] - epsilonSum1(i, j));
        }

        U[j][j] = results.maxValue();
        return j + results.maxValueIndex();
    }

    private void setUValues(int j) {
        for (int i = 0; i < j; i++) {
            U[i][j] = A[i][j] - epsilonSum1(i, j);
        }
    }

    private void setLValues(int j) {
        L[j][j] = 1;
        for (int i = j; i < n; i++) {
            L[i][j] = 1 / U[j][j] * (A[i][j] - epsilonSum2(i, j));
        }
    }

    private double epsilonSum1(int i, int j) {
        double sum = 0;
        for (int k = 0; k < i; k++) {
            sum += L[i][k] * U[k][j];
        }
        return sum;
    }

    private double epsilonSum2(int i, int j) {
        double sum = 0;
        for (int k = 0; k < j; k++) {
            sum += L[i][k] * U[k][j];
        }
        return sum;
    }

    private double epsilonSum3(int i, int j) {
        double sum = 0;
        for (int k = 0; k <= j; k++) {
            sum += A[i][k] * U[k][j];
        }
        return sum;
    }
}

class Resulter extends ArrayList<Double> {
    public Resulter() {
        super();
    }

    public double maxValue() {
        double value = 0;
        int index = 0;
        for (int i = 0; i < this.size(); i++) {
            double curValue = this.get(i);
            if (Math.abs(curValue) > Math.abs(value)) {
                value = curValue;
                index = i;
            }
        }
        return value;
    }

    public int maxValueIndex() {
        double value = 0;
        int index = 0;
        for (int i = 0; i < this.size(); i++) {
            double curValue = this.get(i);
            if (Math.abs(curValue) > Math.abs(value)) {
                value = curValue;
                index = i;
            }
        }
        return index;
    }
}

class LinearEquationSolver {
    private final double[] b;
    private final double[][] L;
    private final double[][] U;
    private double[] y;
    private double[] x;
    private final int n;

    public LinearEquationSolver(double[] b, double[][] L, double[][] U) {
        this.b = b;
        this.L = L;
        this.U = U;
        n = L.length;
        y = null;
        x = null;
    }

    // In the case Ly = b
    private double[] forwardSubstitution() {
        y = new double[n];
        for (int i = 0; i < n; i++) {
            y[i] = b[i] - epsilonSumForwardSub(i);
        }
        return y;
    }

    // In the case Ux = b
    private double[] backwardSubstitution() {
        x = new double[n];
        for (int i = n - 1; i > -1; i--) {
            x[i] = 1/U[i][i] * (y[i] - epsilonSumBackwardSub(i));
        }
        return x;
    }

    private double epsilonSumBackwardSub(int i) {
        double sum = 0;
        for (int j = i+1; j < n; j++) {
            sum += U[i][j] * x[j];
        }
        return sum;
    }

    private double epsilonSumForwardSub(int i) {
        double sum = 0;
        for (int j = 0; j < i; j++) {
            sum += L[i][j] * y[j] / L[i][i];
        }
        return sum;
    }

    public double[] getX() {
        if (x == null) {
            if (y == null) {
                forwardSubstitution();
            }
            backwardSubstitution();
        }
        return x;
    }

    public double[] getY() {
        if (y == null) {
            forwardSubstitution();
        }
        return y;
    }
}




