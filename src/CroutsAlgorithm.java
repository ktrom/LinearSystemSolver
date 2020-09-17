import java.util.ArrayList;
import java.util.Scanner;

public class CroutsAlgorithm {
    public static void main(String[] args) {
        Decomposer p = new Decomposer();
        p.printLUDecomp();
        double[][] L = p.getL();
        double[][] U = p.getU();
        LinearEquationSolver l = new LinearEquationSolver(L, U);
        System.out.println();
        System.out.println(l.toString());
    }
}

class Decomposer {
    double[][] A;
    double[][] L;
    double[][] U;
    int n;

    public Decomposer() {
        getValues();
        LUdecomp();
    }

    public Decomposer(double[][] A) {
        this.A = A;
        n = A.length;
        L = new double[n][n];
        U = new double[n][n];
    }

    public void getValues(){
        Scanner scan = new Scanner(System.in);
        System.out.println("What is the dimension of your matrix?");
        n = scan.nextInt();
        A = new double[n][n];
        for(int i = 0; i < n; i++){
            for(int j = 0; j < n; j++){
                System.out.println("What is the value at location (" + (i+1) + ", " + (j+1) + ")?");
                A[i][j] = scan.nextDouble();
            }
        }

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
    }

    public static String formattedDouble(double d) {
        return String.format("%.4f", d);
    }

    public void printLUDecomp() {
        System.out.println(toString());
    }

    @Override
    public String toString(){
        String str = "";
        for (int j = 0; j < n; j++) {
            str += ("| ");
            for (int i = 0; i < n; i++) {
                str+=(formattedDouble(L[j][i]) + " ");
            }
           str+=("| ");
           str+=("| ");
            for (int i = 0; i < n; i++) {
               str+=(formattedDouble(U[j][i]) + " ");
            }
           str+=("|\n");
        }
        return str;
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
    private double[] b;
    private double[][] L;
    private double[][] U;
    private double[] y;
    private double[] x;
    private int n;

    public LinearEquationSolver(double[][] L, double[][] U){
        this.L = L;
        this.U = U;
        getValues();
    }

    private void getValues(){
        Scanner scan = new Scanner(System.in);
        System.out.println("What is your solution column vector?");
        n = L.length;
        b = new double[n];
        for(int i = 0; i < L.length; i++){
            System.out.println("Dimension " + (i + 1) + ":");
            b[i] = scan.nextDouble();
        }
        scan.close();
    }
    public LinearEquationSolver(double[] b, double[][] L, double[][] U) {
        this.b = b;
        this.L = L;
        this.U = U;
        n = L.length;
        y = null;
        x = null;
    }

    @Override
    public String toString(){
        if(x == null)
            getX();
        String str = "";
        System.out.println("   x " + " " + "    y ");
        for(int i = 0; i < n; i++){
            str+= (Decomposer.formattedDouble(x[i]) + " " + Decomposer.formattedDouble(y[i]) + "\n");
        }
        return str;
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




