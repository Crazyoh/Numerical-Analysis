public class EulerModifiedEulerRungeKutta4  {
    public static void main(String[] args) {
        double x0 = 0; // ค่า x เริ่มต้น
        double y0 = 1;// ค่า y เริ่มต้น
        double h = 0.1; // ขนาดของขั้น
        double xTarget = 1; // ค่า x ปลายทาง

        System.out.println("x\t\tEuler\t\t\tModified Euler\tRunge-Kutta\t\tExact\t\teulerError\t\tmodifiedEulerError" +
                "\t\trungeKuttaError");
        for (double x = x0; x <= xTarget; x += h) {
            double eulerResult = euler(x0, y0, h, x);
            double modifiedEulerResult = modifiedEuler(x0, y0, h, x);
            double rungeKuttaResult = rungeKutta(x0, y0, h, x);
            double exactResult = 1/(x+1); // ผลเฉลยแม่นตรง

            double eulerError = Math.abs(exactResult - eulerResult);
            double modifiedEulerError = Math.abs(exactResult - modifiedEulerResult);
            double rungeKuttaError = Math.abs(exactResult - rungeKuttaResult);

            System.out.printf("%.1f\t\t%.6f\t\t%.6f\t\t%.6f\t\t%.6f\t%.6f\t\t%.6f\t\t\t\t%.6f\t\n", x, eulerResult, modifiedEulerResult,
                    rungeKuttaResult, exactResult, eulerError, modifiedEulerError, rungeKuttaError);
        }
    }
    public static double f(double x, double y) {
        return -y/(x+1); // ตัวอย่าง: f(x, y)
    }
    public static double euler(double x0, double y0, double h, double xTarget) {
        double y = y0;
        double x = x0;
        while (x < xTarget) {
            y = y + h * f(x, y);
            x += h;
        }
        return y;
    }
    public static double modifiedEuler(double x0, double y0, double h, double xTarget) {
        double y = y0;
        double x = x0;
        while (x < xTarget) {
            double k1 = f(x, y);
            double k2 = f(x + h, y + h * k1);
            y = y + h * (k1 + k2) / 2;
            x += h;
        }
        return y;
    }
    public static double rungeKutta(double x0, double y0, double h, double xTarget) {
        double y = y0;
        double x = x0;
        while (x < xTarget) {
            double k1 = h * f(x, y);
            double k2 = h * f(x + 0.5 * h, y + 0.5 * k1);
            double k3 = h * f(x + 0.5 * h, y + 0.5 * k2);
            double k4 = h * f(x + h, y + k3);
            y = y + (k1 + 2 * k2 + 2 * k3 + k4) / 6;
            x += h;
        }
        return y;
    }
}