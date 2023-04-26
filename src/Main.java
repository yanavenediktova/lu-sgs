import static java.lang.Math.*;

public class Main {

    public static void main(String[] args) {
        double Lx = 1;
        double Ly = 1;
        double hx = 0.01;
        double hy = 0.01;
        double kx = 0.2;
        double ky = 0.3;
        double T = 0.0;
        double tau = 0.0001;
        double u0 = 10.0;
        double t = 0.0;

        int Nx = (int) (Lx / hx) + 1;
        int Ny = (int) (Ly / hy) + 1;
        double U[][] = new double[Nx][Ny];

        for (int i = 0; i < Nx; i++) {
            for (int j = 0; j < Ny; j++) {
                U[i][j] = u0;
            }
        }
        for (int i = 0; i < Nx; i++) {
            U[i][0] = 100.0; // Граничное условие u(x,0) = 0
            U[i][Ny - 1] = 0.0; // Граничное условие u(x,Ly) = 0
        }
        for (int j = 0; j < Ny; j++) {
            U[0][j] = 0.0; // Граничное условие u(0,y) = 0
            U[Nx - 1][j] = 10.0; // Граничное условие u(Lx,y) = 0
        }
        double iter = 0.0;
        double norm;
        double eps = 1.e-15;
        do {
            t += tau;

            double L1 = tau * (kx + ky) / (2 * hx * hx);
            double L2 = tau * (kx + ky) / (2 * hy * hy);
            double D = 1 + (2 * Lx + 2 * Ly);
            double U1 = L1;
            double U2 = L2;
            double y[][] = new double[Nx][Ny];
            double x[] = new double[Nx * Ny];
            double sum;

            for (int i = 0; i < Nx; i++) {
                for (int j = 0; j < Ny; j++) {
                    y[i][j] = U[i][j];
                }
            }
            // do {
            iter++;
            norm = 0;
            sum = 0;
            for (int i = 1; i < Nx - 1; i++) {
                // Прямой ход
                for (int j = 1; j < Ny - 1; j++) {
                    sum = y[i][j];
                    sum += L1 * y[i - 1][j] + L2 * y[i][j - 1];
                    y[i][j] = sum / D;
                }
            }
            // Обратный ход
            for (int i = Nx - 2; i > 0; i--) {
                for (int j = Ny - 2; j > 0; j--) {
                    sum = y[i][j];
                    sum += (U1 * U[i + 1][j] + U2 * U[i][j + 1]) / D;
                    y[i][j] = sum;
                    if (norm < (abs(sum - U[i][j]))) {
                        norm = (abs(sum - U[i][j]));
                    }
                }

            }
            double sum1 = 0;
            double sum2 = 0;
            for (int i = Nx - 2; i > 0; i--) {
                for (int j = Ny - 2; j > 0; j--) {
                    sum1 +=y[i][j];
                    sum2 +=U[i][j];
                }
            }
            norm = sqrt(pow(abs(sum1-sum2),2));
            // Обновление значений температуры на сетке
            for (int i = 1; i < Nx - 1; i++) {
                for (int j = 1; j < Ny - 1; j++) {
                    U[i][j] = y[i][j];
                }
            }
            // } while (norm >= eps);
        } while (norm >= eps);
        for (int i = 0; i < Nx; i++) {
            for (int j = 0; j < Ny; j++) {
                System.out.print(U[i][j] + " ");
            }
            System.out.println(" ");
        }
        System.out.println(iter);
    }
}
