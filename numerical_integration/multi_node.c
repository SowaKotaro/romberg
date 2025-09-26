#include <stdio.h>
#include <math.h>
#include <mpi.h>

// --------- 積分対象関数 ---------
double f_true(double x) {
    return sin(x) + cos(x);
}

// 解析的積分
double integral_true(double a, double b) {
    return (-cos(b) + sin(b)) - (-cos(a) + sin(a));
}

// --------- テイラー多項式 P_n(x) （中心 c=0） ---------
double taylor_poly(double x, int n) {
    double sum = 0.0;
    double term;
    int k;

    // cos(x) 展開
    for (k = 0; k <= n/2; k++) {
        term = pow(-1,k) * pow(x, 2*k) / tgamma(2*k+1);
        sum += term;
    }
    // sin(x) 展開
    for (k = 0; k <= (n-1)/2; k++) {
        term = pow(-1,k) * pow(x, 2*k+1) / tgamma(2*k+2);
        sum += term;
    }
    return sum;
}

// --------- 並列台形則積分 (MPI版) ---------
double trapezoidal_mpi(double (*func)(double, int), int n,
                       double a, double b, long long N,
                       int rank, int size, MPI_Comm comm) {
    double h = (b - a) / (double)N;
    long long chunk = N / size;
    long long start = rank * chunk + 1;
    long long end   = (rank == size - 1) ? (N - 1) : (start + chunk - 1);

    double local_sum = 0.0;
    for (long long i = start; i <= end; i++) {
        double x = a + i*h;
        local_sum += func(x, n);
    }

    // 端点処理は rank=0 がまとめて担当
    if (rank == 0) {
        local_sum += 0.5 * (func(a, n) + func(b, n));
    }

    double global_sum;
    MPI_Reduce(&local_sum, &global_sum, 1, MPI_DOUBLE,
               MPI_SUM, 0, comm);

    return global_sum * h;
}

int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    double a = 0.0, b = M_PI;
    double I_true = integral_true(a,b);

    if (rank == 0) {
        printf("真値 I = %.16f\n", I_true);
    }

    int taylor_list[] = {1,3,5,7,9,11,13,15};
    long long N_list[] = {1000,10000,100000,1000000};
    size_t n_t, n_N;

    for (n_t = 0; n_t < sizeof(taylor_list)/sizeof(taylor_list[0]); n_t++) {
        int n = taylor_list[n_t];
        for (n_N = 0; n_N < sizeof(N_list)/sizeof(N_list[0]); n_N++) {
            long long N = N_list[n_N];
            double approx = trapezoidal_mpi(taylor_poly, n, a, b, N,
                                            rank, size, MPI_COMM_WORLD);
            if (rank == 0) {
                double err = approx - I_true;
                printf("Taylor n=%2d, N=%10lld\n approx=%.16f, error=%.16f\n",
                       n, N, approx, err);
            }
        }
    }

    MPI_Finalize();
    return 0;
}
