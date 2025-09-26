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
        // tgamma(2*k+1) は (2*k)! に相当
        term = pow(-1,k) * pow(x, 2*k) / tgamma(2*k+1);
        sum += term;
    }
    // sin(x) 展開
    for (k = 0; k <= (n-1)/2; k++) {
        // tgamma(2*k+2) は (2*k+1)! に相当
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
    // 積分区間 [a, b] を N 等分したとき、i*h の i は 1 から N-1 まで
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

    // ---------------- 測定開始 ----------------
    double start_time = MPI_Wtime();

    // int taylor_list[] = {1,3,5,7,9,11,13,15,17,19,21,23,25,27};
    // long long N_list[] = {10,50,100,500,1000,5000,10000,50000,100000,500000,1000000,5000000,10000000, 50000000, 100000000};
    int taylor_list[] = {25};
    long long N_list[] = {100000000};
    size_t n_t, n_N;

    for (n_t = 0; n_t < sizeof(taylor_list)/sizeof(taylor_list[0]); n_t++) {
        int n = taylor_list[n_t];
        for (n_N = 0; n_N < sizeof(N_list)/sizeof(N_list[0]); n_N++) {
            long long N = N_list[n_N];
            // trapezoidal_mpi内でMPI_Reduceが行われるため、
            // rank 0 の approx の値が全体の積分結果となる
            double approx = trapezoidal_mpi(taylor_poly, n, a, b, N,
                                             rank, size, MPI_COMM_WORLD);
            if (rank == 0) {
                double err = approx - I_true;
                printf("Taylor n=%2d, N=%10lld: approx=%.16f, error=%.16f\n",
                       n, N, approx, err);
            }
        }
    }

    // ---------------- 測定終了 ----------------
    double end_time = MPI_Wtime();
    double elapsed_time = end_time - start_time;

    // 全プロセスの実行時間を集計し、rank 0 で表示
    // 全体の実行時間は、最も遅いプロセス（全ての計算と通信を含む）の時間を採用すべきだが、
    // ここでは単純に rank 0 の時間を表示する（これは計算と通信が比較的均等なら全体の時間に近くなる）
    if (rank == 0) {
        printf("\n全計算時間: %.4f 秒\n", elapsed_time);
    }

    MPI_Finalize();
    return 0;
}
