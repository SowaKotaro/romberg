#include <stdio.h>
#include <math.h>

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

// --------- 台形則積分 ---------
double trapezoidal(double (*func)(double, int), int n, double a, double b, long long N) {
    double h = (b - a) / (double)N;
    double sum = 0.5 * (func(a, n) + func(b, n));
    for (long long i = 1; i < N; i++) {
        double x = a + i*h;
        sum += func(x, n);
    }
    return sum * h;
}

int main(void) {
    double a = 0.0, b = M_PI;
    double I_true = integral_true(a,b);
    printf("真値 I = %.16f\n", I_true);

    int taylor_list[] = {1,3,5,7,9,11,13,15,17,19,21,23};
    long long N_list[] = {10,50,100,500,1000,5000,10000,50000,100000,500000,1000000,5000000,10000000};
    size_t n_t, n_N;

    for (n_t = 0; n_t < sizeof(taylor_list)/sizeof(taylor_list[0]); n_t++) {
        int n = taylor_list[n_t];
        for (n_N = 0; n_N < sizeof(N_list)/sizeof(N_list[0]); n_N++) {
            long long N = N_list[n_N];
            double approx = trapezoidal(taylor_poly, n, a, b, N);
            double err = approx - I_true;
            printf("Taylor n=%2d, N=%5lld\n approx=%.16f, error=%.16f\n",
                   n, N, approx, err);
        }
    }

    return 0;
}
