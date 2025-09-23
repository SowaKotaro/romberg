#include <stdio.h>
#include <math.h>

// 積分したい関数 f(x) = x^3 - 6x^2 + 9x
double f(double x) {
    return x*x*x - 6*x*x + 9*x;
//     return sin(1000*x) * exp(-1 * x * x);
}

// 台形公式による数値積分
double trapezoidal_integral(double a, double b, int n) {
    double h = (b - a) / n;
    double sum = 0.5 * (f(a) + f(b)); // 端の値を半分ずつ加算

    for (int i = 1; i < n; i++) {
        double x = a + i * h;
        sum += f(x);
    }

    return sum * h;
}

double calc_romberg(double nn, double n){
    double better = (4*nn - n) / 3;
    // printf("nn %.32f\n", nn);
    // printf(" n %.32f\n", n);
    // printf(" b %.32f\n", better);
    return better;
}

int main(void) {
    double a, b; // 積分区間 [a, b]
    int n;       // 分割数（大きいほど精度が高くなる）
    double ans = 950.0;
    a = 0;
    b = 10;
    n = 1000;
    int nn = n * 2;
    // 入力
    // printf("積分区間の下限 a を入力してください: ");
    // scanf("%lf", &a);
    // printf("積分区間の上限 b を入力してください: ");
    // scanf("%lf", &b);
    // printf("分割数 n を入力してください（例: 1000）: ");
    // scanf("%d", &n);

    // 台形積分の実行
    double result_n = trapezoidal_integral(a, b, n);



    // 結果の表示
    printf("result_n: %.16f\n", result_n);
    printf("     ans: %.16f\n", ans);
    printf("   diff.: %.16f\n\n", result_n - ans);
    // printf("f(x) = x^3 - 6x^2 + 9x の [%g, %g] における数値積分結果: %.32f\n", a, b, result);
    double result_nn = trapezoidal_integral(a, b, nn);
    double better = calc_romberg(result_nn, result_n);
    printf("original: %.16f\n", result_nn);
    printf("     ans: %.16f\n", ans);
    printf("   diff.: %.16f\n\n", result_nn - ans);
    printf("  better: %.16f\n", better);
    printf("     ans: %.16f\n", ans);
    printf("   diff.: %.16f\n", better - ans);
    
    return 0;
}
