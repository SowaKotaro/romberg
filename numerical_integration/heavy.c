#include <stdio.h>
#include <math.h>

// --- 真値を持つ関数（多項式） ---
double f_true(double x) {
    return x*x*x - 6*x*x + 9*x; // 解析的に積分可能
}

// --- 遅い関数（計算用） ---
double f(double x) {
    double y = f_true(x); // 元の関数ベース
    // 重い処理を追加して意図的に遅くする
    for (int i = 0; i < 500; i++) {
        y = sin(y + x*i) + cos(y - x*i) + exp(-fabs(y)) + log(fabs(y) + 1.0);
        y = sqrt(fabs(y) + 1.0);
    }
    return y;
}

// --- 台形公式 ---
double trapezoidal_integral(double a, double b, long long n) {
    double h = (b - a) / (double)n;
    double sum = 0.5 * (f(a) + f(b));
    for (long long i = 1; i < n; i++) {
        double x = a + i * h;
        sum += f(x);
    }
    return sum * h;
}

// --- 解析的積分値（真値） ---
// f_true(x) = x^3 - 6x^2 + 9x
// ∫ f_true dx = (1/4)x^4 - 2x^3 + (9/2)x^2
double integral_true(double a, double b) {
    double F(double x) { return 0.25*x*x*x*x - 2*x*x*x + 4.5*x*x; }
    return F(b) - F(a);
}

int main(void) {
    double a = 0.0, b = 10.0;
    long long n = 10000000; // ←ここを大きくすると爆発的に時間がかかる

    // 真値
    double ans = integral_true(a, b);

    // 数値積分
    double result = trapezoidal_integral(a, b, n);

    // 出力
    printf("区間: [%g, %g]\n", a, b);
    printf("分割数: %lld\n", n);
    printf("数値積分: %.16f\n", result);
    printf("真値      : %.16f\n", ans);
    printf("誤差      : %.16f\n", result - ans);

    return 0;
}
