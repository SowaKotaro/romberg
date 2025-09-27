#include <stdio.h>
#include <math.h>

// 25次でのsin(x)+cos(x)のテイラー展開
double f_taylor(double x) {
    double sum = 1.0; // 0次の項
    double term = 1.0;
    for (int i = 1; i <= 25; i++) {
        term = term * x / i;
        int sign_selector = i / 2;
        if (sign_selector % 2 != 0) {
             sum -= term;
        } else {
             sum += term;
        }
    }
    return sum;
}

#define TRUE_VALUE 2.0
#define EPS 1e-15
#define MAX_LEVEL 20   // 念のための最大反復回数

// --------- 被積分関数 ---------
double f(double x) {
    return sin(x) + cos(x);
}

// --------- 台形公式 ---------
// n分割台形公式の値を返す
double trapezoidal(double a, double b, int n) {
    double h = (b - a) / n;
    double sum = 0.5 * (f_taylor(a) + f_taylor(b));
    for (int i = 1; i < n; i++) {
        sum += f_taylor(a + i * h);
    }
    return sum * h;
}

int main(void) {
    double a = 0.0, b = M_PI;
    double R[MAX_LEVEL][MAX_LEVEL];  // ロンバーグ表
    int n = 1; // 初期分割数
    int k, j;

    // 最初の値（台形公式）
    R[0][0] = trapezoidal(a, b, n);
    printf("R[0][0] (n=%d) = %.16f, error = %.15f\n", n, R[0][0], fabs(R[0][0] - TRUE_VALUE));

    for (k = 1; k < MAX_LEVEL; k++) {
        n *= 2;
        R[k][0] = trapezoidal(a, b, n);  // 台形則の計算
        printf("R[%d][0] (n=%d) = %.16f, error = %.15f\n", k, n, R[k][0], fabs(R[k][0] - TRUE_VALUE));

        // 加速処理
        for (j = 1; j <= k; j++) {
            R[k][j] = (pow(4, j) * R[k][j-1] - R[k-1][j-1]) / (pow(4, j) - 1);
            printf("R[%d][%d] = %.16f, error = %.15f\n", k, j, R[k][j], fabs(R[k][j] - TRUE_VALUE));
        }

        // 精度判定
        if (fabs(R[k][k] - TRUE_VALUE) < EPS) {
            printf("\n収束しました: R[%d][%d] = %.16f\n", k, k, R[k][k]);
            break;
        }
    }

    return 0;
}
