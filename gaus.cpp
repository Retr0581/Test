#include <iostream>
#include <vector>
#include <random>
#include <ctime>
#include <cmath>
using namespace std;
const double EPS = 1e-12;

void show (vector<double> A, int n, int N){
    printf("\n");
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < N; j++) {
            printf("%lf ",A[i * N + j]);
        }
        printf("\n");
    }
}

void gaus (vector<double>& B, int n, int N){
    for (int i = 0; i < n; ++i) {
        int pivot_row = i;
        double max_val = fabs(B[i * N + i]);

        for (int k = i + 1; k < n; k++) {
            double val = fabs(B[k * N + i]);
            if (val > max_val) {
                max_val = val;
                pivot_row = k;
            }
        }

        if (max_val < EPS){
            cout << "Матрица вырожденная.\n";
            return;
        }

        if (pivot_row != i) {
            for (int j = i; j < N; j++) {
                double temp = B[i * N + j];
                B[i * N + j] = B[pivot_row * N + j];
                B[pivot_row * N + j] = temp;
            }
        }

        double inv_pivot = B[i * N + i];
        for (int j = i; j < N; j++) {
            B[i * N + j] /= inv_pivot;
        }

        for (int k = 0; k < n; k++) {
            double factor = B[k * N + i];
            if ((k != i) && (fabs(factor) > EPS)){
                for (int j = i; j < N; j++) {
                    B[k * N + j] -= factor * B[i * N + j];
                }
            }
        }
    }
}

int main() {
    int n;
    cout << "Введите размерность матрицы: ";
    cin >> n;

    if (n <= 0) {
        cout << "Ошибка: размерность должна быть > 0\n";
        return 1;
    }

    cout << "Генерация матрицы " << n << "x" << n << "...\n";
    vector<double> A(n * n);
    
    // делаем матрицу nxn
    mt19937 gen(time(nullptr));
    uniform_int_distribution<int> dist(-10, 10);

    for (int i = 0; i < n; i++) {
        double row_sum = 0.0;
        for (int j = 0; j < n; j++) {
            //if (i != j) {
            A[i * n + j] = dist(gen);
            //row_sum += fabs(A[i * n + j]);
           // }
        }
        //A[i * n + i] = row_sum + n * 10.0;
    }
    
    cout << "\nИсходная матрица (произвольная):\n";
    //вывод в консоль
    show (A, n, n);

    int N = n + 1;
    vector<double> B(n * N, 0.0);
    /*
    // делаем расширенную (A|b)
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            B[i * N + j] = A[i * n + j];
        }
        B[i * N + n] = 1.0;
    }
    cout << "\nРасширенная матрица (A|b):\n";
    show (B, n, N);

    */

    // делаем расширенную матрицу (A|b)
    for (int i = 0; i < n; i++) {
        double row_sum = 0.0;
        for (int j = 0; j < n; j++) {
            B[i * N + j] = A[i * n + j];
            row_sum += A[i * n + j];  // Сумма элементов строки
        }
        B[i * N + n] = row_sum;  // b_i = сумма строки i
    }
    cout << "\nРасширенная матрица (A|b):\n";
    show(B, n, N);

    // делаем вырожденную матрицу для проверки
    vector<double> C = B;
    for (int i = 0; i < n; i++) {
        C[i * N] = 0.0;
    }





    // Жордан-Гаусс
    cout << "\nМетод Гауса с выбором ведущего эл-та...\n";
    clock_t start = clock();
    
    gaus (B, n, N);

    clock_t finish = clock();
    double t = (finish - start) / (double)CLOCKS_PER_SEC;

    show (B, n, N);

    cout << "\nВычисление завершено.\n";
    cout << "Время выполнения: " << fixed << t << " секунд\n";
    cout << "\nНайденное решение:\n";
    // достаем решение системы из правого столбца
    vector<double> x(n);
    for (int i = 0; i < n; i++) {
            x[i] = B[i * N + n];
            printf("x%d = %lf\n", i+1, x[i]);
    }


    cout << "\nПроверка правильности...\n";
    double max_error = 0.0;
    for (int i = 0; i < n; i++) {
        double sum = 0.0;
        for (int j = 0; j < n; j++) {
            sum += A[i * n + j] * x[j];
        }
        double expected = 1.0;  // потому что b[i] = 1.0
        double error = fabs(sum - expected);
        if (error > max_error) {
            max_error = error;
        }
    }
    cout << "Результат проверки: " << scientific << max_error << endl;

    cout << "\nВырожденная версия исходной матрицы:\n";
    show (C, n, N);
    cout << "\nПробуем ее решить...\n";
    gaus (C, n, N);
    return 0;
}
