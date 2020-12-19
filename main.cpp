#include <iostream>
#include <random>
#include <map>
#include <algorithm>
#include <ctime>
#include <iomanip>

using namespace std;

const double a = -2.;
const double b = 2.;
const double c = -0.5;
const double d = 0.;
const int N = 16;
const double A = 2.;
const double step = (b - a) / (1. * (N - 1));

double myfun(const double &c_coef, const double &d_coef, const double &x) {
    return c_coef * x + d_coef;
}

double Random(double low_limit, double high_limit)
{
    return low_limit + (1.*rand()/RAND_MAX)*(high_limit - low_limit);   // (1.*rand()/RAND_MAX)-генерация числа от 0 до 1
}

auto XandYMapWithNoise(const double &low_limit, const double &noise) {
    map<double, double> XandYmap;
    for (int i = 0; i < N; ++i) {
        double x = low_limit + i * step;
        double y = myfun(c, d, low_limit + i * step) + noise * Random(-0.5, 0.5);
        XandYmap[x] = y;
    }
    return XandYmap;
}

double SumErrors(const map<double, double> &XandYmap, const double &w1, const double &w0) {
    double sum = 0;
    for (auto p : XandYmap) {
        sum += (myfun(w1, w0, p.first) - p.second) * (myfun(w1, w0, p.first) - p.second);
    }
    return sum;
}

int N_count(double q, double P) {
    return (ceil(log(1. - P) / log(1. - q)));
}

double PassiveSearchForW1(const map<double, double> &XandYmap, const double &low_limit, const double &high_limit) { //search for w1
    int i = 1;
    const float eps = 0.01;
    double current_error = 1;
    double x;
    double minx = 0;
    double w0 = 0;
    while (current_error > eps){
        double min = 1.7976931348623158*pow(10, 308);
        for (int k = 1; k < i + 1; k++) {
            x = ((high_limit - low_limit) * k / (i + 1)) + low_limit;
            if (SumErrors(XandYmap,x, w0) < min){
                min = SumErrors(XandYmap,x, w0);
                minx = x;
            }
        }
        current_error = (high_limit - low_limit) / (i + 1);
        i++;
    }
    return minx;
}

double RandomSearchFowW0(const map<double, double> &XandYmap, const double &low_limit, const double &high_limit,
                         const double &w1){
    int z;
    double MAX_DOUBLE = 1.7976931348623158*pow(10, 308);
    double x;
    double q = 0.005;
    double p = 0.9;
    double min;
    double minx;
    for (int i = 0; i < 20; i++) {
        for (int j = 0; j < 10; j++) {
            z = N_count(q, p);
            min = MAX_DOUBLE;
            for (auto k = 0; k < z; k++) {
                x = Random(low_limit, high_limit);
                if (SumErrors(XandYmap, w1, x) < min){
                    min = SumErrors(XandYmap, w1, x);
                    minx = x;
                }
            }
        }
    }
    return minx;
}

void Print(const double &noise) {
    map<double, double> XandYmap = XandYMapWithNoise(a, noise);
    size_t i = 1;
    cout << "_____________________________________\n"
            "|   Номер   | Значение  | Значение  |\n"
            "|   точки   | аргумента | функции   |\n"
            "|    (i)    |   (x)     |    (y)    |\n"
            "-------------------------------------" << endl;
    for (auto p : XandYmap) {
        cout << "|" << setw(11) << i << "|" << setw(11) << p.first << "|"
             << setw(11) << p.second << "|" << endl;
        ++i;
    }
    cout << "-------------------------------------" << endl;

    double Cmin = -1, Cmax = 0;
    double Dmin = -1, Dmax = 1;
    cout << "Cmin = " << Cmin << "\nCmax = " << Cmax << "\nDmin = " << Dmin << "\nDmax = " << Dmax << "\n";
    double w1 = PassiveSearchForW1(XandYmap, Cmin, Cmax);
    cout << "c = " << w1 << "\n";
    double w0 = RandomSearchFowW0(XandYmap, Dmin, Dmax, w1);
    cout << "d = " << w0 << "\n";
    cout << "Error = " << SumErrors(XandYmap, w1, w0) << "\n";
}

int main() {
    srand(time(NULL));
    cout << "Результаты для функции y = " << c << " * x + " << d <<
        " на интервале [ " << a << " , " << b << " ]:\n";
    Print(0.);
    cout <<"Результаты для функции y = " << c << " * x + " << d << "с шумом A = " << A <<
        " на интервале [ " << a << " , " << b << " ]:\n";
    Print(A);
}