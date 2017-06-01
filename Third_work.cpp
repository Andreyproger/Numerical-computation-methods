#define _USE_MATH_DEFINES
#include <cmath>
#include<iostream>
#include <iomanip>
#include <fstream>
#include <vector>
using namespace std;
//START OF SHOOTING
double rungeZ(double y, double z, double h) {
double l1 = h * (y + 10.6);
double l2 = h * (y + 10.987 + l1);
double z1 = z + (l1 + l2) / 2;
return z1;
}
double rungeY(double y, double z, double h) {
double l1 = h * (y + 10.6);
double k1 = h * z;
double k2 = h * (z + l1);
double y1 = y + (k1 + k2) / 2;
return y1;
}
double adamsZ(double x,  double y0, double y1,
double  z, double h) {
double zn = z + h/2 * ( 3 * y1 - y0 + 21.2
+ 4.3*(3*x - 3*x*x - (x - h) + (x-h)*(x-h)) );
return zn;
}
double adamsY(double y, double z0, double z1,
double h) {
double yn = y + h/2*( 3*z1 - z0 );
return yn;
}
double newton(double mu, double y, double z) {
double mun = mu - (z +y - 2*exp(1) - 2.3)
/ (y + z + 10.6);
return mun;
}
//END OF SHOOTING
//START OF SHUTTLE
void   calcMatrix(int   n,   std::vector<double>&
aOut,   std::vector<double>&   bOut,
std::vector<double>& cOut) {
double h = 1.0 / n;
aOut.resize(n + 1);
bOut.resize(n + 1);
cOut.resize(n + 1);
aOut[0] = 0;
bOut[0] = - (h * h + 2 * h + 2);
cOut[0] = 2;
aOut[n] = 2;
bOut[n] = - (2 + 2*h + h * h);
cOut[n] = 0;
for(int i = 1; i <= n - 1; i++) {
aOut[i] = 1;
bOut[i] = - (2 + h * h);
cOut[i] = 1;
}
}
void   calcVector(int   n,   std::vector<double>&
dOut) {
double h = 1.0 / n;
double e = exp(1);
dOut.resize(n + 1);
dOut[0] = 2 * (1 + 4.3) * h * h - 2 * 4.3 *
h;
dOut[n] = 2 * (1 + 4.3) * h * h - 2 * h *
(2*e + 2.3);
for (int i = 1; i <= n - 1; ++i) {
dOut[i] = (2 * 4.3 + 2 + 4.3 * h*i *
(1.0 - h*i)) * h * h;
}
}
double nextLambda(double a, double b, double c,
double lambda) {
return - c / (a * lambda + b);
}
double   nextNu(double   a,   double   b,   double   c,
double d, double lambda, double nu) {
return (d - a * nu) / (a * lambda + b);
}
double prevY(double y, double lambda, double
nu) {
return lambda * y + nu;
}
//END OF SHUTTLE
7
int main() {
double eps = 0.000001;
double n = 20.0;
double h = 1 / n;
double resY[21];
double resZ[21];
double mu0;
double mu1 = 4.3;
ofstream fout("res.txt");
do {
mu0 = mu1;
resY[0] = mu0;
resZ[0] = mu0 - 4.3;
resZ[1] = rungeZ(resY[0], resZ[0],
h); 
resY[1]   =   rungeY(resY[0],
resZ[0], h);
for(int i = 2; i < 21; i++) {
resZ[i]   =   adamsZ(h*i,
resY[i-2], resY[i-1], resZ[i-1], h);
resY[i]   =   adamsY(resY[i-
1], resZ[i-2], resZ[i-1], h);
}
mu1   =   newton(mu0,   resY[20],
resZ[20]);
cout << mu1 << "\n";
} while (abs(mu1 - mu0) > eps);
//================
std::vector<double> lambda(n + 2);
    std::vector<double> nu(n + 2);
std::vector<double> y(n + 1);
std::vector<double> a, b, c, d;
calcMatrix(n, a, b, c);
calcVector(n, d);
lambda[0] = 0;
nu[0] = 0;
for (int i = 0; i < n + 1; ++i) {
lambda[i + 1] = nextLambda(a[i],
b[i], c[i], lambda[i]);
nu[i + 1] = nextNu(a[i], b[i], c[i],
d[i], lambda[i], nu[i]);
}
y[n] = nu[n + 1];
for (int i = n; i >= 1; --i) {
y[i   -   1]   =   prevY(y[i],   lambda[i],
nu[i]);
}
//
fout << "x" << " = [";
for(int i = 0; i <= n; i++) {
fout   <<   fixed   <<   setprecision(5)
<< h*i << "; ";
}
fout << "] \n";
//
fout << "z1" << " = [";
for(int i = 0; i <= n; i++) {
fout   <<   fixed   <<   setprecision(5)
<< resY[i] << "; ";
}
fout << "] \n";
//
fout << "z2" << " = [";
for(int i = 0; i <= n; i++) {
fout   <<   fixed   <<   setprecision(5)
<< y[i] << "; ";
}
fout << "] \n";
//
fout << "t1" << " = [";
for(int i = 0; i <= n; i++) {
fout   <<   fixed   <<   setprecision(5)
<< 
4.3*i*h*i*h-4.3*i*h+exp(-i*h)
+exp(i*h)-2 << "; ";
}
fout << "] \n";
//
fout << "figure\nplot(x" << ",t1,'-ro'," <<
"x" << ",z1,'-.b'," <<
"x" << ",z2,'-*g'," <<
")   \ntitle('[0.0;   1.0],   h   =   "   <<   h   <<   "')
\nlegend('Точное   решение','Метод
стрельбы','Метод
прогонки','location','northwest') \n";
fout.close();
}
