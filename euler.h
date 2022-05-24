#include <cstdio>
#include <iostream>
#include <fstream>
#include <string>

using namespace std;

#define ae 0.271
#define ai 0.999
#define ka 0.0000019
#define ro 0.000000009
#define beta 0.007
#define mu 0.0009
#define gamma 0
#define E0 93
#define R0 96
#define I0 0
#define D0 0

double func_sensitive(double sensitive, double E, double infected, double recovered, double population) {
	return -1*(ai * sensitive * infected / population  + ae * sensitive * E / population) + gamma * recovered;
}

double func_infected(double E, double infected) {
	return ka * E - beta * infected - mu * infected;
}

double func_E(double sensitive, double E, double infected, double population) {
	return (ai * sensitive * infected / population + ae * sensitive * E / population) - (ka + ro) * E;
}

double func_recovered(double E, double infected, double recovered) {
	return beta * infected + ro * E - gamma * recovered;
}

double func_dead(double infected) {
	return mu * infected;
}

void euler_method(double start, double end, double step) {
	FILE *res = fopen("result.txt", "w"); // запись в файл 'result.txt'
	int n = (end - start) / step + 1; // количество промежутков
	
    double *sensitive = new double[n]; // Восприимчивые (незараженные) индивидуумы c 3 лет;
    double *E = new double[n]; // Зараженные или находящиеся в инкубационном периоде индивидуумы;
    double *infected = new double[n]; // Инфицированные с симптомами;
    double *recovered = new double[n]; // Вылеченные;
    double *dead = new double[n]; // Умершие;
    double *population = new double[n]; // Вся популяция;
	
	sensitive[0] = 2798170 - E0 - R0 ;
	E[0] = E0;
	infected[0] = I0;
	recovered[0] = R0;
	dead[0] = D0;
	population[0] = sensitive[0] + E0 + I0 + R0 + D0; // вся популяция
	fprintf(res, "%.0lf %.2lf %.2lf %.2lf %.2lf %.2lf %.2lf\n", start, population[0], sensitive[0], E[0], infected[0], recovered[0], dead[0]); // запись в файл начальных значений
	
	// считаем значения и записываем их в ячейки массивов
	for (int i = 1; i < n; i++) {
		double t_next = start + i * step;
		sensitive[i] = sensitive[i - 1] + step * func_sensitive(sensitive[i - 1], E[i - 1], infected[i - 1], recovered[i - 1], population[i - 1]);
		E[i] = E[i - 1] + step * func_E(sensitive[i - 1], E[i - 1], infected[i - 1], population[i - 1]);
		infected[i] = infected[i - 1] + step * func_infected(E[i - 1], infected[i - 1]);
		recovered[i] = recovered[i - 1] + step * func_recovered(E[i - 1], infected[i - 1], recovered[i - 1]);
		dead[i] = dead[i - 1] + step * func_dead(infected[i - 1]);
		population[i] = sensitive[i] + E[i] + infected[i] + recovered[i] + dead[i];

		sensitive[i] = sensitive[i - 1] + step / 2 * (func_sensitive(sensitive[i - 1], E[i - 1], infected[i - 1], recovered[i - 1], population[i - 1]) + func_sensitive(sensitive[i], E[i], infected[i], recovered[i], population[i]));
		E[i] = E[i - 1] + step / 2 * (func_E(sensitive[i - 1], E[i - 1], infected[i - 1], population[i - 1]) + func_E(sensitive[i], E[i], infected[i], population[i]));
		infected[i] = infected[i - 1] + step / 2 * (func_infected(E[i - 1], infected[i - 1]) + func_infected(E[i], infected[i]));
		recovered[i] = recovered[i - 1] + step / 2 * (func_recovered(E[i - 1], infected[i - 1], recovered[i - 1]) + func_recovered(E[i], infected[i], recovered[i]));
		dead[i] = dead[i - 1] + step / 2 * (func_dead(infected[i - 1]) + func_dead(infected[i]));
		population[i] = sensitive[i] + E[i] + infected[i] + recovered[i] + dead[i];

		fprintf(res, "%.0lf %.2lf %.2lf %.2lf %.2lf %.2lf %.2lf\n", t_next, population[i], sensitive[i], E[i], infected[i], recovered[i], dead[i]); // запись в файл вычисленных значений
		delete[] sensitive;
		delete[] E;
		delete[] infected;
		delete[] recovered;
		delete[] dead;
		delete[] population;
	}
	fclose(res); // закрываем файл 'result.txt'
}
