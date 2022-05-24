#include <cstdio>
#include <fstream>

using namespace std;

#define ae 0.271
#define ai 0.999
#define ka 0.0000019
#define ro 0.000000001
#define beta 0.007
#define mu 0.0009
#define gamma 0
#define E0 93
#define R0 96
#define I0 0
#define D0 0

double deltaSensitive(double sensitive, double E, double infected, double recovered, double population) {
	return -1*(ai * sensitive * infected / population  + ae * sensitive * E / population) + gamma * recovered;
}

double deltaInfected(double E, double infected) {
	return ka * E - beta * infected - mu * infected;
}

double deltaE(double sensitive, double E, double infected, double population) {
	return (ai * sensitive * infected / population + ae * sensitive * E / population) - (ka + ro) * E;
}

double deltaRecovered(double E, double infected, double recovered) {
	return beta * infected + ro * E - gamma * recovered;
}

double deltaDead(double infected) {
	return mu * infected;
}

void euler_method(double start, double end, double step) {
	FILE *res = fopen("result.txt", "w"); // запись в файл 'result.txt'
	int n = (end - start) / step + 1; // количество промежутков
	double sensitive[2] = { 2798170 - E0 - R0 }; // восприимчивые (незараженные) индивидуумы c 3 лет
	double E[2] = { E0 }; // инфицированные без симптомов или находящиеся в инкубационном периоде индивидуумы
	double infected[2] = { I0 }; // инфицированные с симптомами
	double recovered[2] = { R0 }; // вылеченные
	double dead[2] = { D0 }; // умершие
	double population[2] = { sensitive[0] + E0 + I0 + R0 + D0 }; // вся популяция
	fprintf(res, "%.0lf %.2lf %.2lf %.2lf %.2lf %.2lf %.2lf\n", start, population[0], sensitive[0], E[0], infected[0], recovered[0], dead[0]); // запись в файл начальных значений

	// считаем значения и записываем их в ячейки массивов
	for (int i = 1; i < n; i++) {
		double t_prev = start + (i-1) * step;
		double t_next = start + i * step;
		int prev = (i-1) & 1;
		int next = i & 1;
		sensitive[next] = sensitive[prev] + step * deltaSensitive(sensitive[prev], E[prev], infected[prev], recovered[prev], population[prev]);
		E[next] = E[prev] + step * deltaE(sensitive[prev], E[prev], infected[prev], population[prev]);
		infected[next] = infected[prev] + step * deltaInfected(E[prev], infected[prev]);
		recovered[next] = recovered[prev] + step * deltaRecovered(E[prev], infected[prev], recovered[prev]);
		dead[next] = dead[prev] + step * deltaDead(infected[prev]);
		population[next] = sensitive[next] + E[next] + infected[next] + recovered[next] + dead[next];

		sensitive[next] = sensitive[prev] + step / 2 * (deltaSensitive(sensitive[prev], E[prev], infected[prev], recovered[prev], population[prev]) + deltaSensitive(sensitive[next], E[next], infected[next], recovered[next], population[next]));
		E[next] = E[prev] + step / 2 * (deltaE(sensitive[prev], E[prev], infected[prev], population[prev]) + deltaE(sensitive[next], E[next], infected[next], population[next]));
		infected[next] = infected[prev] + step / 2 * (deltaInfected(E[prev], infected[prev]) + deltaInfected(E[next], infected[next]));
		recovered[next] = recovered[prev] + step / 2 * (deltaRecovered(E[prev], infected[prev], recovered[prev]) + deltaRecovered(E[next], infected[next], recovered[next]));
		dead[next] = dead[prev] + step / 2 * (deltaDead(infected[prev]) + deltaDead(infected[next]));
		population[next] = sensitive[next] + E[next] + infected[next] + recovered[next] + dead[next];

		fprintf(res, "%.0lf %.2lf %.2lf %.2lf %.2lf %.2lf %.2lf\n", t_next, population[next], sensitive[next], E[next], infected[next], recovered[next], dead[next]); // запись в файл вычисленных значений
	}
	fclose(res); // закрываем файл 'result.txt'
}
