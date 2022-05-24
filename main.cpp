#include <iostream>
#include "euler.h"

using namespace std;

int main() {
	double start = 0; // начало отрезка
	double end = 90; // конец отрезка
	double step = 1; // шаг

	cout << "Input data:" << endl
			<< "Start of the segment: " << start << endl 
			<< "End of the segment: " << end << endl 
			<< "Calculation step: " << step << endl; 
	euler_method(start, end, step);
	cout << "The output data is written to a file 'result.txt'" << endl;
	
	return 0;
}