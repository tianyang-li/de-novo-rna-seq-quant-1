#include <iostream>
#include <algorithm>
#include <NLF.h>
#include <OptNewton.h>
#include <LinearEquation.h>
#include <CompoundConstraint.h>
#include <BoundConstraint.h>
#include <Constraint.h>

using namespace std;
using namespace NEWMAT;
using namespace OPTPP;

void init_test_calc(int ndim, ColumnVector& x) {
}

void test_calc(int mode, int ndim, const ColumnVector& x, double& fx,
		ColumnVector& gx, SymmetricMatrix& Hx, int& result) {
}

int main() {
	for (int i = 0; i != 1000000; ++i) {
		int ndim = 2;

		Matrix sum_eq1_mat(1, ndim);
		Real *all_1 = new Real[ndim];
		fill(all_1, all_1 + ndim, 1);
		for (int i = 0; i != ndim; ++i) {
			sum_eq1_mat << all_1;
		}
		delete[] all_1;

		ColumnVector sum_eq1_vec(1);
		sum_eq1_vec << 1;

		ColumnVector lower(ndim);
		Real *all_0 = new Real[ndim];
		fill(all_0, all_0 + ndim, 0);
		for (int i = 0; i != ndim; ++i) {
			lower << all_0;
		}
		delete[] all_0;

		Constraint geq0(new BoundConstraint(ndim, lower));
		Constraint sum_eq1(new LinearEquation(sum_eq1_mat, sum_eq1_vec));
		CompoundConstraint *constraints = new CompoundConstraint(geq0, sum_eq1);

		NLF2 calc(ndim, test_calc, init_test_calc, constraints);

		calc.initFcn();

		delete constraints;
	}

	return 0;
}

