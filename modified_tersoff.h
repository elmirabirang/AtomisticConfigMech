#include "atom.h"
#include "bond.h"
#include "point.h"
#include "math.h"
//#include "matrixx.h"
#include <vector>
#include <string>

using namespace std;

template <int dim>
class ModifiedTersoff{
	private:

		Atom <dim> atomi;
		Atom <dim> atomj;
		Atom <dim> atomk;

	public:

		ModifiedTersoff();
		~ModifiedTersoff();

		double ModiTers(Atom <dim> &atom_Alpha, Atom <dim> &atom_beta, Atom <dim> &atom_Gamma);


};

