#include <vector>
//---
#include "atom.h"
#include "boundary.h"

using namespace std;

template <int dim> float triangleArea(double xA, double yA, double xB, double yB, double xC, double yC);
template <int dim> void AssignRegion(Atom <dim>* &atom, vector<Boundary <dim>> boundaries);
