#include <vector>
//---
#include "atom.h"
#include "boundary.h"

template <int dim> float triangleArea(double xA, double yA, double xB, double yB, double xC, double yC);
template <int dim> void AssignRegion(Atoms<dim> *atoms, std::vector<Boundary<dim>> boundaries);
