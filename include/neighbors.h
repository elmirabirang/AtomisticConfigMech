#include <vector>
//---
#include "atom.h"

template<int dim> vector<Atom <dim> *> FindNeighbors(vector <Atom <dim> *> &unrelax_atoms, double cut_radius);
