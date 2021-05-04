#include <vector>
//---
#include "atom.h"

template<int dim> Atoms<dim> *UnrelaxedConfigGenerator(
    int number_x, int number_y, int number_z,
    double lattice_const, int dimension, string type);
