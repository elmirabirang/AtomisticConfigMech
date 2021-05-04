#include "atom.h"

template<int dim> class ModifiedTersoff{
private:
    Atoms<dim> *atoms;

public:
    ModifiedTersoff();
    ~ModifiedTersoff();
    double ModiTers(int i, int j, int k);
};
