#include <vector>
//---
#include "atom.h"

using namespace std;

template<int dim> void writeDataFile(vector <Atom <dim> *> atoms, int load_step, string path);
