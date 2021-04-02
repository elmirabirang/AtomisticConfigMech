#include "atom.h"
#include "bond.h"
#include "point.h"
//#include "matrixx.h"
#include "modified_tersoff.h"
#include "modified_tersoff.cpp"

int main()
{
    Atom <3> atom_alpha;
    Atom <3> atom_beta;
    Atom <3> atom_gamma;
    Point<3> atom_alpha_position(0,0,0);
    Point <3> atom_beta_position(0.65, 0.8, 0);
    Point <3> atom_gamma_position(1.3, 0, 0);
    int atom_alpha_id=0;
    int atom_beta_id=1;
    int atom_gamma_id=2;
    int atom_cell_id=-1;

    atom_alpha.SetSpatialPosition(atom_alpha_position);
    atom_alpha.SetID(atom_alpha_id);
    atom_alpha.SetCellID(atom_cell_id);

    atom_beta.SetSpatialPosition(atom_beta_position);
    atom_beta.SetID(atom_beta_id);
    atom_beta.SetCellID(atom_cell_id);

    atom_gamma.SetSpatialPosition(atom_gamma_position);
    atom_gamma.SetID(atom_gamma_id);
    atom_gamma.SetCellID(atom_cell_id);

    ModifiedTersoff<3> ModT;
    //ModT.ModiTers(atom_alpha, atom_beta, atom_gamma);



   return 0;
}
