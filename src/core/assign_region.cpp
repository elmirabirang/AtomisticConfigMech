/*
 * assign_region.cpp
 *
 *  Created on: Dec 6, 2018
 *      Author: S.Elmira Birang.O
 *      Mechanical Engineering PhD Candidate
 *      Chair of Applied Mechanics
 *      University of Erlangen-Nürnberg
 */

#include <cmath>
#include <iostream>
#include <vector>
//---
#include "atom.h"
#include "boundary.h"
#include "point.h"

template<int dim> float triangleArea(double xA, double yA, double xB, double yB, double xC, double yC) {
   return abs((xA * (yB - yC) + xB * (yC - yA) + xC * (yA - yB)) / 2.0);
}

template<int dim> void AssignRegion(Atoms<dim> *atoms, std::vector<Boundary <dim>> boundaries) {
    LOOP_OVER_ATOMS(atoms, i,
        Point<dim> material_position = atoms->getMaterialPosition(i);
        atoms->setRegion(i, 0);

        if(dim == 2) {
            for(int b = 0; b < boundaries.size(); ++b) {
                int bid = boundaries[i].GetBid();
                if (bid < 10 && b == 0) {
                    double x_zone_min = boundaries[b].GetBoundaryX0();;
                    double y_zone_min = boundaries[b].GetBoundaryY0();;
                    double x_zone_max = boundaries[b].GetBoundaryX1();
                    double y_zone_max = boundaries[b].GetBoundaryY1();
                    double xa = material_position.GetXCoord();
                    double ya = material_position.GetYCoord();

                    if(xa > x_zone_min && xa < x_zone_max && ya > y_zone_min && ya < y_zone_max) {
                        atoms->setZone(i, 1);
                    } else {
                        atoms->setZone(i, 0);
                    }
                } else if (bid < 10 && b != 0) {
                    double x0 = boundaries[b].GetBoundaryX0();
                    double y0 = boundaries[b].GetBoundaryY0();
                    double x1 = boundaries[b].GetBoundaryX1();
                    double y1 = boundaries[b].GetBoundaryY1();
                    double xa = material_position.GetXCoord();
                    double ya = material_position.GetYCoord();

                    if(xa >= x0 && xa <= x1 && ya >= y0 && ya <= y1) {
                        atoms->setRegion(i, bid);
                        break;
                    }
                } else {
                    double xA = boundaries[b].GetBoundaryX0();
                    double yA = boundaries[b].GetBoundaryY0();
                    double xB = boundaries[b].GetBoundaryX1();
                    double yB = boundaries[b].GetBoundaryY1();
                    double xC = boundaries[b].GetBoundaryX2();
                    double yC = boundaries[b].GetBoundaryY2();
                    double xP = material_position.GetXCoord();
                    double yP = material_position.GetYCoord();
                    float area = triangleArea<dim>(xA, yA, xB, yB, xC, yC);  //area of triangle ABC
                    float area1 = triangleArea<dim>(xP, yP, xB, yB, xC, yC); //area of PBC
                    float area2 = triangleArea<dim>(xA, yA, xP, yP, xC, yC); //area of APC
                    float area3 = triangleArea<dim>(xA, yA, xB, yB, xP, yP); //area of ABP
                    if(roundf(area) == roundf(area1 + area2 + area3)) {
                        atoms->setRegion(i, bid);
                        break;
                    }
                }
            }
        } else {
            for(int b = 0; b < boundaries.size(); ++b) {
                int bid = boundaries[b].GetBid();

                if(bid < 10) {
                    double x0 = boundaries[b].GetBoundaryX0();
                    double y0 = boundaries[b].GetBoundaryY0();
                    double z0 = boundaries[b].GetBoundaryZ0();
                    double x1 = boundaries[b].GetBoundaryX1();
                    double y1 = boundaries[b].GetBoundaryY1();
                    double z1 = boundaries[b].GetBoundaryZ1();
                    double xa = material_position.GetXCoord();
                    double ya = material_position.GetYCoord();
                    double za = material_position.GetZCoord();

                    if(xa >= x0 && xa <= x1 && ya >= y0 && ya <= y1 && za >= z0 && za <= z1) {
                        atoms->setRegion(i, bid);
                        break;
                    }
                }
            }
        }
    )
}

template void AssignRegion<2>(Atoms<2> *atoms, std::vector<Boundary<2>> boundaries);
template void AssignRegion<3>(Atoms<3> *atoms, std::vector<Boundary<3>> boundaries);
