/*
 * Stillinger_Weber.h
 *
 *  Created on: Oct 19, 2020
 *      Author: S. Elmira Birang O.
 *      Mechanical Engineering PhD Candidate
 *      University of Erlangen-Nuremberg
 */

#pragma once

#include "atom.h"
#include "point.h"

template<int dim> class StillingerWeber{
public:
    StillingerWeber();
    ~StillingerWeber();

    double StillingerWeberTwoBody_Deformational(
        int i, int j, double sigma_AlphaBeta, double epsilon_AlphaBeta,
        double A_AlphaBeta, double B_AlphaBeta, double I_AlphaBeta, double J_AlphaBeta, double a_AlphaBeta);

    double StillingerWeberThreeBody_Deformational(
        int i, int j, int k, double epsilon_AlphaBetaGamma, double eta_AlphaBetaGamma, double cosine_teta0);

    double StillingerWeberTwoBody_Configurational(
        int i, int j, double sigma_AlphaBeta, double epsilon_AlphaBeta,
        double A_AlphaBeta, double B_AlphaBeta, double I_AlphaBeta, double J_AlphaBeta, double a_AlphaBeta);

    double StillingerWeberThreeBody_Configurational(
        int i, int j, int k, double epsilon_AlphaBetaGamma, double eta_AlphaBetaGamma, double cosine_teta0);

private:
    Atoms<dim> *atoms;
};

