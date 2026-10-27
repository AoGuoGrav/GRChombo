/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef DIAGNOSTICVARIABLES_HPP
#define DIAGNOSTICVARIABLES_HPP

// assign an enum to each variable
enum
{
    c_Ham,

    c_Mom1,
    c_Mom2,
    c_Mom3,

    c_Weyl4_Re,
    c_Weyl4_Im,

    c_rho,

    // energy densities of each field excluding the interaction (see
    // FieldEnergyDensities.hpp)
    c_rho1,
    c_rho2,

    // sqrt(gamma) * rho1, rho2 and rho inside the integration sphere (zero
    // elsewhere), integrated to give data/rho_integrals.dat
    c_rho1_dV,
    c_rho2_dV,
    c_rho_dV,

    NUM_DIAGNOSTIC_VARS
};

namespace DiagnosticVariables
{
static const std::array<std::string, NUM_DIAGNOSTIC_VARS> variable_names = {
    "Ham",

    "Mom1",     "Mom2",     "Mom3",

    "Weyl4_Re", "Weyl4_Im",

    "rho",

    "rho1",     "rho2",

    "rho1_dV",  "rho2_dV",  "rho_dV"};
}

#endif /* DIAGNOSTICVARIABLES_HPP */
