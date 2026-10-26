/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#if !defined(TWOSCALARFIELDS_HPP_)
#error "This file should only be included through TwoScalarFields.hpp"
#endif

#ifndef TWOSCALARFIELDS_IMPL_HPP_
#define TWOSCALARFIELDS_IMPL_HPP_

// Calculate the stress energy tensor elements
template <class potential_t>
template <class data_t, template <typename> class vars_t>
emtensor_t<data_t> TwoScalarFields<potential_t>::compute_emtensor(
    const vars_t<data_t> &vars, const vars_t<Tensor<1, data_t>> &d1,
    const Tensor<2, data_t> &h_UU, const Tensor<3, data_t> &chris_ULL) const
{
    emtensor_t<data_t> out;

    // zero everything, then add in the kinetic part of each field
    out.rho = 0.0;
    FOR(i)
    {
        out.Si[i] = 0.0;
        FOR(j) { out.Sij[i][j] = 0.0; }
    }
    add_field_emtensor(out, vars, vars.Pi1, d1.phi1, h_UU);
    add_field_emtensor(out, vars, vars.Pi2, d1.phi2, h_UU);

    // compute potential and add constributions to EM Tensor (only once)
    data_t V_of_phi = 0.0;
    data_t dVdphi1 = 0.0;
    data_t dVdphi2 = 0.0;
    my_potential.compute_potential(V_of_phi, dVdphi1, dVdphi2, vars);

    out.rho += V_of_phi;
    FOR(i, j) { out.Sij[i][j] += -vars.h[i][j] * V_of_phi / vars.chi; }

    // S = Tr_S_ij of the total S_ij
    out.S = vars.chi * TensorAlgebra::compute_trace(out.Sij, h_UU);

    return out;
}

// Calculate the stress energy tensor elements excluding the potential
template <class potential_t>
template <class data_t, template <typename> class vars_t>
void TwoScalarFields<potential_t>::emtensor_excl_potential(
    emtensor_t<data_t> &out, const vars_t<data_t> &vars,
    const vars_t<Tensor<1, data_t>> &d1, const Tensor<2, data_t> &h_UU,
    const Tensor<3, data_t> &chris_ULL)
{
    out.rho = 0.0;
    FOR(i)
    {
        out.Si[i] = 0.0;
        FOR(j) { out.Sij[i][j] = 0.0; }
    }
    add_field_emtensor(out, vars, vars.Pi1, d1.phi1, h_UU);
    add_field_emtensor(out, vars, vars.Pi2, d1.phi2, h_UU);

    // S = Tr_S_ij
    out.S = vars.chi * TensorAlgebra::compute_trace(out.Sij, h_UU);
}

// Adds the kinetic part of the stress energy tensor of a single field
template <class potential_t>
template <class data_t, template <typename> class vars_t>
void TwoScalarFields<potential_t>::add_field_emtensor(
    emtensor_t<data_t> &out, const vars_t<data_t> &vars, const data_t &Pi,
    const Tensor<1, data_t> &d1_phi, const Tensor<2, data_t> &h_UU)
{
    // Useful quantity Vt
    data_t Vt = -Pi * Pi;
    FOR(i, j) { Vt += vars.chi * h_UU[i][j] * d1_phi[i] * d1_phi[j]; }

    // S_ij = T_ij
    FOR(i, j)
    {
        out.Sij[i][j] +=
            -0.5 * vars.h[i][j] * Vt / vars.chi + d1_phi[i] * d1_phi[j];
    }

    // S_i (note lower index) = - n^a T_ai
    FOR(i) { out.Si[i] += -d1_phi[i] * Pi; }

    // rho = n^a n^b T_ab
    out.rho += Pi * Pi + 0.5 * Vt;
}

// Adds in the RHS for the matter vars
template <class potential_t>
template <class data_t, template <typename> class vars_t,
          template <typename> class diff2_vars_t,
          template <typename> class rhs_vars_t>
void TwoScalarFields<potential_t>::add_matter_rhs(
    rhs_vars_t<data_t> &total_rhs, const vars_t<data_t> &vars,
    const vars_t<Tensor<1, data_t>> &d1,
    const diff2_vars_t<Tensor<2, data_t>> &d2,
    const vars_t<data_t> &advec) const
{
    // call the function for the rhs excluding the potential
    matter_rhs_excl_potential(total_rhs, vars, d1, d2, advec);

    // set the potential values
    data_t V_of_phi = 0.0;
    data_t dVdphi1 = 0.0;
    data_t dVdphi2 = 0.0;
    my_potential.compute_potential(V_of_phi, dVdphi1, dVdphi2, vars);

    // adjust RHS for the potential term
    total_rhs.Pi1 += -vars.lapse * dVdphi1;
    total_rhs.Pi2 += -vars.lapse * dVdphi2;
}

// the RHS excluding the potential terms
template <class potential_t>
template <class data_t, template <typename> class vars_t,
          template <typename> class diff2_vars_t,
          template <typename> class rhs_vars_t>
void TwoScalarFields<potential_t>::matter_rhs_excl_potential(
    rhs_vars_t<data_t> &rhs, const vars_t<data_t> &vars,
    const vars_t<Tensor<1, data_t>> &d1,
    const diff2_vars_t<Tensor<2, data_t>> &d2, const vars_t<data_t> &advec)
{
    using namespace TensorAlgebra;

    // the geometric quantities are shared by both fields
    const auto h_UU = compute_inverse_sym(vars.h);
    const auto chris = compute_christoffel(d1.h, h_UU);

    field_rhs_excl_potential(rhs.phi1, rhs.Pi1, vars.Pi1, d1.phi1, d2.phi1,
                             advec.phi1, advec.Pi1, vars, d1, h_UU, chris.ULL);
    field_rhs_excl_potential(rhs.phi2, rhs.Pi2, vars.Pi2, d1.phi2, d2.phi2,
                             advec.phi2, advec.Pi2, vars, d1, h_UU, chris.ULL);
}

// the RHS of a single field excluding the potential terms
template <class potential_t>
template <class data_t, template <typename> class vars_t>
void TwoScalarFields<potential_t>::field_rhs_excl_potential(
    data_t &rhs_phi, data_t &rhs_Pi, const data_t &Pi,
    const Tensor<1, data_t> &d1_phi, const Tensor<2, data_t> &d2_phi,
    const data_t &advec_phi, const data_t &advec_Pi, const vars_t<data_t> &vars,
    const vars_t<Tensor<1, data_t>> &d1, const Tensor<2, data_t> &h_UU,
    const Tensor<3, data_t> &chris_ULL)
{
    // evolution equations for scalar field and (minus) its conjugate momentum
    rhs_phi = vars.lapse * Pi + advec_phi;
    rhs_Pi = vars.lapse * vars.K * Pi + advec_Pi;

    FOR(i, j)
    {
        // includes non conformal parts of chris not included in chris_ULL
        rhs_Pi += h_UU[i][j] * (-0.5 * d1.chi[j] * vars.lapse * d1_phi[i] +
                                vars.chi * vars.lapse * d2_phi[i][j] +
                                vars.chi * d1.lapse[i] * d1_phi[j]);
        FOR(k)
        {
            rhs_Pi += -vars.chi * vars.lapse * h_UU[i][j] * chris_ULL[k][i][j] *
                      d1_phi[k];
        }
    }
}

#endif /* TWOSCALARFIELDS_IMPL_HPP_ */
