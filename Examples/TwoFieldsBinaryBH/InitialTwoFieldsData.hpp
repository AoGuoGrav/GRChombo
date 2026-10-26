/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef INITIALTWOFIELDSDATA_HPP_
#define INITIALTWOFIELDSDATA_HPP_

#include "Cell.hpp"
#include "Potential.hpp"
#include "Tensor.hpp"
#include "TwoScalarFields.hpp"
#include "UserVariables.hpp" //This files needs NUM_VARS - total no. components
#include "VarsTools.hpp"
#include "simd.hpp"

//! Class which sets uniform initial values of the two scalar fields and the
//! value of K required by the Hamiltonian constraint
/*!
    This must be applied after the vacuum binary BH metric (BinaryBH or
    TwoPunctures) has been set, since it overwrites K.

    Why this is consistent: both BinaryBH and TwoPunctures (Bowen-York) data
    are conformally flat with K = 0 and a tracefree A_ij. The fields are
    uniform, so d_i phi_a = 0, which gives S_i = 0 and a constant energy
    density rho = 1/2 Pi1^2 + 1/2 Pi2^2 + V(phi1_0, phi2_0). Setting K to a
    constant leaves the momentum constraint
      D_j A^ij - (2/3) D^i K = 8 pi G S^i
    equal to the vacuum one, since D^i K = 0 = S^i. The Hamiltonian constraint
      R + (2/3) K^2 - A_ij A^ij = 16 pi G rho
    reduces to the vacuum one, R - A_ij A^ij = 0, because the choice
      K = K_sign * sqrt(24 pi G rho)
    makes (2/3) K^2 cancel 16 pi G rho exactly. Hence TwoPunctures plus
    uniform fields solves the constraints (to the accuracy of the
    TwoPunctures solve), and the approximate BinaryBH data has the same
    constraint violation as the vacuum BinaryBH example.

    K_sign = -1 (the default) corresponds to an expanding background, as in
    GRChombo's convention K < 0 for expansion (see e.g.
    Examples/ScalarFieldCosmo/InitialK.hpp). Note that Examples/ScalarFieldBH
    uses K > 0.
*/
class InitialTwoFieldsData
{
  public:
    //! A structure for the input params for the initial field values
    struct params_t
    {
        double phi1_0; //!< Initial uniform value of phi1
        double phi2_0; //!< Initial uniform value of phi2
        double Pi1_0;  //!< Initial uniform value of Pi1
        double Pi2_0;  //!< Initial uniform value of Pi2
        double K_sign; //!< Sign of K: -1 for expansion, +1 for contraction
    };

    template <class data_t>
    using MatterVars = TwoScalarFields<Potential>::Vars<data_t>;

    //! The constructor
    InitialTwoFieldsData(const params_t a_params, const Potential a_potential,
                         const double a_G_Newton = 1.0)
        : m_params(a_params), m_potential(a_potential), m_G_Newton(a_G_Newton)
    {
    }

    //! The (uniform) energy density of the initial field configuration
    static double compute_rho(const params_t &a_params,
                              const Potential &a_potential)
    {
        MatterVars<double> vars;
        set_matter_vars(vars, a_params);
        double V_of_phi, dVdphi1, dVdphi2;
        a_potential.compute_potential(V_of_phi, dVdphi1, dVdphi2, vars);
        return 0.5 * vars.Pi1 * vars.Pi1 + 0.5 * vars.Pi2 * vars.Pi2 + V_of_phi;
    }

    //! The value of K which solves the Hamiltonian constraint
    static double compute_K(const params_t &a_params,
                            const Potential &a_potential,
                            const double a_G_Newton)
    {
        const double rho = compute_rho(a_params, a_potential);
        return a_params.K_sign * sqrt(24.0 * M_PI * a_G_Newton * rho);
    }

    //! Function to compute the value of the field vars and K on the grid
    template <class data_t> void compute(Cell<data_t> current_cell) const
    {
        MatterVars<data_t> vars;
        set_matter_vars(vars, m_params);
        data_t K = compute_K(m_params, m_potential, m_G_Newton);

        current_cell.store_vars(vars);
        current_cell.store_vars(K, c_K);
    }

  protected:
    const params_t m_params; //!< The matter initial condition params
    const Potential m_potential;
    const double m_G_Newton;

    template <class data_t>
    static void set_matter_vars(MatterVars<data_t> &vars,
                                const params_t &a_params)
    {
        vars.phi1 = a_params.phi1_0;
        vars.Pi1 = a_params.Pi1_0;
        vars.phi2 = a_params.phi2_0;
        vars.Pi2 = a_params.Pi2_0;
    }
};

#endif /* INITIALTWOFIELDSDATA_HPP_ */
