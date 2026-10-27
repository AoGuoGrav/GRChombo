/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef POTENTIAL_HPP_
#define POTENTIAL_HPP_

#include "simd.hpp"

//! The potential V(phi1, phi2) of the two scalar fields
/*!
    The interaction is the cubic coupling L_int = lambda psi^2 phi of
    Zhu, Piao & Zhang, arXiv:2509.10892, where psi is the field that appears
    squared. With interaction_squared_field = 1 (psi = phi1, phi = phi2)

      V = 1/2 mu1^2 phi1^2 + 1/2 mu2^2 phi2^2 - lambda phi1^2 phi2

    and with interaction_squared_field = 2 (psi = phi2, phi = phi1)

      V = 1/2 mu1^2 phi1^2 + 1/2 mu2^2 phi2^2 - lambda phi1 phi2^2

    so that, e.g. for the first case, Box phi1 - mu1^2 phi1 = -2 lambda phi1
    phi2 and Box phi2 - mu2^2 phi2 = -lambda phi1^2 as in their Eqs. (13-14).
    Note that a cubic potential is unbounded below for any lambda != 0.

    To use a different potential, edit compute_potential below (and params_t
    if new parameters are needed, reading them in SimulationParameters.hpp).
    It must set the value V and both derivatives dV/dphi1 and dV/dphi2.
*/
class Potential
{
  public:
    struct params_t
    {
        double scalar_mass1;           //!< mu1, the mass of phi1
        double scalar_mass2;           //!< mu2, the mass of phi2
        double lambda;                 //!< the cubic coupling constant
        int interaction_squared_field; //!< which field is squared (1 or 2)
    };

  private:
    params_t m_params;

  public:
    //! The constructor
    Potential(params_t a_params) : m_params(a_params) {}

    //! Set the potential function for the scalar fields here
    template <class data_t, template <typename> class vars_t>
    void compute_potential(data_t &V_of_phi, data_t &dVdphi1, data_t &dVdphi2,
                           const vars_t<data_t> &vars) const
    {
        const double mu1_sq = m_params.scalar_mass1 * m_params.scalar_mass1;
        const double mu2_sq = m_params.scalar_mass2 * m_params.scalar_mass2;

        // The free (mass) part
        V_of_phi = 0.5 * mu1_sq * vars.phi1 * vars.phi1 +
                   0.5 * mu2_sq * vars.phi2 * vars.phi2;
        dVdphi1 = mu1_sq * vars.phi1;
        dVdphi2 = mu2_sq * vars.phi2;

        // The cubic interaction V_int = -lambda psi^2 phi
        if (m_params.interaction_squared_field == 1)
        {
            V_of_phi += -m_params.lambda * vars.phi1 * vars.phi1 * vars.phi2;
            dVdphi1 += -2.0 * m_params.lambda * vars.phi1 * vars.phi2;
            dVdphi2 += -m_params.lambda * vars.phi1 * vars.phi1;
        }
        else
        {
            V_of_phi += -m_params.lambda * vars.phi1 * vars.phi2 * vars.phi2;
            dVdphi1 += -m_params.lambda * vars.phi2 * vars.phi2;
            dVdphi2 += -2.0 * m_params.lambda * vars.phi1 * vars.phi2;
        }
    }
};

#endif /* POTENTIAL_HPP_ */
