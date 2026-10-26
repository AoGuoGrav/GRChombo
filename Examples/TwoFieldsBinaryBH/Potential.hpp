/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef POTENTIAL_HPP_
#define POTENTIAL_HPP_

#include "simd.hpp"

//! The potential V(phi1, phi2) of the two scalar fields
/*!
    V = 1/2 mu1^2 phi1^2 + 1/2 mu2^2 phi2^2 + 1/2 lambda phi1^2 phi2^2

    To use a different potential, edit compute_potential below (and params_t
    if new parameters are needed, reading them in SimulationParameters.hpp).
    It must set the value V and both derivatives dV/dphi1 and dV/dphi2.
*/
class Potential
{
  public:
    struct params_t
    {
        double scalar_mass1; //!< mu1, the mass of phi1
        double scalar_mass2; //!< mu2, the mass of phi2
        double lambda;       //!< the coupling constant of the interaction
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
        const data_t phi1_sq = vars.phi1 * vars.phi1;
        const data_t phi2_sq = vars.phi2 * vars.phi2;

        // The potential value at (phi1, phi2)
        V_of_phi = 0.5 * mu1_sq * phi1_sq + 0.5 * mu2_sq * phi2_sq +
                   0.5 * m_params.lambda * phi1_sq * phi2_sq;

        // The potential gradients
        dVdphi1 = mu1_sq * vars.phi1 + m_params.lambda * vars.phi1 * phi2_sq;
        dVdphi2 = mu2_sq * vars.phi2 + m_params.lambda * phi1_sq * vars.phi2;
    }
};

#endif /* POTENTIAL_HPP_ */
