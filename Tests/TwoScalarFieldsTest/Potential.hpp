/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef POTENTIAL_HPP_
#define POTENTIAL_HPP_

#include "Tensor.hpp"
#include "simd.hpp"

//! V = 1/2 m^2 phi^2 for the single field ScalarField reference class
class SingleFieldPotential
{
  public:
    SingleFieldPotential(double a_scalar_mass) : m_scalar_mass(a_scalar_mass) {}

    template <class data_t, template <typename> class vars_t>
    void compute_potential(data_t &V_of_phi, data_t &dVdphi,
                           const vars_t<data_t> &vars) const
    {
        V_of_phi = 0.5 * pow(m_scalar_mass * vars.phi, 2.0);
        dVdphi = pow(m_scalar_mass, 2.0) * vars.phi;
    }

  private:
    double m_scalar_mass;
};

//! V = 1/2 m1^2 phi1^2 + 1/2 m2^2 phi2^2 + 1/2 lambda phi1^2 phi2^2
class TwoFieldPotential
{
  public:
    struct params_t
    {
        double scalar_mass1;
        double scalar_mass2;
        double lambda;
    };

    TwoFieldPotential(params_t a_params) : m_params(a_params) {}

    template <class data_t, template <typename> class vars_t>
    void compute_potential(data_t &V_of_phi, data_t &dVdphi1, data_t &dVdphi2,
                           const vars_t<data_t> &vars) const
    {
        const double m1sq = m_params.scalar_mass1 * m_params.scalar_mass1;
        const double m2sq = m_params.scalar_mass2 * m_params.scalar_mass2;
        const data_t phi1sq = vars.phi1 * vars.phi1;
        const data_t phi2sq = vars.phi2 * vars.phi2;

        V_of_phi = 0.5 * m1sq * phi1sq + 0.5 * m2sq * phi2sq +
                   0.5 * m_params.lambda * phi1sq * phi2sq;
        dVdphi1 = m1sq * vars.phi1 + m_params.lambda * vars.phi1 * phi2sq;
        dVdphi2 = m2sq * vars.phi2 + m_params.lambda * phi1sq * vars.phi2;
    }

  private:
    params_t m_params;
};

#endif /* POTENTIAL_HPP_ */
