/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef SINGLEFIELDPOTENTIAL_HPP_
#define SINGLEFIELDPOTENTIAL_HPP_

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

#endif /* SINGLEFIELDPOTENTIAL_HPP_ */
