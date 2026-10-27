/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef FIELDENERGYDENSITIES_HPP_
#define FIELDENERGYDENSITIES_HPP_

#include "CCZ4Geometry.hpp"
#include "Cell.hpp"
#include "Coordinates.hpp"
#include "FourthOrderDerivatives.hpp"
#include "MatterCCZ4.hpp"
#include "Potential.hpp"
#include "Tensor.hpp"
#include "TwoScalarFields.hpp"
#include "UserVariables.hpp" //This files needs NUM_VARS - total no. components
#include "simd.hpp"

//! Calculates the energy density of each field and the integrands used to
//! measure how much of each field accumulates around the binary
/*!
    Stores in the diagnostic variables
    - rho1 = 1/2 Pi1^2 + 1/2 chi h^ij d_i phi1 d_j phi1 + 1/2 mu1^2 phi1^2 and
      the same for rho2: the energy density of each field excluding the
      interaction term V_int (which is only in the total rho);
    - rho1_dV, rho2_dV and rho_dV: sqrt(gamma) times rho1, rho2 and the total
      rho (including V_int), inside the sphere r < m_radius around m_center
      and zero elsewhere. Cells with chi < m_min_chi, i.e. inside the
      horizons, are also excluded since sqrt(gamma) = chi^(-3/2) diverges at
      the punctures.
    Summing rho_dV over the grid with AMRReductions gives
    int_{r < R, outside horizons} rho sqrt(gamma) d^3x.
*/
class FieldEnergyDensities
{
  public:
    typedef TwoScalarFields<Potential> matter_t;

    template <class data_t>
    using Vars = typename MatterCCZ4<matter_t>::template Vars<data_t>;

    //! Constructor
    FieldEnergyDensities(const matter_t &a_matter,
                         const Potential::params_t &a_potential_params,
                         const double a_dx,
                         const std::array<double, CH_SPACEDIM> a_center,
                         const double a_radius, const double a_min_chi)
        : m_matter(a_matter), m_potential_params(a_potential_params),
          m_dx(a_dx), m_deriv(a_dx), m_center(a_center), m_radius(a_radius),
          m_min_chi(a_min_chi)
    {
    }

    template <class data_t> void compute(Cell<data_t> current_cell) const
    {
        const auto vars = current_cell.template load_vars<Vars>();
        const auto d1 = m_deriv.template diff1<Vars>(current_cell);

        using namespace TensorAlgebra;
        const auto h_UU = compute_inverse_sym(vars.h);
        const auto chris = compute_christoffel(d1.h, h_UU);

        // total energy density, including the interaction
        const auto emtensor =
            m_matter.compute_emtensor(vars, d1, h_UU, chris.ULL);

        // energy densities of each field without the interaction
        data_t grad_phi1_sq = 0.;
        data_t grad_phi2_sq = 0.;
        FOR(i, j)
        {
            grad_phi1_sq += vars.chi * h_UU[i][j] * d1.phi1[i] * d1.phi1[j];
            grad_phi2_sq += vars.chi * h_UU[i][j] * d1.phi2[i] * d1.phi2[j];
        }
        const double mu1_sq =
            m_potential_params.scalar_mass1 * m_potential_params.scalar_mass1;
        const double mu2_sq =
            m_potential_params.scalar_mass2 * m_potential_params.scalar_mass2;
        const data_t rho1 = 0.5 * vars.Pi1 * vars.Pi1 + 0.5 * grad_phi1_sq +
                            0.5 * mu1_sq * vars.phi1 * vars.phi1;
        const data_t rho2 = 0.5 * vars.Pi2 * vars.Pi2 + 0.5 * grad_phi2_sq +
                            0.5 * mu2_sq * vars.phi2 * vars.phi2;

        // volume element masked to the integration region
        Coordinates<data_t> coords(current_cell, m_dx, m_center);
        const data_t r = coords.get_radius();
        data_t sqrt_gamma = pow(vars.chi, -1.5);
        sqrt_gamma =
            simd_conditional(simd_compare_lt(r, m_radius), sqrt_gamma, 0.0);
        sqrt_gamma = simd_conditional(simd_compare_gt(vars.chi, m_min_chi),
                                      sqrt_gamma, 0.0);

        current_cell.store_vars(rho1, c_rho1);
        current_cell.store_vars(rho2, c_rho2);
        current_cell.store_vars(sqrt_gamma * rho1, c_rho1_dV);
        current_cell.store_vars(sqrt_gamma * rho2, c_rho2_dV);
        current_cell.store_vars(sqrt_gamma * emtensor.rho, c_rho_dV);
    }

  protected:
    const matter_t m_matter;
    const Potential::params_t m_potential_params;
    const double m_dx;
    const FourthOrderDerivatives m_deriv;
    const std::array<double, CH_SPACEDIM> m_center;
    const double m_radius;  //!< radius of the integration sphere
    const double m_min_chi; //!< cells with chi below this are excluded
};

#endif /* FIELDENERGYDENSITIES_HPP_ */
