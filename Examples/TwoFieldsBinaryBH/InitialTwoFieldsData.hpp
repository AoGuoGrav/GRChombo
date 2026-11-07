/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef INITIALTWOFIELDSDATA_HPP_
#define INITIALTWOFIELDSDATA_HPP_

#include "Cell.hpp"
#include "Coordinates.hpp"
#include "Potential.hpp"
#include "Tensor.hpp"
#include "TwoScalarFields.hpp"
#include "UserVariables.hpp" //This files needs NUM_VARS - total no. components
#include "VarsTools.hpp"
#include "simd.hpp"
#include <random>
#include <vector>

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

    Optionally a small random seed is added to phi2 (phi2_noise_amplitude >
    0). With V_int = -lambda phi1 phi2^2 and phi2 = 0 the field phi2 stays
    exactly zero, so the decay phi1 -> phi2 phi2 (parametric resonance of
    phi2 driven by the oscillation of phi1) needs a seed to start from. The
    seed is a sum of num_modes waves
      delta phi2 = A (2 / sqrt(N)) sum_n cos(kx_n x + ky_n y + theta_n)
                   cos(kz_n z)
    with random directions, |k| uniform in [k_min, k_max] and random phases
    theta_n, so that its rms value is A, it is a smooth function of position
    (the same on every level and MPI rank, given the seed) and it is even in
    z (consistent with the reflective boundary at z = 0). It is not included
    in rho and K, which is fine as long as mu2^2 A^2 << rho.
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
        double phi2_noise_amplitude; //!< rms of the seed added to phi2
        double phi2_noise_k_min;     //!< smallest wavenumber of the seed
        double phi2_noise_k_max;     //!< largest wavenumber of the seed
        int phi2_noise_num_modes;    //!< number of waves in the seed
        int phi2_noise_seed;         //!< seed of the random number generator
    };

    template <class data_t>
    using MatterVars = TwoScalarFields<Potential>::Vars<data_t>;

    //! The constructor
    InitialTwoFieldsData(const params_t a_params, const Potential a_potential,
                         const double a_G_Newton, const double a_dx,
                         const std::array<double, CH_SPACEDIM> a_center)
        : m_params(a_params), m_potential(a_potential), m_G_Newton(a_G_Newton),
          m_dx(a_dx), m_center(a_center)
    {
        if (m_params.phi2_noise_amplitude > 0.0)
        {
            // the same sequence on every rank, since the seed is the same
            std::mt19937 generator(m_params.phi2_noise_seed);
            std::uniform_real_distribution<double> uniform(0.0, 1.0);
            for (int n = 0; n < m_params.phi2_noise_num_modes; ++n)
            {
                const double k =
                    m_params.phi2_noise_k_min +
                    (m_params.phi2_noise_k_max - m_params.phi2_noise_k_min) *
                        uniform(generator);
                const double cos_theta = 2.0 * uniform(generator) - 1.0;
                const double sin_theta = sqrt(1.0 - cos_theta * cos_theta);
                const double azimuth = 2.0 * M_PI * uniform(generator);
                m_kx.push_back(k * sin_theta * cos(azimuth));
                m_ky.push_back(k * sin_theta * sin(azimuth));
                m_kz.push_back(k * cos_theta);
                m_phase.push_back(2.0 * M_PI * uniform(generator));
            }
        }
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

        // add the seed of phi2
        const int num_modes = m_kx.size();
        if (num_modes > 0)
        {
            Coordinates<data_t> coords(current_cell, m_dx, m_center);
            data_t noise = 0.;
            for (int n = 0; n < num_modes; ++n)
            {
                noise +=
                    cos(m_kx[n] * coords.x + m_ky[n] * coords.y + m_phase[n]) *
                    cos(m_kz[n] * coords.z);
            }
            vars.phi2 += m_params.phi2_noise_amplitude * 2.0 /
                         sqrt(static_cast<double>(num_modes)) * noise;
        }

        current_cell.store_vars(vars);
        current_cell.store_vars(K, c_K);
    }

  protected:
    const params_t m_params; //!< The matter initial condition params
    const Potential m_potential;
    const double m_G_Newton;
    const double m_dx;
    const std::array<double, CH_SPACEDIM> m_center;
    //! wavevectors and phases of the waves in the seed of phi2
    std::vector<double> m_kx, m_ky, m_kz, m_phase;

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
