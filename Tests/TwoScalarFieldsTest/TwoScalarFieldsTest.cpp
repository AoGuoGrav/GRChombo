/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

// Tests the TwoScalarFields matter class against the single field
// ScalarField class:
//  (a) with lambda = 0 the TwoScalarFields RHS of (phi1, Pi1) and (phi2, Pi2)
//      equals the ScalarField RHS run on each field alone
//  (b) the matter contribution to the metric RHS (MatterCCZ4RHS minus the
//      vacuum CCZ4RHS) and to the constraints (MatterConstraints minus
//      Constraints) equals the sum of the two single field contributions
//  (c) with the cubic coupling V_int = -lambda psi^2 phi of the
//      TwoFieldsBinaryBH example (psi = phi1 or phi2) the changes in rhs.Pi1,
//      rhs.Pi2 and Ham are -lapse dV_int/dphi1, -lapse dV_int/dphi2 and
//      -16 pi G V_int

// Chombo includes
#include "FArrayBox.H"

// Other includes
#ifdef _OPENMP
#include <omp.h>
#endif

#include "BoxLoops.hpp"
#include "CCZ4RHS.hpp"
#include "MatterCCZ4RHS.hpp"
#include "NewMatterConstraints.hpp"
// the two field potential tested is the one of the TwoFieldsBinaryBH example
#include "../../Examples/TwoFieldsBinaryBH/Potential.hpp"
#include "ScalarField.hpp"
#include "SingleFieldPotential.hpp"
#include "TwoScalarFields.hpp"
#include "UserVariables.hpp"
#include <iomanip>
#include <iostream>
#include <string>

// Chombo namespace
#include "UsingNamespace.H"

// relative tolerance for agreement to round-off
static const double TOLERANCE = 1e-12;

// largest relative error found in the checks of the current section
static double max_rel_err = 0.;

// max over the box of |a(comp_a) - b(comp_b)|
double max_abs_diff(const FArrayBox &a, int comp_a, const FArrayBox &b,
                    int comp_b, const Box &box)
{
    double max_diff = 0.;
    for (BoxIterator bit(box); bit.ok(); ++bit)
    {
        const IntVect iv = bit();
        max_diff = std::max(max_diff, std::abs(a(iv, comp_a) - b(iv, comp_b)));
    }
    return max_diff;
}

// check that the error is at round-off level relative to the scale of the
// quantities involved and print a message if not
bool check(const std::string &name, double err, double scale)
{
    const double rel_err = err / std::max(scale, 1.0);
    max_rel_err = std::max(max_rel_err, rel_err);
    if (rel_err > TOLERANCE)
    {
        std::cout << "  " << name << " DOES NOT AGREE: max error = " << err
                  << ", scale = " << scale << ", relative error = " << rel_err
                  << std::endl;
        return false;
    }
    return true;
}

void set_initial_data(FArrayBox &in_fab, const int N_GRID, const double dx)
{
    for (int zz = -3; zz < N_GRID + 3; ++zz)
    {
        const double z = zz * dx;
        for (int yy = -3; yy < N_GRID + 3; ++yy)
        {
            const double y = yy * dx;
            for (int xx = -3; xx < N_GRID + 3; ++xx)
            {
                const double x = xx * dx;
                const IntVect iv(xx, yy, zz);

                double g[3][3];
                double g_UU[3][3];
                double chi;
                {
                    g[0][0] = 1.36778 + 2.39731 * x + 4.53541 * x * x +
                              19.9771 * x * y * y * y + 6.13801 * y * z +
                              5.65185 * z * z + 9.35842 * z * z * z * z;
                    g[0][1] = -0.07646 - 0.48786 * x - 0.75098 * x * x -
                              1.73683 * x * y * y * y + 1.71676 * y * z +
                              1.03662 * z * z + 0.35630 * z * z * z * z;
                    g[0][2] = -0.10083 + 0.12445 * x - 1.26649 * x * x -
                              1.95052 * x * y * y * y + 0.73091 * y * z -
                              1.49835 * z * z - 2.39024 * z * z * z * z;
                    g[1][1] = 0.84072 + 2.31163 * x + 3.32275 * x * x +
                              15.1662 * x * y * y * y + 8.48730 * y * z +
                              3.05098 * z * z + 17.8448 * z * z * z * z;
                    g[1][2] = -0.42495 - 0.33464 * x - 0.47012 * x * x -
                              7.38477 * x * y * y * y + 0.41896 * y * z -
                              1.36394 * z * z + 5.25894 * z * z * z * z;
                    g[2][2] = 0.60995 + 1.30428 * x + 3.86237 * x * x +
                              22.7614 * x * y * y * y + 6.93818 * y * z +
                              4.39250 * z * z + 19.0244 * z * z * z * z;
                    g[1][0] = g[0][1];
                    g[2][0] = g[0][2];
                    g[2][1] = g[1][2];

                    const double detg = g[0][0] * g[1][1] * g[2][2] +
                                        2 * g[0][1] * g[0][2] * g[1][2] -
                                        g[0][0] * g[1][2] * g[1][2] -
                                        g[1][1] * g[0][2] * g[0][2] -
                                        g[2][2] * g[0][1] * g[0][1];
                    g_UU[0][0] = (g[1][1] * g[2][2] - g[1][2] * g[1][2]) / detg;
                    g_UU[0][1] = (g[0][2] * g[1][2] - g[0][1] * g[2][2]) / detg;
                    g_UU[0][2] = (g[0][1] * g[1][2] - g[0][2] * g[1][1]) / detg;
                    g_UU[1][1] = (g[0][0] * g[2][2] - g[0][2] * g[0][2]) / detg;
                    g_UU[1][2] = (g[0][1] * g[0][2] - g[0][0] * g[1][2]) / detg;
                    g_UU[2][2] = (g[0][0] * g[1][1] - g[0][1] * g[0][1]) / detg;
                    g_UU[1][0] = g_UU[0][1];
                    g_UU[2][0] = g_UU[0][2];
                    g_UU[2][1] = g_UU[1][2];

                    chi = pow(fabs(detg), -1.0 / GR_SPACEDIM);
                    in_fab(iv, c_chi) = chi;
                    in_fab(iv, c_h11) = chi * g[0][0];
                    in_fab(iv, c_h12) = chi * g[0][1];
                    in_fab(iv, c_h13) = chi * g[0][2];
                    in_fab(iv, c_h22) = chi * g[1][1];
                    in_fab(iv, c_h23) = chi * g[1][2];
                    in_fab(iv, c_h33) = chi * g[2][2];
                }

                {
                    double K[3][3];
                    K[0][0] = -0.16238 - 0.74295 * x + 0.51595 * x * x -
                              6.60239 * x * y * y * y - 0.76401 * y * z -
                              1.81131 * z * z - 3.88228 * z * z * z * z;
                    K[0][1] = 0.15054 - 0.60088 * x - 0.15428 * x * x +
                              3.16779 * x * y * y * y - 2.00687 * y * z -
                              1.35442 * z * z - 0.67601 * z * z * z * z;
                    K[0][2] = -0.02174 - 0.36243 * x + 0.81531 * x * x +
                              4.34918 * x * y * y * y + 0.90419 * y * z -
                              0.85088 * z * z - 6.45097 * z * z * z * z;
                    K[1][1] = -0.47653 - 0.43889 * x + 0.87342 * x * x +
                              4.24684 * x * y * y * y + 0.26290 * y * z +
                              1.90095 * z * z + 3.69515 * z * z * z * z;
                    K[1][2] = 0.37472 + 0.03657 * x - 0.10327 * x * x -
                              0.95744 * x * y * y * y - 1.20800 * y * z -
                              0.43064 * z * z - 0.25419 * z * z * z * z;
                    K[2][2] = 0.34184 + 0.21495 * x - 0.73195 * x * x +
                              7.81626 * x * y * y * y + 2.48359 * y * z +
                              1.89657 * z * z - 4.10980 * z * z * z * z;
                    K[1][0] = K[0][1];
                    K[2][0] = K[0][2];
                    K[2][1] = K[1][2];

                    double trK = 0;
                    FOR(a, b) { trK += g_UU[a][b] * K[a][b]; }

                    in_fab(iv, c_K) = trK;
                    in_fab(iv, c_A11) =
                        chi * (K[0][0] - trK * g[0][0] / GR_SPACEDIM);
                    in_fab(iv, c_A12) =
                        chi * (K[0][1] - trK * g[0][1] / GR_SPACEDIM);
                    in_fab(iv, c_A13) =
                        chi * (K[0][2] - trK * g[0][2] / GR_SPACEDIM);
                    in_fab(iv, c_A22) =
                        chi * (K[1][1] - trK * g[1][1] / GR_SPACEDIM);
                    in_fab(iv, c_A23) =
                        chi * (K[1][2] - trK * g[1][2] / GR_SPACEDIM);
                    in_fab(iv, c_A33) =
                        chi * (K[2][2] - trK * g[2][2] / GR_SPACEDIM);
                }

                in_fab(iv, c_Theta) = 0.12345 + 0.04321 * x - 0.30402 * y * y +
                                      0.21345 * x * y * z;

                in_fab(iv, c_Gamma1) =
                    -0.49482 + 0.89227 * x + 0.05571 * x * x -
                    5.38570 * x * y * y * y + 0.13979 * y * z -
                    0.68588 * z * z - 4.39964 * z * z * z * z;
                in_fab(iv, c_Gamma2) =
                    -0.09082 - 0.31017 * x + 1.06980 * x * x +
                    7.81524 * x * y * y * y - 1.65016 * y * z -
                    0.53352 * z * z - 3.20997 * z * z * z * z;
                in_fab(iv, c_Gamma3) =
                    -0.42367 + 0.03891 * x - 0.87898 * x * x +
                    6.67657 * x * y * y * y - 3.44662 * y * z -
                    0.19655 * z * z + 2.97524 * z * z * z * z;

                in_fab(iv, c_lapse) = 0.73578 + 0.36898 * x + 0.64348 * x * x +
                                      9.33487 * x * y * y * y +
                                      0.99469 * y * z + 0.20515 * z * z +
                                      8.88385 * z * z * z * z;
                in_fab(iv, c_shift1) = 0.00000 + 0.18795 * x - 0.52389 * x * x -
                                       4.14079 * x * y * y * y +
                                       0.73135 * y * z - 0.27057 * z * z +
                                       3.24187 * z * z * z * z;
                in_fab(iv, c_shift2) = 0.00000 - 0.30316 * x - 0.15184 * x * x -
                                       0.48815 * x * y * y * y +
                                       2.45991 * y * z - 0.79248 * z * z +
                                       7.14007 * z * z * z * z;
                in_fab(iv, c_shift3) = 0.00000 + 0.68835 * x - 0.52219 * x * x -
                                       7.50449 * x * y * y * y -
                                       2.35372 * y * z - 0.21476 * z * z +
                                       4.36363 * z * z * z * z;
                in_fab(iv, c_B1) = -0.26928 + 0.35045 * x - 0.48884 * x * x +
                                   2.72465 * x * y * y * y - 2.59022 * y * z -
                                   0.27384 * z * z + 0.38748 * z * z * z * z;
                in_fab(iv, c_B2) = 0.40234 + 0.26741 * x + 1.94822 * x * x -
                                   0.78276 * x * y * y * y + 2.12346 * y * z +
                                   0.69086 * z * z - 4.47639 * z * z * z * z;
                in_fab(iv, c_B3) = 0.40313 + 0.00569 * x - 1.12452 * x * x -
                                   5.49255 * x * y * y * y - 2.21932 * y * z +
                                   0.49523 * z * z + 1.29460 * z * z * z * z;

                in_fab(iv, c_phi1) = 0.34578 + 0.26898 * x + 0.54348 * x * x +
                                     0.33487 * x * y * y * y + 0.79469 * y * z +
                                     0.30515 * z * z + 1.88385 * z * z * z * z;
                in_fab(iv, c_Pi1) = 0.65668 + 0.20188 * x + 0.34348 * x * x +
                                    0.31787 * x * y * y * y + 0.88469 * y * z +
                                    0.10515 * z * z + 1.88385 * z * z * z * z;
                in_fab(iv, c_phi2) = -0.21734 + 0.41257 * y - 0.62019 * y * y +
                                     0.51983 * x * x * z - 0.37741 * x * z +
                                     0.92813 * y * y * z * z -
                                     1.20155 * x * x * x;
                in_fab(iv, c_Pi2) = 0.11972 - 0.57390 * z + 0.29871 * z * z +
                                    0.83651 * x * y - 0.44019 * x * x * y +
                                    0.61732 * y * y * y - 0.98123 * x * z * z;
            }
        }
    }
}

int main()
{
#ifdef _OPENMP
    std::cout << "#threads = " << omp_get_max_threads() << std::endl;
#endif

    const int N_GRID = 16;
    Box box(IntVect(0, 0, 0), IntVect(N_GRID - 1, N_GRID - 1, N_GRID - 1));
    Box ghosted_box(IntVect(-3, -3, -3),
                    IntVect(N_GRID + 2, N_GRID + 2, N_GRID + 2));

    const double dx = 0.5 / N_GRID;

    FArrayBox in_fab(ghosted_box, NUM_VARS);
    in_fab.setVal(0.);
    set_initial_data(in_fab, N_GRID, dx);

    // copy of the data with field 2 stored in the slots of field 1 so that the
    // single field ScalarField class (which uses c_phi = c_phi1, c_Pi = c_Pi1)
    // can be run on field 2
    FArrayBox in_fab_field2(ghosted_box, NUM_VARS);
    in_fab_field2.copy(in_fab);
    in_fab_field2.copy(in_fab, c_phi2, c_phi1, 1);
    in_fab_field2.copy(in_fab, c_Pi2, c_Pi1, 1);

    CCZ4_params_t<MovingPunctureGauge::params_t> params;
    params.kappa1 = 0.1;
    params.kappa2 = -0.3;
    params.kappa3 = 0.5;
    params.covariantZ4 = true;
    params.shift_Gamma_coeff = 0.75;
    params.lapse_advec_coeff = 1.0;
    params.lapse_power = 1.0;
    params.lapse_coeff = 2.0;
    params.shift_advec_coeff = 1.0;
    params.eta = 1.0;

    const double scalar_mass1 = 1.1;
    const double scalar_mass2 = 0.7;
    const double lambda = 2.3;
    const double G_Newton = 0.8;
    const double sigma = 0.1;

    typedef ScalarField<SingleFieldPotential> SingleField;
    typedef TwoScalarFields<Potential> TwoFields;

    SingleField field1(SingleFieldPotential{scalar_mass1});
    SingleField field2(SingleFieldPotential{scalar_mass2});
    TwoFields two_fields_free(Potential({scalar_mass1, scalar_mass2, 0., 1}));

    bool passed = true;

    for (int formulation : {CCZ4RHS<>::USE_CCZ4, CCZ4RHS<>::USE_BSSN})
    {
        std::cout << "Formulation: "
                  << (formulation == CCZ4RHS<>::USE_CCZ4 ? "CCZ4" : "BSSN")
                  << std::endl;

        // the kappas must vanish for BSSN
        if (formulation == CCZ4RHS<>::USE_BSSN)
        {
            params.kappa1 = 0.;
            params.kappa2 = 0.;
            params.kappa3 = 0.;
        }

        FArrayBox out_two(box, NUM_VARS);
        FArrayBox out_one1(box, NUM_VARS);
        FArrayBox out_one2(box, NUM_VARS);
        FArrayBox out_vac(box, NUM_VARS);
        for (FArrayBox *fab : {&out_two, &out_one1, &out_one2, &out_vac})
            fab->setVal(0.);

        // RHS
        BoxLoops::loop(
            MatterCCZ4RHS<TwoFields, MovingPunctureGauge,
                          FourthOrderDerivatives>(two_fields_free, params, dx,
                                                  sigma, formulation, G_Newton),
            in_fab, out_two);
        BoxLoops::loop(MatterCCZ4RHS<SingleField, MovingPunctureGauge,
                                     FourthOrderDerivatives>(
                           field1, params, dx, sigma, formulation, G_Newton),
                       in_fab, out_one1);
        BoxLoops::loop(MatterCCZ4RHS<SingleField, MovingPunctureGauge,
                                     FourthOrderDerivatives>(
                           field2, params, dx, sigma, formulation, G_Newton),
                       in_fab_field2, out_one2);
        BoxLoops::loop(CCZ4RHS<MovingPunctureGauge, FourthOrderDerivatives>(
                           params, dx, sigma, formulation),
                       in_fab, out_vac);

        // Constraints
        BoxLoops::loop(MatterConstraints<TwoFields>(two_fields_free, dx,
                                                    G_Newton, c_Ham,
                                                    Interval(c_Mom1, c_Mom3)),
                       in_fab, out_two);
        BoxLoops::loop(MatterConstraints<SingleField>(field1, dx, G_Newton,
                                                      c_Ham,
                                                      Interval(c_Mom1, c_Mom3)),
                       in_fab, out_one1);
        BoxLoops::loop(MatterConstraints<SingleField>(field2, dx, G_Newton,
                                                      c_Ham,
                                                      Interval(c_Mom1, c_Mom3)),
                       in_fab_field2, out_one2);
        BoxLoops::loop(Constraints(dx, c_Ham, Interval(c_Mom1, c_Mom3)), in_fab,
                       out_vac);

        // (a) field RHS equal to the single field ones when lambda = 0
        bool passed_a = true;
        max_rel_err = 0.;
        passed_a &= check("rhs phi1",
                          max_abs_diff(out_two, c_phi1, out_one1, c_phi, box),
                          out_one1.norm(0, c_phi, 1));
        passed_a &=
            check("rhs Pi1", max_abs_diff(out_two, c_Pi1, out_one1, c_Pi, box),
                  out_one1.norm(0, c_Pi, 1));
        passed_a &= check("rhs phi2",
                          max_abs_diff(out_two, c_phi2, out_one2, c_phi, box),
                          out_one2.norm(0, c_phi, 1));
        passed_a &=
            check("rhs Pi2", max_abs_diff(out_two, c_Pi2, out_one2, c_Pi, box),
                  out_one2.norm(0, c_Pi, 1));
        std::cout << "  (a) single field RHS comparison "
                  << (passed_a ? "passed" : "FAILED") << " (max relative error "
                  << max_rel_err << ")" << std::endl;

        // (b) matter contributions to the metric RHS and constraints add up
        // (two - vac) - (one1 - vac) - (one2 - vac) = two - one1 - one2 + vac
        bool passed_b = true;
        max_rel_err = 0.;
        std::vector<int> metric_and_constraint_comps;
        for (int icomp = 0; icomp < NUM_CCZ4_VARS; ++icomp)
            metric_and_constraint_comps.push_back(icomp);
        for (int icomp = c_Ham; icomp <= c_Mom3; ++icomp)
            metric_and_constraint_comps.push_back(icomp);
        for (int icomp : metric_and_constraint_comps)
        {
            FArrayBox residual(box, 1);
            residual.copy(out_two, icomp, 0, 1);
            residual.minus(out_one1, icomp, 0, 1);
            residual.minus(out_one2, icomp, 0, 1);
            residual.plus(out_vac, icomp, 0, 1);
            const double scale = std::max(
                {out_two.norm(0, icomp, 1), out_one1.norm(0, icomp, 1),
                 out_one2.norm(0, icomp, 1), out_vac.norm(0, icomp, 1)});
            passed_b &= check("metric/constraint component " +
                                  UserVariables::variable_names[icomp],
                              residual.norm(0, 0, 1), scale);
        }
        std::cout << "  (b) additivity of matter terms "
                  << (passed_b ? "passed" : "FAILED") << " (max relative error "
                  << max_rel_err << ")" << std::endl;

        // (c) cubic interaction terms in the Pi RHS and in Ham, for both
        // choices of the squared field
        bool passed_c = true;
        max_rel_err = 0.;
        for (int squared_field : {1, 2})
        {
            TwoFields two_fields_coupled(
                Potential({scalar_mass1, scalar_mass2, lambda, squared_field}));
            FArrayBox out_two_coupled(box, NUM_VARS);
            out_two_coupled.setVal(0.);
            BoxLoops::loop(MatterCCZ4RHS<TwoFields, MovingPunctureGauge,
                                         FourthOrderDerivatives>(
                               two_fields_coupled, params, dx, sigma,
                               formulation, G_Newton),
                           in_fab, out_two_coupled);
            BoxLoops::loop(
                MatterConstraints<TwoFields>(two_fields_coupled, dx, G_Newton,
                                             c_Ham, Interval(c_Mom1, c_Mom3)),
                in_fab, out_two_coupled);

            double err_Pi1 = 0., err_Pi2 = 0., err_Ham = 0.;
            double scale_Pi1 = 0., scale_Pi2 = 0., scale_Ham = 0.;
            for (BoxIterator bit(box); bit.ok(); ++bit)
            {
                const IntVect iv = bit();
                const double lapse = in_fab(iv, c_lapse);
                const double phi1 = in_fab(iv, c_phi1);
                const double phi2 = in_fab(iv, c_phi2);
                // psi is the squared field, phi the other one
                const double psi = (squared_field == 1) ? phi1 : phi2;
                const double phi = (squared_field == 1) ? phi2 : phi1;
                const double V_int = -lambda * psi * psi * phi;
                const double dVint_dpsi = -2.0 * lambda * psi * phi;
                const double dVint_dphi = -lambda * psi * psi;
                const double dVint_dphi1 =
                    (squared_field == 1) ? dVint_dpsi : dVint_dphi;
                const double dVint_dphi2 =
                    (squared_field == 1) ? dVint_dphi : dVint_dpsi;

                const double dPi1 =
                    out_two_coupled(iv, c_Pi1) - out_two(iv, c_Pi1);
                const double dPi2 =
                    out_two_coupled(iv, c_Pi2) - out_two(iv, c_Pi2);
                const double dHam =
                    out_two_coupled(iv, c_Ham) - out_two(iv, c_Ham);
                err_Pi1 =
                    std::max(err_Pi1, std::abs(dPi1 + lapse * dVint_dphi1));
                err_Pi2 =
                    std::max(err_Pi2, std::abs(dPi2 + lapse * dVint_dphi2));
                err_Ham = std::max(
                    err_Ham, std::abs(dHam + 16. * M_PI * G_Newton * V_int));
                scale_Pi1 = std::max(scale_Pi1, std::abs(out_two(iv, c_Pi1)));
                scale_Pi2 = std::max(scale_Pi2, std::abs(out_two(iv, c_Pi2)));
                scale_Ham = std::max(scale_Ham, std::abs(out_two(iv, c_Ham)));
            }
            const std::string label =
                " (squared field " + std::to_string(squared_field) + ")";
            passed_c &= check("interaction term in rhs Pi1" + label, err_Pi1,
                              scale_Pi1);
            passed_c &= check("interaction term in rhs Pi2" + label, err_Pi2,
                              scale_Pi2);
            passed_c &=
                check("interaction term in Ham" + label, err_Ham, scale_Ham);
            // the field RHS and Mom must not depend on lambda
            for (int icomp : {c_phi1, c_phi2, c_Mom1, c_Mom2, c_Mom3})
                passed_c &= check(
                    UserVariables::variable_names[icomp] +
                        " independent of lambda" + label,
                    max_abs_diff(out_two_coupled, icomp, out_two, icomp, box),
                    out_two.norm(0, icomp, 1));
        }
        std::cout << "  (c) cubic interaction terms "
                  << (passed_c ? "passed" : "FAILED") << " (max relative error "
                  << max_rel_err << ")" << std::endl;

        passed &= passed_a && passed_b && passed_c;
    }

    if (passed)
        std::cout << "TwoScalarFields test passed..." << std::endl;
    else
        std::cout << "TwoScalarFields test failed..." << std::endl;

    return passed ? 0 : -1;
}
