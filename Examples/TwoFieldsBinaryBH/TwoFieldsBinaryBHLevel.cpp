/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#include "TwoFieldsBinaryBHLevel.hpp"
#include "AMRReductions.hpp"
#include "BinaryBH.hpp"
#include "BoxLoops.hpp"
#include "ChiExtractionTaggingCriterion.hpp"
#include "ChiPunctureExtractionTaggingCriterion.hpp"
#include "ComputePack.hpp"
#include "NanCheck.hpp"
#include "PositiveChiAndAlpha.hpp"
#include "PunctureTracker.hpp"
#include "SetValue.hpp"
#include "SixthOrderDerivatives.hpp"
#include "SmallDataIO.hpp"
#include "TraceARemoval.hpp"
#include "TwoPuncturesInitialData.hpp"
#include "WeylExtraction.hpp"

// Matter includes
#include "EMTensor.hpp"
#include "InitialTwoFieldsData.hpp"
#include "MatterCCZ4RHS.hpp"
#include "MatterWeyl4.hpp"
#include "NewMatterConstraints.hpp"

// Things to do during the advance step after RK4 steps
void TwoFieldsBinaryBHLevel::specificAdvance()
{
    // Enforce the trace free A_ij condition and positive chi and alpha
    BoxLoops::loop(make_compute_pack(TraceARemoval(), PositiveChiAndAlpha()),
                   m_state_new, m_state_new, INCLUDE_GHOST_CELLS);

    // Check for nan's
    if (m_p.nan_check)
        BoxLoops::loop(
            NanCheck(m_dx, m_p.center, "NaNCheck in specific Advance"),
            m_state_new, m_state_new, EXCLUDE_GHOST_CELLS, disable_simd());
}

// The vacuum binary BH metric is set first, then the uniform scalar fields
// and the value of K which solves the Hamiltonian constraint
void TwoFieldsBinaryBHLevel::initialData()
{
    CH_TIME("TwoFieldsBinaryBHLevel::initialData");
    if (m_verbosity)
        pout() << "TwoFieldsBinaryBHLevel::initialData " << m_level << endl;
#ifdef USE_TWOPUNCTURES
    TwoPuncturesInitialData two_punctures_initial_data(
        m_dx, m_p.center, m_tp_amr.m_two_punctures);
    // Can't use simd with this initial data
    BoxLoops::loop(two_punctures_initial_data, m_state_new, m_state_new,
                   INCLUDE_GHOST_CELLS, disable_simd());
#else
    // Set up the compute class for the BinaryBH initial data
    BinaryBH binary(m_p.bh1_params, m_p.bh2_params, m_dx);

    // First set everything to zero (to avoid undefinded values in constraints)
    // then calculate initial data
    BoxLoops::loop(make_compute_pack(SetValue(0.), binary), m_state_new,
                   m_state_new, INCLUDE_GHOST_CELLS);
#endif

    // Now set the uniform fields and K (this overwrites K, which is zero in
    // the vacuum data above)
    Potential potential(m_p.potential_params);
    InitialTwoFieldsData initial_fields(m_p.initial_params, potential,
                                        m_p.G_Newton);
    BoxLoops::loop(initial_fields, m_state_new, m_state_new,
                   INCLUDE_GHOST_CELLS);
}

// Things to do when restarting from a checkpoint, including a restart from
// the output of an initial condition solver (e.g. GRTresna)
void TwoFieldsBinaryBHLevel::postRestart()
{
    // On restart calculate the constraints on every level
    fillAllGhosts();
    Potential potential(m_p.potential_params);
    TwoFieldsWithPotential two_fields(potential);
    BoxLoops::loop(
        MatterConstraints<TwoFieldsWithPotential>(
            two_fields, m_dx, m_p.G_Newton, c_Ham, Interval(c_Mom1, c_Mom3)),
        m_state_new, m_state_diagnostics, EXCLUDE_GHOST_CELLS);

    // restart works from level 0 to highest level, so want this to happen last
    // on finest level
    if (m_level == m_p.max_level)
    {
        AMRReductions<VariableType::diagnostic> amr_reductions(m_gr_amr);
        double L2_Ham = amr_reductions.norm(c_Ham);
        double L2_Mom = amr_reductions.norm(Interval(c_Mom1, c_Mom3));
        pout() << "The initial norm of the constraint vars on restart is "
               << L2_Ham << " for the Hamiltonian constraint and " << L2_Mom
               << " for the momentum constraints" << endl;
    }
}

// Calculate RHS during RK4 substeps
void TwoFieldsBinaryBHLevel::specificEvalRHS(GRLevelData &a_soln,
                                             GRLevelData &a_rhs,
                                             const double a_time)
{
    // Enforce positive chi and alpha and trace free A
    BoxLoops::loop(make_compute_pack(TraceARemoval(), PositiveChiAndAlpha()),
                   a_soln, a_soln, INCLUDE_GHOST_CELLS);

    // Calculate MatterCCZ4 right hand side with matter_t = TwoScalarFields
    Potential potential(m_p.potential_params);
    TwoFieldsWithPotential two_fields(potential);
    if (m_p.max_spatial_derivative_order == 4)
    {
        MatterCCZ4RHS<TwoFieldsWithPotential, MovingPunctureGauge,
                      FourthOrderDerivatives>
            my_ccz4_matter(two_fields, m_p.ccz4_params, m_dx, m_p.sigma,
                           m_p.formulation, m_p.G_Newton);
        BoxLoops::loop(my_ccz4_matter, a_soln, a_rhs, EXCLUDE_GHOST_CELLS);
    }
    else if (m_p.max_spatial_derivative_order == 6)
    {
        MatterCCZ4RHS<TwoFieldsWithPotential, MovingPunctureGauge,
                      SixthOrderDerivatives>
            my_ccz4_matter(two_fields, m_p.ccz4_params, m_dx, m_p.sigma,
                           m_p.formulation, m_p.G_Newton);
        BoxLoops::loop(my_ccz4_matter, a_soln, a_rhs, EXCLUDE_GHOST_CELLS);
    }
}

// enforce trace removal during RK4 substeps
void TwoFieldsBinaryBHLevel::specificUpdateODE(GRLevelData &a_soln,
                                               const GRLevelData &a_rhs,
                                               Real a_dt)
{
    // Enforce the trace free A_ij condition
    BoxLoops::loop(TraceARemoval(), a_soln, a_soln, INCLUDE_GHOST_CELLS);
}

void TwoFieldsBinaryBHLevel::preTagCells()
{
    // We only use chi in the tagging criterion so only fill the ghosts for chi
    fillAllGhosts(VariableType::evolution, Interval(c_chi, c_chi));
}

// specify the cells to tag
void TwoFieldsBinaryBHLevel::computeTaggingCriterion(
    FArrayBox &tagging_criterion, const FArrayBox &current_state,
    const FArrayBox &current_state_diagnostics)
{
    if (m_p.track_punctures)
    {
        std::vector<double> puncture_masses;
#ifdef USE_TWOPUNCTURES
        // use calculated bare masses from TwoPunctures
        puncture_masses = {m_tp_amr.m_two_punctures.mm,
                           m_tp_amr.m_two_punctures.mp};
#else
        puncture_masses = {m_p.bh1_params.mass, m_p.bh2_params.mass};
#endif /* USE_TWOPUNCTURES */
        auto puncture_coords =
            m_bh_amr.m_puncture_tracker.get_puncture_coords();
        BoxLoops::loop(ChiPunctureExtractionTaggingCriterion(
                           m_dx, m_level, m_p.max_level, m_p.extraction_params,
                           puncture_coords, m_p.activate_extraction,
                           m_p.track_punctures, puncture_masses),
                       current_state, tagging_criterion);
    }
    else
    {
        BoxLoops::loop(ChiExtractionTaggingCriterion(m_dx, m_level,
                                                     m_p.extraction_params,
                                                     m_p.activate_extraction),
                       current_state, tagging_criterion);
    }
}

void TwoFieldsBinaryBHLevel::specificPostTimeStep()
{
    CH_TIME("TwoFieldsBinaryBHLevel::specificPostTimeStep");

    bool first_step =
        (m_time == 0.); // this form is used when 'specificPostTimeStep' was
                        // called during setup at t=0 from Main
    // bool first_step = (m_time == m_dt); // if not called in Main

    Potential potential(m_p.potential_params);
    TwoFieldsWithPotential two_fields(potential);

    if (m_p.activate_extraction == 1)
    {
        int min_level = m_p.extraction_params.min_extraction_level();
        bool calculate_weyl = at_level_timestep_multiple(min_level);
        if (calculate_weyl)
        {
            // Populate the Weyl Scalar values on the grid
            fillAllGhosts();
            BoxLoops::loop(MatterWeyl4<TwoFieldsWithPotential>(
                               two_fields, m_p.extraction_params.center, m_dx,
                               m_p.formulation, m_p.G_Newton),
                           m_state_new, m_state_diagnostics,
                           EXCLUDE_GHOST_CELLS);

            // Do the extraction on the min extraction level
            if (m_level == min_level)
            {
                CH_TIME("WeylExtraction");
                // Now refresh the interpolator and do the interpolation
                // fill ghosts manually to minimise communication
                bool fill_ghosts = false;
                m_gr_amr.m_interpolator->refresh(fill_ghosts);
                m_gr_amr.fill_multilevel_ghosts(
                    VariableType::diagnostic, Interval(c_Weyl4_Re, c_Weyl4_Im),
                    min_level);
                WeylExtraction my_extraction(m_p.extraction_params, m_dt,
                                             m_time, first_step,
                                             m_restart_time);
                my_extraction.execute_query(m_gr_amr.m_interpolator);
            }
        }
    }

    if (m_p.calculate_constraint_norms)
    {
        fillAllGhosts();
        BoxLoops::loop(MatterConstraints<TwoFieldsWithPotential>(
                           two_fields, m_dx, m_p.G_Newton, c_Ham,
                           Interval(c_Mom1, c_Mom3)),
                       m_state_new, m_state_diagnostics, EXCLUDE_GHOST_CELLS);
        if (m_level == 0)
        {
            AMRReductions<VariableType::diagnostic> amr_reductions(m_gr_amr);
            double L2_Ham = amr_reductions.norm(c_Ham);
            double L2_Mom = amr_reductions.norm(Interval(c_Mom1, c_Mom3));
            SmallDataIO constraints_file(m_p.data_path + "constraint_norms",
                                         m_dt, m_time, m_restart_time,
                                         SmallDataIO::APPEND, first_step);
            constraints_file.remove_duplicate_time_data();
            if (first_step)
            {
                constraints_file.write_header_line({"L^2_Ham", "L^2_Mom"});
            }
            constraints_file.write_time_data_line({L2_Ham, L2_Mom});
        }
    }

    // do puncture tracking on requested level
    if (m_p.track_punctures && m_level == m_p.puncture_tracking_level)
    {
        CH_TIME("PunctureTracking");
        // only do the write out for every coarsest level timestep
        int coarsest_level = 0;
        bool write_punctures = at_level_timestep_multiple(coarsest_level);
        m_bh_amr.m_puncture_tracker.execute_tracking(m_time, m_restart_time,
                                                     m_dt, write_punctures);
    }

#ifdef USE_AHFINDER
    if (m_p.AH_activate && m_level == m_p.AH_params.level_to_run)
    {
        if (m_p.AH_set_origins_to_punctures && m_p.track_punctures)
        {
            m_bh_amr.m_ah_finder.set_origins(
                m_bh_amr.m_puncture_tracker.get_puncture_coords());
        }
        m_bh_amr.m_ah_finder.solve(m_dt, m_time, m_restart_time);
    }
#endif
}

#ifdef CH_USE_HDF5
// Things to do before a plot level - need to calculate the Weyl scalars,
// the constraints and the energy density
void TwoFieldsBinaryBHLevel::prePlotLevel()
{
    fillAllGhosts();
    Potential potential(m_p.potential_params);
    // EMTensor stores a reference to the matter object so it must stay in
    // scope for the whole loop
    TwoFieldsWithPotential two_fields(potential);
    MatterConstraints<TwoFieldsWithPotential> constraints(
        two_fields, m_dx, m_p.G_Newton, c_Ham, Interval(c_Mom1, c_Mom3));
    EMTensor<TwoFieldsWithPotential> emtensor(two_fields, m_dx, c_rho);
    if (m_p.activate_extraction == 1)
    {
        BoxLoops::loop(
            make_compute_pack(MatterWeyl4<TwoFieldsWithPotential>(
                                  two_fields, m_p.extraction_params.center,
                                  m_dx, m_p.formulation, m_p.G_Newton),
                              constraints, emtensor),
            m_state_new, m_state_diagnostics, EXCLUDE_GHOST_CELLS);
    }
    else
    {
        BoxLoops::loop(make_compute_pack(constraints, emtensor), m_state_new,
                       m_state_diagnostics, EXCLUDE_GHOST_CELLS);
    }
}
#endif /* CH_USE_HDF5 */
