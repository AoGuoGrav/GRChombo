/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef TWOFIELDSBINARYBHLEVEL_HPP_
#define TWOFIELDSBINARYBHLEVEL_HPP_

#include "DefaultLevelFactory.hpp"
#include "GRAMRLevel.hpp"
// TPAMR.hpp includes BHAMR.hpp
#include "TPAMR.hpp"
// Problem specific includes
#include "Potential.hpp"
#include "TwoScalarFields.hpp"

//!  A class for the evolution of two interacting scalar fields around a
//!  binary black hole, minimally coupled to gravity
/*!
     The binary BH metric is set either with the approximate BinaryBH initial
     data or with TwoPunctures, then the two fields are set to uniform values
     and K is adjusted so that the constraints are satisfied (see
     InitialTwoFieldsData). Everything is then evolved with the CCZ4 equations
     including the matter terms. \sa TwoScalarFields(), MatterCCZ4RHS(),
     MatterConstraints()
*/
class TwoFieldsBinaryBHLevel : public GRAMRLevel
{
    friend class DefaultLevelFactory<TwoFieldsBinaryBHLevel>;
    // Inherit the contructors from GRAMRLevel
    using GRAMRLevel::GRAMRLevel;

    BHAMR &m_bh_amr = dynamic_cast<BHAMR &>(m_gr_amr);
#ifdef USE_TWOPUNCTURES
    TPAMR &m_tp_amr = dynamic_cast<TPAMR &>(m_gr_amr);
#endif /* USE_TWOPUNCTURES */

    // Typedef for the scalar fields
    typedef TwoScalarFields<Potential> TwoFieldsWithPotential;

    /// Things to do at every full timestep
    ///(might include several substeps, e.g. in RK4)
    virtual void specificAdvance() override;

    /// Initial data calculation
    virtual void initialData() override;

    /// Things to do when restarting from a checkpoint (or from the output of
    /// an initial condition solver)
    virtual void postRestart() override;

    /// Calculation of the right hand side for the time stepping
    virtual void specificEvalRHS(GRLevelData &a_soln, GRLevelData &a_rhs,
                                 const double a_time) override;

    /// Things to do after dt*rhs has been added to the solution
    virtual void specificUpdateODE(GRLevelData &a_soln,
                                   const GRLevelData &a_rhs,
                                   Real a_dt) override;

    /// Things to do before tagging cells (i.e. filling ghosts)
    virtual void preTagCells() override;

    /// Identify and tag the cells that need higher resolution
    virtual void computeTaggingCriterion(
        FArrayBox &tagging_criterion, const FArrayBox &current_state,
        const FArrayBox &current_state_diagnostics) override;

    // to do post each time step on every level
    virtual void specificPostTimeStep() override;

#ifdef CH_USE_HDF5
    /// Any actions that should happen just before plot files output
    virtual void prePlotLevel() override;
#endif /* CH_USE_HDF5 */
};

#endif /* TWOFIELDSBINARYBHLEVEL_HPP_ */
