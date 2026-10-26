/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef TWOSCALARFIELDS_HPP_
#define TWOSCALARFIELDS_HPP_

#include "CCZ4Geometry.hpp"
#include "DimensionDefinitions.hpp"
#include "FourthOrderDerivatives.hpp"
#include "Tensor.hpp"
#include "TensorAlgebra.hpp"
#include "UserVariables.hpp" //This files needs NUM_VARS, total num of components
#include "VarsTools.hpp"

//!  Calculates the matter type specific elements for two interacting real
//!  scalar fields, minimally coupled to gravity
/*!
     This class is a matter_t object (see ScalarField for the single field
     version, with which it shares its interface) for two real scalar fields
     phi1 and phi2 described by the action

       S = \int d^4x \sqrt{-g} [ - 1/2 g^{ab} d_a phi1 d_b phi1
                                 - 1/2 g^{ab} d_a phi2 d_b phi2
                                 - V(phi1, phi2) ]

     with signature (-,+,+,+), so that each field obeys
     Box phi_a = dV/dphi_a. The fields interact only through V(phi1, phi2).

     The evolved momenta follow the ScalarField convention:
       Pi_a = n^mu d_mu phi_a = (d_t phi_a - beta^i d_i phi_a) / lapse,
     i.e. d_t phi_a = lapse * Pi_a + beta^i d_i phi_a. This is (minus) the
     conjugate momentum in the sign convention where Pi = -n^mu d_mu phi.

     The class is templated over a potential function potential_t which the
     user must specify in a class with the member function
       compute_potential(V, dVdphi1, dVdphi2, vars)
     which sets the value of the potential and its two derivatives.
     \sa ScalarField(), MatterCCZ4RHS(), MatterConstraints()
*/
template <class potential_t> class TwoScalarFields
{
  protected:
    //! The local copy of the potential
    potential_t my_potential;

  public:
    //!  Constructor of class TwoScalarFields, inputs are the matter parameters.
    TwoScalarFields(const potential_t a_potential) : my_potential(a_potential)
    {
    }

    //! Structure containing the rhs variables for the matter fields
    template <class data_t> struct Vars
    {
        data_t phi1;
        data_t Pi1;
        data_t phi2;
        data_t Pi2;

        /// Defines the mapping between members of Vars and Chombo grid
        /// variables (enum in User_Variables)
        template <typename mapping_function_t>
        void enum_mapping(mapping_function_t mapping_function)
        {
            VarsTools::define_enum_mapping(mapping_function, c_phi1, phi1);
            VarsTools::define_enum_mapping(mapping_function, c_Pi1, Pi1);
            VarsTools::define_enum_mapping(mapping_function, c_phi2, phi2);
            VarsTools::define_enum_mapping(mapping_function, c_Pi2, Pi2);
        }
    };

    //! Structure containing the rhs variables for the matter fields requiring
    //!  2nd derivs
    template <class data_t> struct Diff2Vars
    {
        data_t phi1;
        data_t phi2;

        /// Defines the mapping between members of Vars and Chombo grid
        ///  variables (enum in User_Variables)
        template <typename mapping_function_t>
        void enum_mapping(mapping_function_t mapping_function)
        {
            VarsTools::define_enum_mapping(mapping_function, c_phi1, phi1);
            VarsTools::define_enum_mapping(mapping_function, c_phi2, phi2);
        }
    };

    //! The function which calculates the EM Tensor, given the vars and
    //! derivatives, including the potential
    template <class data_t, template <typename> class vars_t>
    emtensor_t<data_t> compute_emtensor(
        const vars_t<data_t> &vars,          //!< the value of the variables
        const vars_t<Tensor<1, data_t>> &d1, //!< the value of the 1st derivs
        const Tensor<2, data_t> &h_UU, //!< the inverse metric (raised indices)
        const Tensor<3, data_t> &chris_ULL)
        const; //!< the conformal christoffel symbol

    //! The function which calculates the EM Tensor, given the vars and
    //! derivatives, excluding the potential
    template <class data_t, template <typename> class vars_t>
    static void emtensor_excl_potential(
        emtensor_t<data_t> &out,             //!< the em tensor output
        const vars_t<data_t> &vars,          //!< the value of the variables
        const vars_t<Tensor<1, data_t>> &d1, //!< the value of the first derivs
        const Tensor<2, data_t> &h_UU, //!< the inverse metric (raised indices).
        const Tensor<3, data_t>
            &chris_ULL); //!< the conformal christoffel symbol

    //! The function which adds in the RHS for the matter field vars,
    //! including the potential
    template <class data_t, template <typename> class vars_t,
              template <typename> class diff2_vars_t,
              template <typename> class rhs_vars_t>
    void add_matter_rhs(
        rhs_vars_t<data_t> &total_rhs,       //!< value of the RHS for all vars
        const vars_t<data_t> &vars,          //!< value of the variables
        const vars_t<Tensor<1, data_t>> &d1, //!< value of the 1st derivs
        const diff2_vars_t<Tensor<2, data_t>> &d2, //!< value of the 2nd derivs
        const vars_t<data_t> &advec)
        const; //!< the value of the advection terms

    //! The function which calculates the RHS for the matter field vars
    //! excluding the potential
    template <class data_t, template <typename> class vars_t,
              template <typename> class diff2_vars_t,
              template <typename> class rhs_vars_t>
    static void matter_rhs_excl_potential(
        rhs_vars_t<data_t> &rhs, //!< the value of the RHS terms for the sf vars
        const vars_t<data_t> &vars, //!< the values of all the variables
        const vars_t<Tensor<1, data_t>> &d1, //!< the value of the 1st derivs
        const diff2_vars_t<Tensor<2, data_t>> &d2, //!< value of the 2nd derivs
        const vars_t<data_t> &advec);

  protected:
    //! Adds the kinetic contribution of a single field to Sij, Si and rho
    //! (same expressions as ScalarField::emtensor_excl_potential). S is not
    //! touched.
    template <class data_t, template <typename> class vars_t>
    static void add_field_emtensor(emtensor_t<data_t> &out,
                                   const vars_t<data_t> &vars, const data_t &Pi,
                                   const Tensor<1, data_t> &d1_phi,
                                   const Tensor<2, data_t> &h_UU);

    //! Sets the RHS of a single field and its momentum excluding the potential
    //! (same expressions as ScalarField::matter_rhs_excl_potential)
    template <class data_t, template <typename> class vars_t>
    static void field_rhs_excl_potential(
        data_t &rhs_phi, data_t &rhs_Pi, const data_t &Pi,
        const Tensor<1, data_t> &d1_phi, const Tensor<2, data_t> &d2_phi,
        const data_t &advec_phi, const data_t &advec_Pi,
        const vars_t<data_t> &vars, const vars_t<Tensor<1, data_t>> &d1,
        const Tensor<2, data_t> &h_UU, const Tensor<3, data_t> &chris_ULL);
};

#include "TwoScalarFields.impl.hpp"

#endif /* TWOSCALARFIELDS_HPP_ */
