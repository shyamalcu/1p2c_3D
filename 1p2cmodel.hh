// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the            *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \file
 *
 * \brief Base class for all models which use the single-phase,
 *        two-component fully implicit model.
 *        Adaption of the fully implicit scheme to the one-phase two-component flow model.
 */

#ifndef DUMUX_ONEP_TWOC_MODEL_HH
#define DUMUX_ONEP_TWOC_MODEL_HH

#include <dumux/implicit/common/implicitvelocityoutput.hh>
#include "1p2cproperties.hh"

namespace Dumux
{

/*!
 * \ingroup OnePTwoCBoxModel
 * \brief Adaption of the fully implicit scheme to the one-phase two-component flow model.
 *
 * This model implements a one-phase flow of a compressible fluid, that consists of two components,
 * using a standard Darcy
 * approach as the equation for the conservation of momentum:
 \f[
 v = - \frac{\textbf K}{\mu}
 \left(\textbf{grad}\, p - \varrho {\textbf g} \right)
 \f]
 *
 * Gravity can be enabled or disabled via the property system.
 * By inserting this into the continuity equation, one gets
 \f[
 \phi\frac{\partial \varrho}{\partial t} - \text{div} \left\{
   \varrho \frac{\textbf K}{\mu}  \left(\textbf{grad}\, p - \varrho {\textbf g} \right)
 \right\} = q \;,
 \f]
 *
 * The transport of the components \f$\kappa \in \{ w, a \}\f$ is described by the following equation:
 \f[
 \phi \frac{ \partial \varrho X^\kappa}{\partial t}
 - \text{div} \left\lbrace \varrho X^\kappa \frac{{\textbf K}}{\mu} \left( \textbf{grad}\, p -
 \varrho {\textbf g} \right)
 + \varrho D^\kappa_\text{pm} \frac{M^\kappa}{M_\alpha} \textbf{grad} x^\kappa \right\rbrace = q.
 \f]
 *
 * All equations are discretized using a vertex-centered finite volume (box)
 * or cell-centered finite volume scheme as spatial
 * and the implicit Euler method as time discretization.
 * The model is able to use either mole or mass fractions. The property useMoles can be set to either true or false in the
 * problem file. Make sure that the according units are used in the problem setup. useMoles is set to true by default.
 *
 * The primary variables are the pressure \f$p\f$ and the mole or mass fraction of dissolved component \f$x\f$.
 */

template<class TypeTag >
class OnePTwoCModel : public GET_PROP_TYPE(TypeTag, BaseModel)
{
    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables) ElementVolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, SolutionVector) SolutionVector;

    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    enum {
        dim = GridView::dimension
    };
    typedef typename GridView::template Codim<0>::Iterator ElementIterator;

    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
    enum { phaseIdx = Indices::phaseIdx };

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;

    enum { isBox = GET_PROP_VALUE(TypeTag, ImplicitIsBox) };
    enum { dofCodim = isBox ? dim : 0 };

public:
    /*!
     * \brief \copybrief ImplicitModel::addOutputVtkFields
     *
     * Specialization for the OnePTwoCModel, adding pressure,
     * mass and mole fractions, and the process rank to the VTK writer.
     */
    template<class MultiWriter>
    void addOutputVtkFields(const SolutionVector &sol,
                            MultiWriter &writer)
    {
        typedef Dune::BlockVector<Dune::FieldVector<double, 1> > ScalarField;
        typedef Dune::BlockVector<Dune::FieldVector<double, dim> > VectorField;

        // create the required scalar fields
        unsigned numDofs = this->numDofs();
        ScalarField &pressure = *writer.allocateManagedBuffer(numDofs);
        ScalarField &delp = *writer.allocateManagedBuffer(numDofs);
        ScalarField &moleFraction0 = *writer.allocateManagedBuffer(numDofs);
        ScalarField &moleFraction1 = *writer.allocateManagedBuffer(numDofs);
        ScalarField &massFraction0 = *writer.allocateManagedBuffer(numDofs);
        ScalarField &massFraction1 = *writer.allocateManagedBuffer(numDofs);
        ScalarField &rho = *writer.allocateManagedBuffer(numDofs);
        ScalarField &mu = *writer.allocateManagedBuffer(numDofs);
        ScalarField &dispersivity_x = *writer.allocateManagedBuffer(numDofs);
        ScalarField &dispersivity_y = *writer.allocateManagedBuffer(numDofs);
        ScalarField &permeability = *writer.allocateManagedBuffer(numDofs);
        ScalarField &porosity = *writer.allocateManagedBuffer(numDofs);

        ScalarField &concentration0  = *writer.allocateManagedBuffer(numDofs);
        ScalarField &concentration1  = *writer.allocateManagedBuffer(numDofs);
        VectorField *velocity = writer.template allocateManagedBuffer<double, dim>(numDofs);
        ImplicitVelocityOutput<TypeTag> velocityOutput(this->problem_());

        if (velocityOutput.enableOutput())
        {
            // initialize velocity field
            for (unsigned int i = 0; i < numDofs; ++i)
            {
                (*velocity)[i] = Scalar(0);
            }
        }

        unsigned numElements = this->gridView_().size(0);
        ScalarField &rank = *writer.allocateManagedBuffer(numElements);
        ScalarField &porosityElem = *writer.allocateManagedBuffer(numElements);
        ScalarField &permElem = *writer.allocateManagedBuffer(numElements);


        ElementIterator eIt = this->gridView_().template begin<0>();
        ElementIterator eEndIt = this->gridView_().template end<0>();
        for (; eIt != eEndIt; ++eIt)
        {
            int eIdx = this->problem_().model().elementMapper().map(*eIt);
            rank[eIdx] = this->gridView_().comm().rank();

            FVElementGeometry fvGeometry;
            fvGeometry.update(this->gridView_(), *eIt);

            ElementVolumeVariables elemVolVars;
            elemVolVars.update(this->problem_(),
                               *eIt,
                               fvGeometry,
                               false /* oldSol? */);


            for (int scvIdx = 0; scvIdx < fvGeometry.numScv; ++scvIdx)
            {
                int globalIdx = this->dofMapper().map(*eIt, scvIdx, dofCodim);

                pressure[globalIdx] = elemVolVars[scvIdx].pressure();
                delp[globalIdx] = elemVolVars[scvIdx].pressure() - 1e5;
                moleFraction0[globalIdx] = elemVolVars[scvIdx].moleFraction(0);
                moleFraction1[globalIdx] = elemVolVars[scvIdx].moleFraction(1);
                massFraction0[globalIdx] = elemVolVars[scvIdx].massFraction(0);
                massFraction1[globalIdx] = elemVolVars[scvIdx].massFraction(1);
                concentration0[globalIdx] = elemVolVars[scvIdx].concentration(0);
                concentration1[globalIdx] = elemVolVars[scvIdx].concentration(1);
                rho[globalIdx] = elemVolVars[scvIdx].density();
                mu[globalIdx] = elemVolVars[scvIdx].viscosity();
                dispersivity_x[globalIdx] = elemVolVars[scvIdx].dispersivity()[0];
                dispersivity_y[globalIdx] = elemVolVars[scvIdx].dispersivity()[1];
                permeability[globalIdx] = this->problem_().spatialParams()
                		.intrinsicPermeability(*eIt, fvGeometry, scvIdx);
                porosity[globalIdx] = this->problem_().spatialParams()
                		.porosity(*eIt, fvGeometry, scvIdx);
                porosityElem[eIdx] = this->problem_().spatialParams()
                            		.porosity(*eIt, fvGeometry, scvIdx);
                permElem[eIdx] = this->problem_().spatialParams()
                            		.intrinsicPermeability(*eIt, fvGeometry, scvIdx);


            }

            // velocity output
            velocityOutput.calculateVelocity(*velocity, elemVolVars, fvGeometry, *eIt, phaseIdx);
        }

        writer.attachDofData(pressure, "P", isBox);
        writer.attachDofData(delp, "delp", isBox);
        if (velocityOutput.enableOutput())
        {
            writer.attachDofData(*velocity,  "velocity", isBox, dim);
        }
        char nameMoleFraction0[42], nameMoleFraction1[42];
        snprintf(nameMoleFraction0, 42, "x_%s", FluidSystem::componentName(0));
        snprintf(nameMoleFraction1, 42, "x_%s", FluidSystem::componentName(1));
        writer.attachDofData(moleFraction0, nameMoleFraction0, isBox);
        writer.attachDofData(moleFraction1, nameMoleFraction1, isBox);

        char nameMassFraction0[42], nameMassFraction1[42];
        snprintf(nameMassFraction0, 42, "X_%s", FluidSystem::componentName(0));
        snprintf(nameMassFraction1, 42, "X_%s", FluidSystem::componentName(1));
        writer.attachDofData(massFraction0, nameMassFraction0, isBox);
        writer.attachDofData(massFraction1, nameMassFraction1, isBox);
        writer.attachDofData(rho, "rho", isBox);
        writer.attachDofData(mu, "mu", isBox);
        writer.attachDofData(dispersivity_x, "dispersivity_x", isBox);
        writer.attachDofData(dispersivity_y, "dispersivity_y", isBox);
        writer.attachDofData(permeability, "permeability", isBox);
        writer.attachDofData(porosity, "porosity", isBox);
        writer.attachDofData(concentration0, "concentration0", isBox);
        writer.attachDofData(concentration1, "concentration1", isBox);
        writer.attachCellData(rank, "process rank");
        writer.attachCellData(porosityElem, "porosity");
        writer.attachCellData(permElem, "permeability");
    }
};
}

#include "1p2cpropertydefaults.hh"

#endif
