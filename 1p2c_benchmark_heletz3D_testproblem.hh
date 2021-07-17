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
/**
 * \file
 * \brief Definition of a problem, for the 1p2c problem:
 * Component transport of nitrogen dissolved in the water phase.
 */
#ifndef DUMUX_1P2C_BENCHMARK_PROBLEM_HH
#define DUMUX_1P2C_BENCHMARK_PROBLEM_HH

#if HAVE_UG
#include <dune/grid/io/file/dgfparser/dgfug.hh>
#endif
#include <dune/grid/io/file/dgfparser/dgfs.hh>
#include <dune/grid/io/file/dgfparser/dgfyasp.hh>

#include "../1p2ctracermodel/1p2cmodel.hh"
#include <dumux/implicit/common/implicitporousmediaproblem.hh>

#include <dumux/material/fluidsystems/h2on2liquidphasefluidsystem.hh>
#include "1p2c_benchmark_heletz3D_spatialparams.hh"

namespace Dumux
{

template <class TypeTag>
class OnePTwoCOutflowProblem;

namespace Properties
{
NEW_TYPE_TAG(OnePTwoCOutflowProblem, INHERITS_FROM(OnePTwoC));
NEW_TYPE_TAG(OnePTwoCOutflowBoxProblem, INHERITS_FROM(BoxModel, OnePTwoCOutflowProblem));
NEW_TYPE_TAG(OnePTwoCOutflowCCProblem, INHERITS_FROM(CCModel, OnePTwoCOutflowProblem));

// Set the grid type
SET_PROP(OnePTwoCOutflowProblem, Grid)
{
#if HAVE_UG
    typedef Dune::UGGrid<3> type;
#else
    typedef Dune::SGrid<3, 3> type;
    //typedef Dune::YaspGrid<2> type;
#endif
};

// Set the problem property
SET_PROP(OnePTwoCOutflowProblem, Problem)
{
    typedef Dumux::OnePTwoCOutflowProblem<TypeTag> type;
};

// Set fluid configuration
SET_PROP(OnePTwoCOutflowProblem, FluidSystem)
{ private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
public:
    typedef Dumux::FluidSystems::H2ON2LiquidPhase<Scalar, false> type;
};

// Set the spatial parameters
SET_TYPE_PROP(OnePTwoCOutflowProblem,
              SpatialParams,
              Dumux::OnePTwoCOutflowSpatialParams<TypeTag>);

// Define whether mole(true) or mass (false) fractions are used
SET_BOOL_PROP(OnePTwoCOutflowProblem, UseMoles, false);

// Enable velocity output
SET_BOOL_PROP(OnePTwoCOutflowProblem, VtkAddVelocity, true);

// Disable gravity
SET_BOOL_PROP(OnePTwoCOutflowProblem, ProblemEnableGravity, false);
}


/*!
 * \ingroup OnePTwoCBoxModel
 * \ingroup ImplicitTestProblems
 *
 * \brief Definition of a problem, for the 1p2c problem:
 * Nitrogen is dissolved in the water phase and
 * is transported with the water flow from the left side to the right.
 *
 * The model domain is 1 m times 1 m with a discretization length of 0.05 m
 * and homogeneous soil properties (\f$ \mathrm{K=10e-10, \Phi=0.4, \tau=0.28}\f$).
 * Initially the domain is filled with pure water.
 *
 * At the left side, a Dirichlet condition defines a nitrogen mole fraction
 * of 0.3 mol/mol.
 * The water phase flows from the left side to the right due to the applied pressure
 * gradient of 1e5 Pa/m. The nitrogen is transported with the water flow
 * and leaves the domain at the right boundary
 * where an outflow boundary condition is applied.
 * 
 * The model is able to use either mole or mass fractions. The property useMoles can be set to either true or false in the
 * problem file. Make sure that the according units are used in the problem setup. The default setting for useMoles is true.
 *
 * This problem uses the \ref OnePTwoCBoxModel model.
 *
 * To run the simulation execute the following line in shell:
 * <tt>./test_box1p2c -parameterFile ./test_box1p2c.input</tt> or 
 * <tt>./test_cc1p2c -parameterFile ./test_cc1p2c.input</tt>
 */
template <class TypeTag>
class OnePTwoCOutflowProblem : public ImplicitPorousMediaProblem<TypeTag>
{
    typedef ImplicitPorousMediaProblem<TypeTag> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, BoundaryTypes) BoundaryTypes;
    typedef typename GET_PROP_TYPE(TypeTag, TimeManager) TimeManager;

    // copy some indices for convenience
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
    enum {
        // world dimension
        dimWorld = GridView::dimensionworld
    };
    enum {
        // indices of the primary variables
        pressureIdx = Indices::pressureIdx,
        massOrMoleFracIdx = Indices::massOrMoleFracIdx
    };
    enum {
        // index of the transport equation
    	conti0EqIdx = Indices::conti0EqIdx,
        transportEqIdx = Indices::transportEqIdx
    };


    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GridView::Intersection Intersection;

    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;

    //! property that defines whether mole or mass fractions are used
        static const bool useMoles = GET_PROP_VALUE(TypeTag, UseMoles);

public:
    OnePTwoCOutflowProblem(TimeManager &timeManager, const GridView &gridView)
        : ParentType(timeManager, gridView)
        , eps_(1e-6)
    {
        //initialize fluid system
        FluidSystem::init();
        episodeEnd_ = 86400;

        name_ 		= GET_RUNTIME_PARAM_FROM_GROUP(TypeTag,
                                             std::string, 
                                             Problem, 
                                             Name);
        episodeEnd_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag,
        										Scalar,
        										Problem,
        										EpisodeEnd);
        temperature_			= GET_RUNTIME_PARAM(TypeTag, Scalar,
        						InitialConditions.temperature);
        distanceBetweenWells_	= GET_RUNTIME_PARAM(TypeTag, Scalar,
        						InitialConditions.distanceBetweenWells);
        injectionTime_H_ 		= GET_RUNTIME_PARAM(TypeTag, Scalar,
        						InitialConditions.injectionTime);
        waterInjRateVolH_ 		= GET_RUNTIME_PARAM(TypeTag, Scalar,
        						InitialConditions.waterInjRate);
        lensA_width_			= GET_RUNTIME_PARAM(TypeTag, Scalar,
        						SpatialParameters.lensA_width);
        lensW_width_			= GET_RUNTIME_PARAM(TypeTag, Scalar,
        						SpatialParameters.lensW_width);
        lensK_width_			= GET_RUNTIME_PARAM(TypeTag, Scalar,
        						SpatialParameters.lensK_width);
        lensShaleUp_width_		= GET_RUNTIME_PARAM(TypeTag, Scalar,
        						SpatialParameters.lensShaleUp_width);
        lensShaleDown_width_	= GET_RUNTIME_PARAM(TypeTag, Scalar,
        						SpatialParameters.lensShaleDown_width);


        permSoilSandStoneAvg_ = GET_RUNTIME_PARAM(TypeTag, Scalar,
        						SpatialParameters.permSoilSandStoneAvg);
        permSoilShaleAvg_ = GET_RUNTIME_PARAM(TypeTag, Scalar,
        						SpatialParameters.permSoilShaleAvg);
        porositySoilSandStoneAvg_ = GET_RUNTIME_PARAM(TypeTag, Scalar,
        						SpatialParameters.porositySoilSandStoneAvg);
        porosiySoilShaleAvg_ = GET_RUNTIME_PARAM(TypeTag, Scalar,
        						SpatialParameters.porositySoilShaleAvg);

        bBoxMin_ = this->bBoxMin();
        bBoxMax_ = this->bBoxMax();

        this->spatialParams().setDomainProperties(bBoxMin_,bBoxMax_);
        //stateing in the console whether mole or mass fractions are used
        if(!useMoles)
        {
        	std::cout<<"problem uses mass-fractions"<<std::endl;
        }
        else
        {
        	std::cout<<"problem uses mole-fractions"<<std::endl;
        }
        this->timeManager().startNextEpisode(episodeEnd_);
    }

    /*!
     * \name Problem parameters
     */
    // \{

    /*!
     * \brief The problem name.
     *
     * This is used as a prefix for files generated by the simulation.
     */
    const char *name() const
    {
        return name_.c_str();
    }

    
    void episodeEnd()
	{
		// Start new episode if episode is over and assign new boundary conditions
		//if(this->timeManager().episodeIndex() ==1 )

		if (this->timeManager().time()<100*86400)
		{
			this->timeManager().startNextEpisode(episodeEnd_);
			this->timeManager().setTimeStepSize(episodeEnd_);
		}
		else
		{
			this->timeManager().startNextEpisode(episodeEnd_);
			this->timeManager().setTimeStepSize(episodeEnd_);
		}


	}
    
    /*!
     * \brief Returns the temperature within the domain.
     *
     * This problem assumes a temperature of 20 degrees Celsius.
     */
    Scalar temperature() const
    { return 273.15 + 20; }; // in [K]

    // \}

    /*!
     * \name Boundary conditions
     */
    // \{

    /*!
     * \brief Specifies which kind of boundary condition should be
     *        used for which equation on a given boundary segment.
     *
     * \param values The boundary types for the conservation equations
     * \param globalPos The position for which the bc type should be evaluated
     */
    void boundaryTypesAtPos(BoundaryTypes &values, 
                            const GlobalPosition &globalPos) const
    {
    	values.setAllNeumann();
        if(
        		onLeftBoundary_(globalPos)
        		|| onRightBoundary_(globalPos)
        		|| onUpperBoundary_(globalPos)
        		)
        {
            values.setAllDirichlet();
        }
        else if (injectionWellHorizontal_(globalPos))
        {
        	values.setDirichlet(massOrMoleFracIdx);
        	values.setNeumann(pressureIdx);
        }
        else
        {
            values.setAllNeumann();
        }
        
        //outflow condition for the transport equation at right boundary
        if (
        		extractionWellHorizontal_(globalPos)
//        		|| onRightBoundary_(globalPos)
//        		|| onLeftBoundary_(globalPos)
//        		|| onUpperBoundary_(globalPos)
        		)
        {
            values.setOutflow(transportEqIdx);
        }
    }

    /*!
     * \brief Evaluate the boundary conditions for a dirichlet
     *        boundary segment.
     *
     * \param values The dirichlet values for the primary variables
     * \param globalPos The position for which the bc type should be evaluated
     *
     * For this method, the \a values parameter stores primary variables.
     */
    void dirichletAtPos(PrimaryVariables &values, const GlobalPosition &globalPos) const
    {
        initial_(values, globalPos);
        Scalar simTime  = this->timeManager().time();
        if (injectionWellHorizontal_(globalPos))
        {
//        	values[massOrMoleFracIdx] = 1.0e-3 / (1 * 0.2); /// 1 * 1e-3 [g/l]/(1 [g/l] * 0.2)
            if (simTime < injectionTime_ )
            	values[massOrMoleFracIdx] = 1;
            else
            	values[massOrMoleFracIdx] = 0;
        }

    }

    /*!
     * \brief Evaluate the boundary conditions for a Neumann
     *        boundary segment.
     *
     * For this method, the \a priVars parameter stores the mass flux
     * in normal direction of each component. Negative values mean
     * influx.
     *
     * The units must be according to either using mole or mass fractions. (mole/(m^2*s) or kg/(m^2*s))
     */
    void neumann(PrimaryVariables &priVars,
                 const Element &element,
                 const FVElementGeometry &fvGeometry,
                 const Intersection &intersection,
                 const int scvIdx,
                 const int boundaryFaceIdx) const
    {


    	const GlobalPosition &globalPos=element.geometry().corner(scvIdx);
        priVars = 0;
        if (injectionWellHorizontal_(globalPos))
        {
            GlobalPosition centerBoundaryVector = intersection.geometry().center()
                                        - element.geometry().center();
            Scalar injLength = 2 * fabs(centerBoundaryVector[1]);

    		priVars[pressureIdx] = -0.2 / injLength /2; //0.2 l/s*m = 0.2 kg/s (rho=1000) 17.29 m^3/day // To be divided by the injection length and 2 semi-problem
//     		priVars[massOrMoleFracIdx] = 0;
        }
    }


    // \}

    /*!
     * \name Volume terms
     */
    // \{

    /*!
     * \brief Evaluate the source term for all phases within a given
     *        sub-control-volume.
     *
     * For this method, the \a priVars parameter stores the rate mass
     * of a component is generated or annihilate per volume
     * unit. Positive values mean that mass is created, negative ones
     * mean that it vanishes.
     *
     * The units must be according to either using mole or mass fractions. (mole/(m^3*s) or kg/(m^3*s))
     */
    void sourceAtPos(PrimaryVariables &priVars,
                     const GlobalPosition &globalPos) const
    {
    	priVars = Scalar(0.0);

        Scalar waterInjRateMass = waterInjRateVol_ * 1000 /*volVars.density()*/; //TODO implementvolVars

        Scalar simTime        = this->timeManager().time();
        Scalar waitingTime 	  = 4 * 7 * 24 * 3600; /*sec*/ /*2weeks*/
        Scalar extractionTime = 6 * 7 * 24 * 3600; /*sec*/ /*2weeks*/


        if (injectionWellHorizontal_(globalPos)) {
        	priVars[conti0EqIdx] = waterInjRateMass ;// / totalInjectionVolume_; // [kg/(s*m^3)]
        	priVars[transportEqIdx] = tracerInjRate_; // / totalInjectionVolume_; // [mol/(s*m^3)]

            if (simTime < injectionTime_ )
            {// source term of the total mass
            	priVars[conti0EqIdx] = waterInjRateMass ; // / totalInjectionVolume_; // [kg/(s*m^3)]
            	priVars[transportEqIdx] = tracerInjRate_ ; // [mol/(s*m^3)] or [kg/(s*m^3)]
            }
//            else if (injectionTime <= simTime && simTime <= waitingTime )
//			{// source term of the total mass
//				priVars[conti0EqIdx] = waterInjRateMass / totalInjectionVolume_; // [kg/(s*m^3)]
//				priVars[transportEqIdx] = tracerInjRate_ / totalInjectionVolume_; // [mol/(s*m^3)]
//			}
//            else if (injectionTime <= waitingTime && simTime <= extractionTime )
//            {// source term of the total mass
//				priVars[conti0EqIdx] = waterInjRateMass / totalInjectionVolume_; // [kg/(s*m^3)]
//				priVars[transportEqIdx] = tracerInjRate_ / totalInjectionVolume_; // [mol/(s*m^3)]
//			}
            else
            {// source term of the total mass
//            	std::cout<<"simTime = " <<simTime<< "injTime = "<<injectionTime << std::endl;
				priVars[conti0EqIdx] = waterInjRateMass ; // / totalInjectionVolume_; // [kg/(s*m^3)]
				priVars[transportEqIdx] = 0 ;// / totalInjectionVolume_; // [mol/(s*m^3)]
			}

        }

        else if (extractionWellHorizontal_(globalPos))
        {
        	if (enableExtraction_){
        		priVars[conti0EqIdx] = -waterInjRateMass ; // / totalInjectionVolume_; // [kg/(s*m^3)]

//            priVars[conti0EqIdx] = - waterInjRateMass / totalInjectionVolume_; // [kg/(s*m^3)]
////            priVars[transportEqIdx] = - volVars.fluidState().massFraction(0, x1Idx)
////            						* volVars.density()
////            						* waterInjRateMass /totalInjectionVolume_; // [kg/(s*m^3)] //TODO injection volume is not the extraction volume

        	}
        	else
        	{
        		priVars[conti0EqIdx] = 0;
        		priVars[transportEqIdx] = 0;
        	}
        }

        else
        	priVars = 0;
    }

    /*!
     * \brief Evaluate the initial value for a control volume.
     *
     * \param values The initial values for the primary variables
     * \param globalPos The position for which the initial condition should be evaluated
     *
     * For this method, the \a values parameter stores primary
     * variables.
     */
    void initialAtPos(PrimaryVariables &values, const GlobalPosition &globalPos) const
    {
        initial_(values, globalPos);
    }

    // \}

    /*!
	 * \brief Called directly after the time integration.
	 */
    void computeInjectionArea() const
    {

    }

private:
    // the internal method for the initial condition
    void initial_(PrimaryVariables &priVars,
                  const GlobalPosition &globalPos) const
    {
        priVars[pressureIdx] = 1.0e5; // - 1e5*globalPos[0]; // initial condition for the pressure
        priVars[massOrMoleFracIdx] = 0.0;  // initial condition for the N2 molefraction
//    	if (injectionLocationMiddle(globalPos))
//    	{
//    		priVars[massOrMoleFracIdx] = 1.0e-3 * 1000*0.2; // 1 mg/l
//    	}
    }

    bool injectionLocationMiddle (const GlobalPosition &globalPos) const
    {
    	Scalar midPoint_X_ = (this->bBoxMax()[0] - this->bBoxMin()[0])/2;
    	Scalar midPoint_Y_ = (this->bBoxMax()[1] - this->bBoxMin()[1])/2;
    	return ((globalPos[0]> midPoint_X_ - eps_) && (globalPos[0]< midPoint_X_ + eps_)
    			&& (globalPos[1]> midPoint_Y_ - eps_) && (globalPos[1]< midPoint_Y_ + eps_));
    }


    // position of the injection well in the plane (bird eye) view
    bool injectionWellHorizontal_(const GlobalPosition &globalPos) const
    {
    	Scalar well_width_x_ = 0; // inject as point or as line
    	Scalar xPosInjectionWell = x_length_ / 2 - distanceBetweenWells_/2;
    	Scalar yPosInjectionWell = 0.0;
    	return
        	xPosInjectionWell - well_width_x_ - eps_< globalPos[0]
         && globalPos[0] < xPosInjectionWell + well_width_x_ + eps_
         && yPosInjectionWell - eps_ < globalPos[1]
         && globalPos[1] < yPosInjectionWell + eps_ ;
    }

    bool extractionWellHorizontal_(const GlobalPosition &globalPos) const
    {
    	Scalar well_width_x_ = 0; // inject as point or as line
    	Scalar xPosExtractionWell = 1;
    	Scalar yPosExtractionWell = 0.0;

        return
    			globalPos[0] > this->bBoxMin()[0] + xPosExtractionWell - well_width_x_ - eps_
    		 && globalPos[0] < this->bBoxMin()[0] + this->bBoxMin()[0]+ xPosExtractionWell + well_width_x_ + eps_
        	 && this->bBoxMin()[1] + yPosExtractionWell - eps_ < globalPos[1]
        	 && globalPos[1] < this->bBoxMin()[1] + yPosExtractionWell + eps_ ;
    };

    bool onLeftBoundary_(const GlobalPosition &globalPos) const
    {
        return globalPos[0] < this->bBoxMin()[0] + eps_;
    }

    bool onRightBoundary_(const GlobalPosition &globalPos) const
    {
        return globalPos[0] > this->bBoxMax()[0] - eps_;
    }

    bool onLowerBoundary_(const GlobalPosition &globalPos) const
    {
	    return globalPos[1] < this->bBoxMin()[1] + eps_;
    }

	bool onUpperBoundary_(const GlobalPosition &globalPos) const
	{
	    return globalPos[1] > this->bBoxMax()[1] - eps_;
	}

    Scalar eps_;
    std::string name_;
    Scalar episodeEnd_;

	Scalar temperature_;
	Scalar lensA_width_;
	Scalar lensW_width_;
	Scalar lensK_width_;
	Scalar permSoilSandStoneAvg_;
	Scalar permSoilShaleAvg_;
	Scalar porositySoilSandStoneAvg_;
	Scalar porosiySoilShaleAvg_;
	Scalar lensShaleUp_width_;
	Scalar lensShaleDown_width_;

    Scalar x_length_;
	Scalar y_length_;
	Scalar initialConcentration_;
	Scalar totalInjectionVolume_;
	Scalar injectionTime_;
	Scalar injectionTime_H_;
	Scalar distanceBetweenWells_;
	Scalar waterInjRateVol_, waterInjRateVolH_;
	Scalar tracerInjRate_;
	Scalar tracerInjConcentration_;

	bool continuousInjection_;
	bool enableExtraction_;
	Scalar globalIdxExtWell_;

	GlobalPosition bBoxMin_;
	GlobalPosition bBoxMax_;


};

} //end namespace
#endif
