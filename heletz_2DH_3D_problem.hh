// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2009 by Katherina Baber
 *   Copyright (C) 2009 by Karin Erbertseder                                 *
 *   Copyright (C) 2009 by Andreas Lauser                                    *
 *   Copyright (C) 2008 by Bernd Flemisch                                    *
 *   Institute for Modelling Hydraulic and Environmental Systems             *
 *   University of Stuttgart, Germany                                        *
 *   email: <givenname>.<name>@iws.uni-stuttgart.de                          *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/**
 * \file
 * \brief Definition of a problem, for the 1p2c box problem:
 * Component transport of nitrogen dissolved in the water phase.
 */
#ifndef DUMUX_1P2C_OUTFLOW_PROBLEM_HH
#define DUMUX_1P2C_OUTFLOW_PROBLEM_HH

#if HAVE_UG
#include <dune/grid/io/file/dgfparser/dgfug.hh>
#endif
#include <dune/grid/io/file/dgfparser/dgfs.hh>
#include <dune/grid/io/file/dgfparser/dgfyasp.hh>

#include <dumux/boxmodels/1p2c/1p2cmodel.hh>
#include <dumux/boxmodels/common/porousmediaboxproblem.hh>

#include <dumux/material/fluidsystems/h2on2liquidphasefluidsystem.hh>
#include "heletz_2DH_spatialparameters.hh"

namespace Dumux
{

template <class TypeTag>
class OnePTwoCOutflowProblem;

namespace Properties
{
NEW_TYPE_TAG(OnePTwoCOutflowProblem, INHERITS_FROM(BoxOnePTwoC));

// Set the grid type
SET_PROP(OnePTwoCOutflowProblem, Grid)
{
#if HAVE_UG
    typedef Dune::UGGrid<3> type;
#else
    typedef Dune::SGrid<2, 2> type;
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

//Define whether mole(true) or mass (false) fractions are used
SET_BOOL_PROP(OnePTwoCOutflowProblem, UseMoles, false);

// Disable gravity
SET_BOOL_PROP(OnePTwoCOutflowProblem, EnableGravity, false);
}


/*!
 * \ingroup OnePTwoCBoxModel
 * \ingroup BoxTestProblems
 *
 * \brief Definition of a problem, for the 1p2c box problem:
 * Nitrogen is dissolved in the water phase and
 * is transported with the water flow from the left side to the right.
 *
 * The model domain is 1m times 1m with a discretization length of 0.05m
 * and homogeneous soil properties (\f$ \mathrm{K=10e-10, \Phi=0.4}\f$).
 * Initially the domain is filled with pure water.
 *
 * At the left side, a Dirichlet condition defines a nitrogen mole fraction
 * of 0.3 mol/mol.
 * The water phase flows from the left side to the right due to the applied pressure
 * gradient of 1e5Pa/m. The nitrogen is transported with the water flow
 * and leaves the domain at the right boundary
 * where an outflow boundary condition is applied.
 * This problem uses the \ref OnePTwoCBoxModel.
 *
 * To run the simulation execute the following line in shell:
 * <tt>./test_1p2c -parameterFile ./test_1p2c.input</tt>
 */
template <class TypeTag>
class OnePTwoCOutflowProblem : public PorousMediaBoxProblem<TypeTag>
{
    typedef PorousMediaBoxProblem<TypeTag> ParentType;

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
        // Grid and world dimension
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld,
    };
    enum {
        // indices of the primary variables
        pressureIdx = Indices::pressureIdx,
        massOrMoleFracIdx = Indices::massOrMoleFracIdx,
    };
    enum {
        // indices of the equations
        conti0EqIdx = Indices::conti0EqIdx,
        transportEqIdx = Indices::transportEqIdx
    };


    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GridView::template Codim<0>::Iterator ElementIterator;
    typedef typename GridView::template Codim<dim>::Entity Vertex;
    typedef typename GridView::Intersection Intersection;

    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;

public:
    OnePTwoCOutflowProblem(TimeManager &timeManager, const GridView &gridView)
        : ParentType(timeManager, gridView)
        , eps_(1e-2)
    {

    	try
    	{
            temperature_	= GET_RUNTIME_PARAM(TypeTag, Scalar, InitialConditions.temperature);
            distanceBetweenWells_	= GET_RUNTIME_PARAM(TypeTag, Scalar, InitialConditions.distanceBetweenWells);
            injectionTime_H_ = GET_RUNTIME_PARAM(TypeTag, Scalar, InitialConditions.injectionTime);
            waterInjRateVolH_ = GET_RUNTIME_PARAM(TypeTag, Scalar, InitialConditions.waterInjRate);
            lensA_width_	= GET_RUNTIME_PARAM(TypeTag, Scalar, SpatialParameters.lensA_width);
            lensW_width_	= GET_RUNTIME_PARAM(TypeTag, Scalar, SpatialParameters.lensW_width);
            lensK_width_	= GET_RUNTIME_PARAM(TypeTag, Scalar, SpatialParameters.lensK_width);
            lensShaleUp_width_	= GET_RUNTIME_PARAM(TypeTag, Scalar, SpatialParameters.lensShaleUp_width);
            lensShaleDown_width_	= GET_RUNTIME_PARAM(TypeTag, Scalar, SpatialParameters.lensShaleDown_width);


            permSoilSandStoneAvg_ = GET_RUNTIME_PARAM(TypeTag, Scalar, SpatialParameters.permSoilSandStoneAvg);
            permSoilShaleAvg_ = GET_RUNTIME_PARAM(TypeTag, Scalar, SpatialParameters.permSoilShaleAvg);
            porositySoilSandStoneAvg_ = GET_RUNTIME_PARAM(TypeTag, Scalar, SpatialParameters.porositySoilSandStoneAvg);
            porosiySoilShaleAvg_ = GET_RUNTIME_PARAM(TypeTag, Scalar, SpatialParameters.porositySoilShaleAvg);

    	}
        catch (Dumux::ParameterException &e) {
            std::cerr << e << ". Abort!\n";
            exit(1) ;
        }
        //initialize fluid system
        FluidSystem::init();
        // calculate the injection volume
        totalInjectionVolume_ = 0;
//        distanceBetweenWells_ = 50; /*m*/
//        injectionTime_H_ 	  = 24; /*h*/
        injectionTime_ = injectionTime_H_ * 3600;

        bboxMin_ = this->bboxMin();
        bboxMax_ = this->bboxMax();

        std::cout<<"bboxMin_[2] = "<<bboxMin_[2]<<"\nbboxMax_[2]="<<bboxMax_[2]<<"\n";
    	x_length_ = bboxMax_[0] - bboxMin_[0];
    	y_length_ = bboxMax_[1] - bboxMin_[1];
        // total volumetric injection rate in m^3/h
//        waterInjRateVolH_  = 2.5;
    	waterInjRateVolH_ /=2; // only half of the domain is solved at a time
    	double injectionLength;
    	injectionLength = lensA_width_ + lensW_width_ + lensK_width_ +
    			lensShaleUp_width_ + lensShaleDown_width_; // width of the Heletz reservoir

        injectionVolume(gridView); // compute the Total Injection Volume = sum of SCV at the injection location
    	waterInjRateVol_ = waterInjRateVolH_ ;
    	waterInjRateVol_ /= totalInjectionVolume_;

    	if (dim == 2)
    	{
    		waterInjRateVol_ /= injectionLength; // / 10.6 integrated over the width
    	}
        waterInjRateVol_ *= 1.0 / 3600;  //transformation in seconds [m^3 / s]

        // tracer concentration in injected fluid in [mol/ml] or [g/l]
        tracerInjConcentration_  = 0.1; //kg/m^3
        tracerInjRate_ 		 = tracerInjConcentration_ * waterInjRateVol_; //[kg/s]

        // continuousInjection - tracer is injected during the whole simulation
        continuousInjection_ = false;
        // enable extraction at the extraction well
        enableExtraction_ 	 = false;


        /* output information about the control volume of the injection and length*/
        std::cout<<"totalInjectionVolume_ = "<<totalInjectionVolume_ <<"\n";
        std::cout<<"injectionLength = "<<injectionLength<<"\n";
        this->spatialParameters().setDomainProperties(bboxMin_,
        		bboxMax_, lensA_width_, lensW_width_, lensK_width_,
        		lensShaleUp_width_, lensShaleDown_width_,
        		permSoilSandStoneAvg_,permSoilShaleAvg_,
        		porositySoilSandStoneAvg_, porosiySoilShaleAvg_ );
        // Duration of the injection
        if (!continuousInjection_){
			std::cout<<"Tracer injection of for "<<	injectionTime_H_ << "hours\n";
			this->timeManager().startNextEpisode(/*length=*/injectionTime_);
        }
        else
        {
        	injectionTime_=1e10;
        	this->timeManager().startNextEpisode(/*length=*/injectionTime_);
        }
    }


    /*!
     * \brief Called by Dumux::TimeManager in order to do a time
     *        integration on the model.
     */
    void timeIntegration()
    {
        const int maxFails = 10;
        for (int i = 0; i < maxFails; ++i) {
            if (i > 0 && this->gridView().comm().rank() == 0)
                std::cout << "Newton solver did not converge. Retrying with 	time step of "
                          << this->timeManager().timeStepSize() << "sec\n";

            if (this->model().update(this->newtonMethod(),this-> newtonController()))
                return;

            // update failed
            Scalar dt = this->timeManager().timeStepSize();
            Scalar nextDt = dt / 2;
            this->timeManager().setTimeStepSize(nextDt);
        }

        DUNE_THROW(Dune::MathError,
                   "Newton solver didn't converge after "
                   << maxFails
                   << " time-step divisions. dt="
                   << this->timeManager().timeStepSize());
    }

    // Calculate the mass balance in the domain
    void postTimeStep()
    {
        PrimaryVariables tmp;
        this->model().globalStorage(tmp);
        std::cout << "Amount of conserved quantities in domain: " << tmp << "\n";
    }


    void injectionVolume(const GridView &gridView)
    {
        FVElementGeometry fvGeom;
        ElementIterator elemIt = gridView.template begin<0>();
        const ElementIterator endIt = gridView.template end<0>();
        for (; elemIt != endIt; ++ elemIt)
        {
            fvGeom.update(gridView, *elemIt);
            for (int i = 0; i < fvGeom.numVertices; ++i)
            {
                const GlobalPosition &pos = fvGeom.subContVol[i].global;
                if (injectionWellHorizontal_(pos))
                    totalInjectionVolume_ += fvGeom.subContVol[i].volume;
            };
        }
        totalInjectionVolume_ = this->gridView().comm().sum(totalInjectionVolume_);
    }

    void episodeEnd()
    {
    	this->timeManager().startNextEpisode(/*length=*/injectionTime_*2);
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
    	std::string out_baseName = "out_heletz1p2c_horizontal_dipole_";
    	std::stringstream ss;
    	Scalar wR = waterInjRateVolH_ * 2.0;

    	ss << wR;
    	std::string out_q = ss.str();
    	ss.str(""); //reset position in the stringstream
    	ss << tracerInjConcentration_;
    	std::string out_Cinit = ss.str();
    	ss.str(""); //reset position in the stringstream
    	ss << injectionTime_H_;
    	std::string out_T = ss.str();
    	std::string fileOutputName = out_baseName + "Q"+out_q
    			+ "_Cin" + out_Cinit +"_T" + out_T;

    	const char *charOutputName = fileOutputName.c_str();
    	return charOutputName;
//    	return "out_heletz1p2c_horizontal_dipole";
    }

    /*!
     * \brief Returns the temperature within the domain.
     *
     * This problem assumes a temperature of 20 degrees Celsius.
     */
    Scalar temperature() const
    {
    	return temperature_;
    }; // in [K]

    // \}

    /*!
     * \name Boundary conditions
     */
    // \{

    /*!
     * \brief Specifies which kind of boundary condition should be
     *        used for which equation on a given boundary segment.
     */
    void boundaryTypes(BoundaryTypes &values, const Vertex &vertex) const
    {
        const GlobalPosition globalPos = vertex.geometry().center();

        if(
        		onLeftBoundary_(globalPos)
        		|| onRightBoundary_(globalPos)
        		|| onUpperBoundary_(globalPos)
        		)
        {
//        	values.setOutflow(transportEqIdx);
            values.setAllDirichlet();
        }
        else if (injectionWellHorizontal_(globalPos))
        {
//        	values.setAllNeumann();
        	values.setDirichlet(transportEqIdx);
        }
        else
            values.setAllNeumann();

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
     * \brief Evaluate the boundary conditions for a Dirichlet
     *        boundary segment.
     *
     * For this method, the \a priVars parameter stores primary variables.
     */
    void dirichlet(PrimaryVariables &priVars, const Vertex &vertex) const
    {
        const GlobalPosition globalPos = vertex.geometry().center();
        Scalar simTime        = this->timeManager().time();
        initial_(priVars, globalPos);
        if (injectionWellHorizontal_(globalPos))
        {
            if (simTime < injectionTime_ )
            	priVars[massOrMoleFracIdx] = 1;
            else
            	priVars[massOrMoleFracIdx] = 0;
        }
        //condition for the N2 molefraction at left boundary
//        if(globalPos[0] < eps_)
//            priVars[massOrMoleFracIdx] = 2.0e-5;
    }

    /*!
     * \brief Evaluate the boundary conditions for a Neumann
     *        boundary segment.
     *
     * For this method, the \a priVars parameter stores the mass flux
     * in normal direction of each component. Negative values mean
     * influx.
     */
    void neumann(PrimaryVariables &priVars,
                 const Element &element,
                 const FVElementGeometry &fvGeometry,
                 const Intersection &is,
                 const int scvIdx,
                 const int boundaryFaceIdx) const
    {
        //const GlobalPosition &globalPos = element.geometry().corner(scvIdx);
        priVars = 0;
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
     * For this method, the \a priVars parameter stores primary
     * variables.
     */
    void initial(PrimaryVariables &priVars,
                 const Element &element,
                 const FVElementGeometry &fvGeometry,
                 const int scvIdx) const
    {
        const GlobalPosition &globalPos
            = element.geometry().corner(scvIdx);

        initial_(priVars, globalPos);
    }

    // \}

private:
    // the internal method for the initial condition
    void initial_(PrimaryVariables &priVars,
                  const GlobalPosition &globalPos) const
    {
        priVars[pressureIdx] = 1630*9.81*1e3; //- 1e5 * globalPos[0]/x_length_;//0.0; //initial condition for the pressure
        priVars[massOrMoleFracIdx] = 0.0; //initial condition for the N2 molefraction
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
    };

    bool injectionWellVertical_(const GlobalPosition &globalPos) const
	{
		Scalar xPosInjectionWell = x_length_ / 2;
		Scalar yPosInjectionWell = y_length_;

		return
				this->bboxMin()[0] + xPosInjectionWell - eps_< globalPos[0]
			            && globalPos[0] < this->bboxMin()[0] + xPosInjectionWell + eps_
//			            && yPosInjectionWell - eps_ < globalPos[1]
//						&& globalPos[1] < yPosInjectionWell + eps_
						;
	};


    bool extractionWellHorizontal_(const GlobalPosition &globalPos) const
    {
    	Scalar well_width_x_ = 0; // inject as point or as line
    	Scalar xPosExtractionWell = x_length_ / 2 + distanceBetweenWells_/2;
    	Scalar yPosExtractionWell = 0.0;
    	return
    			this->bboxMin()[0] + xPosExtractionWell - well_width_x_ - eps_< globalPos[0]
    		 && globalPos[0] < this->bboxMin()[0] + this->bboxMin()[0]+ xPosExtractionWell + well_width_x_ + eps_
        	 && this->bboxMin()[1] + yPosExtractionWell - eps_ < globalPos[1]
        	 && globalPos[1] < this->bboxMin()[1] + yPosExtractionWell + eps_ ;
    };


    bool onLeftBoundary_(const GlobalPosition &globalPos) const
    {
        return globalPos[0] < this->bboxMin()[0] + eps_;
    }

    bool onRightBoundary_(const GlobalPosition &globalPos) const
    {
        return globalPos[0] > this->bboxMax()[0] - eps_;
    }

    bool onLowerBoundary_(const GlobalPosition &globalPos) const
    {
	    return globalPos[1] < this->bboxMin()[1] + eps_;
    }

	bool onUpperBoundary_(const GlobalPosition &globalPos) const
	{
	    return globalPos[1] > this->bboxMax()[1] - eps_;
	}


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
	GlobalPosition bboxMin_;
	GlobalPosition bboxMax_;
	bool continuousInjection_;
	bool enableExtraction_;
	Scalar globalIdxExtWell_;
    const Scalar eps_;
};

} //end namespace
#endif
