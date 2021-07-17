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
 * \brief Definition of the spatial parameters for the 1p2c
 *        outlfow problem.
 */
#ifndef DUMUX_1P2C_OUTFLOW_SPATIAL_PARAMS_HH
#define DUMUX_1P2C_OUTFLOW_SPATIAL_PARAMS_HH

#include <dumux/material/spatialparams/implicitspatialparams1p.hh>
#include <dumux/material/fluidmatrixinteractions/2p/linearmaterial.hh>
#include <dumux/material/fluidmatrixinteractions/2p/regularizedbrookscorey.hh>
#include <dumux/material/fluidmatrixinteractions/2p/efftoabslaw.hh>

namespace Dumux
{

/*!
 * \ingroup OnePTwoCBoxModel
 * \ingroup ImplicitTestProblems
 *
 * \brief Definition of the spatial parameters for the 1p2c
 *        outflow problem.
 */
template<class TypeTag>
class OnePTwoCOutflowSpatialParams : public ImplicitSpatialParamsOneP<TypeTag>
{
    typedef ImplicitSpatialParamsOneP<TypeTag> ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, Grid) Grid;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;

    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GridView::Intersection Intersection;
    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;

    enum {
        // world dimension
        dimWorld = GridView::dimensionworld
    };

    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;

    typedef typename GET_PROP_TYPE(TypeTag, SolutionVector) SolutionVector;





    //typedef LinearMaterial<Scalar> EffMaterialLaw;
public:
    OnePTwoCOutflowSpatialParams(const GridView &gridView)
        : ParentType(gridView)
    {
//        Scalar Kh = 9.975e-5;
//        permeability_ = Kh * 0.001/(1000*9.81);

    	eps_ = 1e-6;
    	permeability_ = 1e-12;
        porosity_ = 0.2;
        dispersivity_ = 10.0;
    	permeability_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag,
    			Scalar,
    			Problem,
    			Permeability);


        dispersivity_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag,
        		Scalar,
        		Problem,
        		Dispersivity);

        porosity_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag,
        		Scalar,
        		Problem,
        		Porosity);
    	bool homogeneousCase = false;

        permSoilShaleAvg_ 		= GET_RUNTIME_PARAM(TypeTag, Scalar, SpatialParameters.permSoilShaleAvg);
        permSoilSandStoneAvg_ 	= GET_RUNTIME_PARAM(TypeTag, Scalar, SpatialParameters.permSoilSandStoneAvg);
        porosityShale_ 			= GET_RUNTIME_PARAM(TypeTag, Scalar, SpatialParameters.porositySoilShaleAvg);
        porositySoilSandStone_ 	= GET_RUNTIME_PARAM(TypeTag, Scalar, SpatialParameters.porositySoilSandStoneAvg);
        tortuosityShale_ 		= 0; //0.706;   // ???
        tortuositySandStone_ 	= 0; //0.280;   // ???
        lensA_width_ 			= GET_RUNTIME_PARAM(TypeTag, Scalar,SpatialParameters.lensA_width);
        lensW_width_ 			= GET_RUNTIME_PARAM(TypeTag, Scalar,SpatialParameters.lensW_width);
        lensK_width_ 			= GET_RUNTIME_PARAM(TypeTag, Scalar,SpatialParameters.lensK_width);  //1.2;
        lensWidthShaleDown_ 	= GET_RUNTIME_PARAM(TypeTag, Scalar, SpatialParameters.lensShaleDown_width);
        lensWidthShaleUp_   	= GET_RUNTIME_PARAM(TypeTag, Scalar, SpatialParameters.lensShaleUp_width);

//    	dispersivity_[0] = 1;
//    	dispersivity_[1] = 0.1;
    }

    ~OnePTwoCOutflowSpatialParams()
    {}


    /*!
     * \brief Update the spatial parameters with the flow solution
     *        after a timestep.
     *
     * \param globalSolution the global solution vector
     */
    void update(const SolutionVector &globalSolution)
    {
    };

    /*!
     * \brief Define the intrinsic permeability \f$\mathrm{[m^2]}\f$.
     *
     * \param element The current finite element
     * \param fvGeometry The current finite volume geometry of the element
     * \param scvIdx The index of the sub-control volume
     */
    const Scalar intrinsicPermeability(const Element &element,
                                       const FVElementGeometry &fvGeometry,
                                       const int scvIdx) const
    {
        const GlobalPosition &pos = fvGeometry.subContVol[scvIdx].global;
        if (isLowerShaleLens_(pos) || isUpperShaleLens_(pos))
            return permSoilShaleAvg_;
        else
            return permSoilSandStoneAvg_;
    }

    /*!
     * \brief Define the porosity \f$\mathrm{[-]}\f$.
     *
     * \param element The finite element
     * \param fvGeometry The finite volume geometry
     * \param scvIdx The local index of the sub-control volume where
     */
    double porosity(const Element &element,
                    const FVElementGeometry &fvGeometry,
                    const int scvIdx) const
    {
        const GlobalPosition &pos = fvGeometry.subContVol[scvIdx].global;
        if (isLowerShaleLens_(pos) || isUpperShaleLens_(pos))
            return porosityShale_;
        else
            return porositySoilSandStone_ ;
    }

    /*!
     * \brief Define the dispersivity.
     *
     * \param element The finite element
     * \param fvGeometry The finite volume geometry
     * \param scvIdx The local index of the sub-control volume where
     */
    double dispersivity(const Element &element,
                    const FVElementGeometry &fvGeometry,
                    const int scvIdx) const
    {
    	const GlobalPosition &pos = fvGeometry.subContVol[scvIdx].global;
        if (isLowerShaleLens_(pos) || isUpperShaleLens_(pos))
            return 1;
        else
            return dispersivity_;
    }

    bool useTwoPointGradient(const Element &element,
                             const int vertexI,
                             const int vertexJ) const
    {
        return false;
    }


    void setDomainProperties(GlobalPosition &bboxMin, GlobalPosition &bboxMax)
    {
    	bboxMin_ = bboxMin;
    	bboxMax_ = bboxMax;
    }

private:
    //Vertical Lenses of Shale Figure 2. D062
    bool isLowerShaleLens_(const GlobalPosition &globalPos) const
	{
    	Scalar lensPositionUp = bboxMin_[2] + lensA_width_ + lensWidthShaleDown_ ; /*m*/
    	Scalar lensPositionDown = bboxMin_[2] + lensA_width_ ; /*m*/
		if ( /*bboxMin_[0]-eps_ < globalPos[0] && globalPos[0] < bboxMax_[0]+eps_  &&*/
				lensPositionDown - eps_< globalPos[2] && globalPos[2] < lensPositionUp +eps_)
			return true;
		return false;
	}
    bool isUpperShaleLens_(const GlobalPosition &globalPos) const
	{
    	Scalar lensPositionUp = bboxMin_[2] + lensA_width_ + lensWidthShaleDown_ + lensK_width_ + lensWidthShaleUp_ ; /*m*/
    	Scalar lensPositionDown = bboxMin_[2] + lensA_width_ + lensWidthShaleDown_ + lensK_width_; /*m*/
		if (/*bboxMin_[0]-eps_ < globalPos[0] && globalPos[0] < bboxMax_[0]+eps_  &&*/
				lensPositionDown - eps_< globalPos[2] && globalPos[2] < lensPositionUp +eps_)
			return true;
		return false;
	}
    Scalar permeability_;
    Scalar porosity_;
    Scalar dispersivity_;

    Scalar lensWidthShaleDown_;
    Scalar lensWidthShaleUp_;
    Scalar lensA_width_;
    Scalar lensW_width_;
    Scalar lensK_width_;
    Scalar lensShaleUp_width_;
    Scalar lensShaleDown_width_;
    Scalar permSoilShaleAvg_;
    Scalar permSoilSandStoneAvg_;
    Scalar porosityShale_;
    Scalar porositySoilSandStone_ ;
    Scalar tortuosityShale_;
    Scalar tortuositySandStone_;

	GlobalPosition bboxMin_;
	GlobalPosition bboxMax_;
    Scalar eps_;
};

}

#endif
