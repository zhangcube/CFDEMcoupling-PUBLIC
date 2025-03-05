/*---------------------------------------------------------------------------*\
    CFDEMcoupling - Open Source CFD-DEM coupling

    CFDEMcoupling is part of the CFDEMproject
    www.cfdem.com
                                Christoph Goniva, christoph.goniva@cfdem.com
                                Copyright 2009-2012 JKU Linz
                                Copyright 2012-     DCS Computing GmbH, Linz
-------------------------------------------------------------------------------
License
    This file is part of CFDEMcoupling.

    CFDEMcoupling is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 3 of the License, or (at your
    option) any later version.

    CFDEMcoupling is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with CFDEMcoupling; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Description
    This code is designed to realize coupled CFD-DEM simulations using LIGGGHTS
    and OpenFOAM(R). Note: this code is not part of OpenFOAM(R) (see DISCLAIMER).
\*---------------------------------------------------------------------------*/

#include "error.H"

#include "SeepageForce.H"
#include "addToRunTimeSelectionTable.H"
    //添加渗透力
    //#include "capillarityModel.H"
    #include "relativePermeabilityModel.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(SeepageForce, 0);

addToRunTimeSelectionTable
(
    forceModel,
    SeepageForce,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
SeepageForce::SeepageForce
(
    const dictionary& dict,
    cfdemCloud& sm
)
:
    forceModel(dict,sm),
    propsDict_(dict.subDict(typeName + "Props")),
    VOFvoidfractionFieldName_(propsDict_.lookup("VOFvoidfractionFieldName")),
    alpha_(sm.mesh().lookupObject<volScalarField> (VOFvoidfractionFieldName_)),
    mu1(readScalar(propsDict_.lookup("mu1"))),
    mu2(readScalar(propsDict_.lookup("mu2"))),
    rho1(readScalar(propsDict_.lookup("rho1"))),
    rho2(readScalar(propsDict_.lookup("rho2"))),
    epslim_(readScalar(propsDict_.lookup("epslim"))),
    
    velFieldName_(propsDict_.lookup("velFieldName")),
    U_(sm.mesh().lookupObject<volVectorField> (velFieldName_)),
    voidfractionFieldName_(propsDict_.lookup("voidfractionFieldName")),
    voidfraction_(sm.mesh().lookupObject<volScalarField> (voidfractionFieldName_)),
    UsFieldName_(propsDict_.lookupOrDefault("granVelFieldName",word("Us"))),
    UsField_(sm.mesh().lookupObject<volVectorField> (UsFieldName_)),
    //relativepermeability
    KpermeabilityName_(propsDict_.lookup("KpermeabilityName")),
    Kpermeability(sm.mesh().lookupObject<volScalarField> (KpermeabilityName_)),

    interpolation_(false)
{


    // init force sub model
    setForceSubModels(propsDict_);

    // define switches which can be read from dict
    forceSubM(0).setSwitchesList(0,true); // activate search for treatExplicit switch
    forceSubM(0).setSwitchesList(2,true); // activate search for implDEM switch
    forceSubM(0).setSwitchesList(3,true); // activate search for verbose switch
    forceSubM(0).setSwitchesList(4,true); // activate search for interpolate switch
    forceSubM(0).setSwitchesList(7,true); // activate implForceDEMacc switch
    forceSubM(0).setSwitchesList(8,true); // activate scalarViscosity switch

    // read those switches defined above, if provided in dict
    for (int iFSub=0;iFSub<nrForceSubModels();iFSub++)
        forceSubM(iFSub).readSwitches();

    particleCloud_.checkCG(true);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

SeepageForce::~SeepageForce()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void SeepageForce::setForce() const
{
    //update force submodels to prepare for loop
    for (int iFSub=0;iFSub<nrForceSubModels();iFSub++)
        forceSubM(iFSub).preParticleLoop(forceSubM(iFSub).verbose());

    vector position(0,0,0);
    scalar voidfraction(1);
    vector Ufluid(0,0,0);
    label cellI=0;

    vector Us(0,0,0);
    vector Ur(0,0,0);
    scalar ds(0);
    scalar dParcel(0);
	scalar Vs(0);
    scalar Vcell(0);
    scalar Drag_(0);
    scalar piBySix(M_PI/6);

    //Updating Relative Permeability Model
volScalarField k0("k0", Kpermeability);

k0 = Kpermeability*(voidfraction*voidfraction*voidfraction)/(max((1-voidfraction)*(1-voidfraction),SMALL));

autoPtr<relativePermeabilityModel> krModel = relativePermeabilityModel::New("krModel",propsDict_,alpha_);
krModel->correct(); 
volScalarField kr1 = krModel->krb();
volScalarField kr2 = krModel->kra();

//Mobilities
volScalarField M1 ("M1",k0*kr1/mu1);	
volScalarField M2 ("M2",k0*kr2/mu2);	
volScalarField M ("M",M1+M2);


//Drag Coefficient Calculation
volScalarField Drag ("Drag", 1/M);

krModel->correct();
kr1=krModel->krb();
kr2=krModel->kra();

// Updating Mobilities
M1 = k0*kr1/mu1; 	
M2 = k0*kr2/mu2;
M = M1+M2;
Drag = 1/M;     

    #include "resetVoidfractionInterpolator.H"
    #include "resetUInterpolator.H"

    for(int index = 0;index <  particleCloud_.numberOfParticles(); index++)
    {
            scalar dp = 2*particleCloud_.radius(index);
            vector position = particleCloud_.position(index);
            cellI = particleCloud_.cellIDs()[index][0];
            Vs = 0;
            Vcell = 0;
            Ufluid =vector(0,0,0);
            voidfraction=0;
            Drag_=0;
            int n(0);

            if (cellI > -1) // particle Found
            {
                vector SeepageForceCell = Foam::vector(0,0,0);
                vector SeepageForce = Foam::vector(0,0,0);
                if(forceSubM(0).interpolation())
                {
	                position = particleCloud_.position(index);
                    voidfraction = voidfractionInterpolator_().interpolate(position,cellI);
                    Ufluid = UInterpolator_().interpolate(position,cellI);
                    Us = UsField_[cellI];
                    Drag_= Drag[cellI];
                }else
                {
					voidfraction = voidfraction_[cellI];
                    Ufluid = U_[cellI];
                    Us = UsField_[cellI];   //网格中颗粒速度
                    Drag_= Drag[cellI];
                }

                ds = particleCloud_.d(index);
                dParcel = ds;
                forceSubM(0).scaleDia(ds,index); //caution: this fct will scale ds!

                Vcell = particleCloud_.mesh().V()[cellI];
                //Us = particleCloud_.velocity(index);  //DEM中颗粒速度
                Ur = Ufluid-voidfraction*(1-voidfraction)*Us;
                Vs = ds*ds*ds*piBySix;

            if ( voidfraction < epslim_)
            {
                    SeepageForceCell = Drag_ * Vcell * Ur;
                    n = round((1-voidfraction)*Vcell/Vs);
                    SeepageForce = SeepageForceCell/n;
            }
            else
            {
                    SeepageForce = Foam::vector(0,0,0);
            }

               // write particle based data to global array
               forceSubM(0).partToArray(index,SeepageForce,vector::zero);
               forceSubM(0).passSeepageForce(index,SeepageForce);

        } // end if particle found on proc domain
        //}// end if in mask
    }// end loop particles
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
