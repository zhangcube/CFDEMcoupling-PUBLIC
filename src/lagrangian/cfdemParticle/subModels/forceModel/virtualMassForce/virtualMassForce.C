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

#include "virtualMassForce.H"
#include "addToRunTimeSelectionTable.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

#define NOTONCPU 9999

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(virtualMassForce, 0);

addToRunTimeSelectionTable
(
    forceModel,
    virtualMassForce,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
virtualMassForce::virtualMassForce
(
    const dictionary& dict,
    cfdemCloud& sm
)
:
    forceModel(dict,sm),
    propsDict_(dict.subDict(typeName + "Props")),
    velFieldName_(propsDict_.lookup("velFieldName")),
    U_(sm.mesh().lookupObject<volVectorField> (velFieldName_)),
    voidfractionFieldName_(propsDict_.lookup("voidfractionFieldName")),
    voidfraction_(sm.mesh().lookupObject<volScalarField> (voidfractionFieldName_)),
    UsFieldName_(propsDict_.lookup("granVelFieldName")),
    Us_(sm.mesh().lookupObject<volVectorField> (UsFieldName_)),
    phiFieldName_(propsDict_.lookup("phiFieldName")),
    phi_(sm.mesh().lookupObject<surfaceScalarField> (phiFieldName_)),
    UrelOldRegName_(typeName + "UrelOld"),
    UrelOld_(NULL),
    splitUrelCalculation_(propsDict_.lookupOrDefault<bool>("splitUrelCalculation",false)),
    useUs_(propsDict_.lookupOrDefault<bool>("useUs",false)),
    useFelderhof_(propsDict_.lookupOrDefault<bool>("useFelderhof",false)),
    Cadd_(0.5),
    Crho(1),
    DDtUrel_
    (   IOobject
        (
            "DDtUrel",
            sm.mesh().time().timeName(),
            sm.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        sm.mesh(),
        dimensionedVector("zero", dimensionSet(0,1,-2,0,0,0,0), vector(0,0,0))
    ),
    smoothingModel_
    (
        smoothingModel::New
        (
            propsDict_,
            sm
        )
    )
{
    // allocate UrelOld only if needed
    int UrelOldSize = 0;
    if(!splitUrelCalculation_ || !useUs_ )
        UrelOldSize = 3;
    
    // init force sub model
    setForceSubModels(propsDict_);
    // define switches which can be read from dict
    forceSubM(0).setSwitchesList(0,true); // activate treatExplicit switch
    forceSubM(0).setSwitchesList(1,true); // activate treatForceDEM switch (DEM side only treatment)
    forceSubM(0).setSwitchesList(4,true); // activate search for interpolate switch

    //set default switches (hard-coded default = false)
    forceSubM(0).setSwitches(0,true);  // will treat forces explicitly on CFD side - IMPORTANT!

    for (int iFSub=0;iFSub<nrForceSubModels();iFSub++)
        forceSubM(iFSub).readSwitches();

    //Extra switches/settings 
    if(propsDict_.found("Cadd"))
    {
        Cadd_ = readScalar(propsDict_.lookup("Cadd")); 
        Info << "Virtual mass model: using non-standard Cadd = " << Cadd_ << endl;
    }
    if(propsDict_.found("Crho"))
    {
        Crho = readScalar(propsDict_.lookup ("Crho"));
    }
    else
    {
        FatalError << "Virtual mass model: requires particle density Crho=particle density/1000" << abort(FatalError);
    }
    
    if(useUs_)
    {
        Info << "Virtual mass model: using averaged Us \n";
        Info << "WARNING: ignoring virtual mass of relative particle motion/collisions \n";
        
        if(splitUrelCalculation_)
        {
            FatalError << "Virtual mass model: useUs=true requires splitUrelCalculation_=false" << abort(FatalError);
        }
    }

    particleCloud_.checkCG(true);

    //Append the field names to be probed
    particleCloud_.probeM().initialize(typeName, typeName+".logDat");
    particleCloud_.probeM().vectorFields_.append("virtualMassForce"); //first entry must the be the force
    particleCloud_.probeM().vectorFields_.append("UrelOld");
    particleCloud_.probeM().scalarFields_.append("Vs");
    particleCloud_.probeM().scalarFields_.append("rho");
    particleCloud_.probeM().writeHeader();

    if(!splitUrelCalculation_)
        FatalError << "Virtual mass model: you have set 'splitUrelCalculation' to false, but this is not implemented. use true!" << abort(FatalError);

    if (particleCloud_.dataExchangeM().maxNumberOfParticles() > 0 && !splitUrelCalculation_)
    {
        // get memory for 2d array
        particleCloud_.dataExchangeM().allocateArray(UrelOld_,NOTONCPU,3);
        Info << "**Virtual mass model: allocating UrelOld " << endl;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

virtualMassForce::~virtualMassForce()
{
    if(!splitUrelCalculation_)
        particleCloud_.dataExchangeM().destroy(UrelOld_,3);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void virtualMassForce::setForce() const
{
    reAllocArrays();

    scalar dt = U_.mesh().time().deltaT().value();

    //Compute acceleration field
    if(splitUrelCalculation_)
        DDtUrel_ = fvc::ddt(U_) + fvc::div(phi_, U_); //Total Derivative of fluid velocity
    else if(useUs_)
        DDtUrel_ = fvc::ddt(U_) + fvc::div(phi_, U_) - fvc::ddt(Us_); //Total Derivative of fluid velocity minus average particle velocity
    
    // smoothen
    smoothingM().smoothen(DDtUrel_);

    interpolationCellPoint<vector> UInterpolator_(U_);
    interpolationCellPoint<vector> DDtUrelInterpolator_(DDtUrel_);
    interpolationCellPoint<scalar> voidfractionInterpolator_(voidfraction_);

    #include "setupProbeModel.H"

    bool haveUrelOld_(false); 

    for(int index = 0;index <  particleCloud_.numberOfParticles(); index++)
    {
            vector virtualMassForce(0,0,0);
        vector position(0,0,0);
        vector DDtUrel(0,0,0);

        scalar voidfraction(1);
        scalar sg(1);

            label cellI = particleCloud_.cellIDs()[index][0];

            if (cellI > -1) // particle Found
            {
            // particle position
                if(forceSubM(0).interpolation()) 
                position = particleCloud_.position(index);

            //********* acceleration value *********//
            if(splitUrelCalculation_ || useUs_ )
            {
                // DDtUrel from acceleration field 
                if(forceSubM(0).interpolation())
                    DDtUrel = DDtUrelInterpolator_.interpolate(position,cellI);
                else
                    DDtUrel = DDtUrel_[cellI];
                }
                else
                {
                // DDtUrel from UrelOld
                vector Ufluid(0,0,0);

                // relative velocity
                if(forceSubM(0).interpolation())
                    Ufluid = UInterpolator_.interpolate(position,cellI);
                else
                    Ufluid = U_[cellI];

                vector Us   = particleCloud_.velocity(index);
                vector Urel = Ufluid - Us;

                    //Check of particle was on this CPU the last step
                    if(UrelOld_[index][0]==NOTONCPU) //use 1. element to indicate that particle was on this CPU the last time step
                        haveUrelOld_ = false;
                    else
                        haveUrelOld_ = true;

                vector UrelOld(0.,0.,0.);
                vector ddtUrel(0.,0.,0.);
                    for(int j=0;j<3;j++)
                    {
                        UrelOld[j]         = UrelOld_[index][j];
                    UrelOld_[index][j] = Urel[j];
                    }
                    if(haveUrelOld_ ) //only compute force if we have old (relative) velocity
                    ddtUrel = (Urel-UrelOld)/dt;
                }

            //********* Cadd value *********//
            scalar rho  = forceSubM(0).rhoField()[cellI];
                scalar ds = 2*particleCloud_.radius(index);
                scalar Vs = ds*ds*ds*M_PI/6;
            
            scalar Cadd;
            
            if (useFelderhof_)
            {                    
                scalar epsilons(0);
                scalar logsg(0);

                if(forceSubM(0).interpolation())
                    voidfraction = voidfractionInterpolator_.interpolate(position,cellI);
                else
                    voidfraction = voidfraction_[cellI];

                sg       = Crho*1000 / rho;
                logsg    = log(sg);
                epsilons = 1-voidfraction;

                Cadd = 0.5
                        + ( 0.047*logsg + 0.13)*epsilons
                        + (-0.066*logsg - 0.58)*epsilons*epsilons
                        + (               1.42)*epsilons*epsilons*epsilons;
            }
            else
            {
                // use predefined value
                Cadd = Cadd_;
            }

            //********* calculate force *********//
            virtualMassForce = Cadd_ * rho * Vs * DDtUrel;

                if( forceSubM(0).verbose() ) //&& index>100 && index < 105)
                {
                    Pout << "index / cellI = " << index << tab << cellI << endl;
                    Pout << "position = " << particleCloud_.position(index) << endl;
                }

                //Set value fields and write the probe
                if(probeIt_)
                {
                    #include "setupProbeModelfields.H"
                    // Note: for other than ext one could use vValues.append(x)
                    // instead of setSize
                    vValues.setSize(vValues.size()+1, virtualMassForce);           //first entry must the be the force 
                    sValues.setSize(sValues.size()+1, Vs);
                    sValues.setSize(sValues.size()+1, rho);
                    particleCloud_.probeM().writeProbe(index, sValues, vValues);
                }
            }
    
            // write particle based data to global array
            forceSubM(0).partToArray(index,virtualMassForce,vector::zero);
            forceSubM(0).passvirtualMassForce(index,virtualMassForce);
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
void Foam::virtualMassForce::reAllocArrays() const
{
    if (!splitUrelCalculation_)
    {
        if(particleCloud_.numberOfParticlesChanged())
        {
            particleCloud_.dataExchangeM().allocateArray(UrelOld_,NOTONCPU,3);
            Info << "**Virtual mass model: allocating UrelOld " << endl;
        }

        // get DEM data
        particleCloud_.dataExchangeM().getData("UrelOld", "vector-atom", UrelOld_);
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
