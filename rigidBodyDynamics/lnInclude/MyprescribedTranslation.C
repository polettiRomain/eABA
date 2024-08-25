/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018-2021 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "MyprescribedTranslation.H"
#include "rigidBodyModel.H"
#include "addToRunTimeSelectionTable.H"
#include "Time.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace RBD
{
namespace restraints
{
    defineTypeNameAndDebug(MyprescribedTranslation, 0);

    addToRunTimeSelectionTable
    (
        restraint,
        MyprescribedTranslation,
        dictionary
    );
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::RBD::restraints::MyprescribedTranslation::MyprescribedTranslation
(
    const word& name,
    const dictionary& dict,
    const rigidBodyModel& model
)
:
    restraint(name, dict, model),

    omega_(Zero),
    oldForce_(Zero),
    error0_(Zero),
    integral0_(Zero)
{
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::RBD::restraints::MyprescribedTranslation::~MyprescribedTranslation()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::RBD::restraints::MyprescribedTranslation::restrain
(
    scalarField& tau,
    Field<spatialVector>& fx,
    const rigidBodyModelState& state
) const
{

    vector omega = model_.v(model_.master(bodyID_)).w();

    //scalar Inertia = mag(model_.I(model_.master(bodyID_)).Ic());

    // from the definition of the angular momentum:
    // moment = Inertia*ddt(omega)

   scalar Xref=-A_*cos(2*M_PI*model_.time().value())+A_;
   // Xref=(A/asin(K_))*asin(K_*sin(2*pi*model_.time().value()+1.5*M_PI))+A_;
   vector pos= model_.X0(bodyID_).r();

    vector error = Xref*vector(1,0,0)-pos;
    vector integral = integral0_ + error;
    vector derivative = (error - error0_);

    vector force = (p_*error + i_*integral + d_*derivative);
    force = relax_*force+(1-relax_)*oldForce_;

    force = (force&vector(1,0,0))*vector(1,0,0);

    // Accumulate the force for the restrained body
    fx[bodyIndex_] +=  spatialVector(Zero, force);

    oldForce_ = force;
    error0_ = error;
    integral0_ = integral;

    Info << "pos" << model_.X0(bodyID_).r() << endl;
    Info << "fx_MyprescribeTranslation = " << fx << endl;
    vector Orientation = model_.X0(bodyID_).E() & vector(1,0,0);
    scalar angle = acos(Orientation[0]);
    Info << "angle(" << angle << endl;

}


bool Foam::RBD::restraints::MyprescribedTranslation::read
(
    const dictionary& dict
)
{
    restraint::read(dict);

    refQ_ = coeffs_.getOrDefault<tensor>("referenceOrientation", I);

    if (mag(mag(refQ_) - sqrt(3.0)) > ROOTSMALL)
    {
        FatalErrorInFunction
            << "referenceOrientation " << refQ_ << " is not a rotation tensor. "
            << "mag(referenceOrientation) - sqrt(3) = "
            << mag(refQ_) - sqrt(3.0) << nl
            << exit(FatalError);
    }





    coeffs_.readEntry("relax", relax_);

    coeffs_.readEntry("p", p_);
    coeffs_.readEntry("i", i_);
    coeffs_.readEntry("d", d_);
    coeffs_.readEntry("K", K_);
    coeffs_.readEntry("Amplitude", A_);
    coeffs_.readEntry("phase", phi);



    // Read the actual entry


    return true;
}


void Foam::RBD::restraints::MyprescribedTranslation::write
(
    Ostream& os
) const
{
    restraint::write(os);

    os.writeEntry("referenceOrientation", refQ_);
  //  os.writeEntry("axis", axis_);

}


// ************************************************************************* //
