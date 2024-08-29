/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016 OpenFOAM Foundation
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

#include "wangParametrization.H"
#include "joint.H"
#include "addToRunTimeSelectionTable.H"
#include "Time.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace RBD
{
namespace ImposedMotion
{
    defineTypeNameAndDebug(wangParametrization, 0);

    addToRunTimeSelectionTable
    (
        imposedmotion,
        wangParametrization,
        dictionary
    );
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::RBD::ImposedMotion::wangParametrization::wangParametrization
(
    const word& name,
    const dictionary& dict,
    const joint& joint
)
:
    imposedmotion(name, dict, joint)
{
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::RBD::ImposedMotion::wangParametrization::~wangParametrization()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::RBD::ImposedMotion::wangParametrization::imposedmotio
(
  //  scalarField& tau,
    Field<spatialVector>& q_doubleDot,
    Field<spatialVector>& FD_type,
    const rigidBodyModelState& state
) const
{


scalar ti;
//ti = model_.time().value();
// SET ACCELERATION
scalar q_doubleDot_alpha = 4 * M_PI*M_PI * f2*f2 *(A2/tanh(factortanh))*factortanh* 1/sqr(cosh(factortanh*cos(2*M_PI*f2*ti))) * (2*factortanh*sqr(sin(2*M_PI*f2*ti)) * tanh(factortanh*cos(2*M_PI*f2*ti)) + cos(2*M_PI*f2*ti));

scalar q_doubleDot_phi = -(A1 * factortanh_phi * 4*M_PI*M_PI*f1*f1 * sin(2*M_PI*f1*ti)* (sqr(factortanh_phi)*sqr(sin(2*M_PI*f1*ti)) + sqr(factortanh_phi)*sqr(cos(2*M_PI*f1*ti)) - 1)) / (asin(factortanh_phi) * pow((1 - sqr(factortanh_phi) * sqr(sin(2*M_PI*f1*ti))),1.5) );

vector twoD_acceleration = vector(q_doubleDot_phi,q_doubleDot_alpha,0);          
// SET POSITION
scalar q_alpha = -A2/tanh(factortanh)*tanh(factortanh*cos(2*M_PI*f2*ti)) + A2;
scalar q_phi = -A1/asin(factortanh_phi) * asin(factortanh_phi*sin(2*M_PI*f1*ti)); 
vector twoD_position = vector(q_phi,q_alpha,0);
// SET Velocity
scalar qDot_alpha = 2*M_PI*f2*(A2/tanh(factortanh)) *factortanh* sin(2*M_PI*f2*ti)/sqr(cosh(factortanh*sin(2*M_PI*f2*ti)));
scalar qDot_phi =  -A1/asin(factortanh_phi) * 2*M_PI*f1*factortanh_phi * cos(2*M_PI*f1*ti) / (sqrt(1-(factortanh_phi*factortanh_phi)*(sqr(sin(2*M_PI*f1*ti)))));
vector twoD_velocity = vector(qDot_phi,qDot_alpha,0);

twoD_acceleration[2]=qDot_phi;
twoD_position[2]=qDot_alpha;

//Info << "spat vect bef " << q_doubleDot << endl;

Info << "t = " << model_.time().value() << endl;
/*Info << "qphi test " << q_phi << "--" << q_alpha << endl;
Info << "qphi dot  test " << qDot_phi << "--" << qDot_alpha  << endl;
Info << "qphi  ddot test " << q_doubleDot_phi << "--" <<q_doubleDot_alpha << endl;

Info << "2D acc" << twoD_acceleration << endl;*/
q_doubleDot  += spatialVector(twoD_acceleration, twoD_position);



//q_doubleDot  += spatialVector(Zero, Zero); // I could put += instead
// of = if I make many wangParametrizations ad for restraints. In that case the force is a sum of
//all the restraints, here I just have one wangParametrization so += means =, it is something to
// keep in mind if I make mani wangParametrizations
Info << "q_doubleDot wangParametrization" << q_doubleDot << endl;
//q_doubleDot[bodyIndex_] with Bodyindex I have a structure of three or four spatial spatialVector
// with only the last vector, the one of rigidBodyMotion, used. The others are empty

//Info << "Fdt before " << FD_type << endl;
FD_type = spatialVector(FDtype,NfreeJoints,NcontrolledJoints,0,0,0);
//Info << "Fdt after " << FD_type << endl;


//value usefull for plot
  {  vector pos= model_.X0(bodyID_).r();
    Info << "pos" << model_.X0(bodyID_).r() << endl;
    vector Flapping_phi = model_.X0(bodyID_).E() & vector(0,1,0);
    vector Pitching_alpha = model_.X0(bodyID_).E() & vector(0,0,1);
    scalar FlappingAngle = acos(Flapping_phi[1]);
    scalar PitchingAlpha = acos(Pitching_alpha[2]);
    Info << "PitchingAngle(" << PitchingAlpha << endl;
    Info << "FlappingAngle(" << FlappingAngle << endl;
  }

}


bool Foam::RBD::ImposedMotion::wangParametrization::read
(
    const dictionary& dict
)
{
    imposedmotion::read(dict);
    
    coeffs_.readEntry("motionType", motionType);
    coeffs_.readEntry("factortanh_phi", factortanh_phi);
    coeffs_.readEntry("factortanh", factortanh);

    coeffs_.readEntry("Amplitude1", A1);
    coeffs_.readEntry("Frequency1", f1);
    coeffs_.readEntry("K1", K1);
    coeffs_.readEntry("phase1", phi1);

    coeffs_.readEntry("Amplitude2", A2);
    coeffs_.readEntry("Frequency2", f2);
    coeffs_.readEntry("K2", K2);
    coeffs_.readEntry("phase2", phi2);
    coeffs_.readEntry("FD_type", FDtype);
    coeffs_.readEntry("NfreeJoints", NfreeJoints);
    coeffs_.readEntry("NcontrolledJoints", NcontrolledJoints);

    return true;
}


void Foam::RBD::ImposedMotion::wangParametrization::write
(
    Ostream& os
) const
{
    imposedmotion::write(os);

    os.writeEntry("Amplitude1", A1);
    os.writeEntry("Frequency1", f1);
    os.writeEntry("K1", K1);
    os.writeEntry("phase1", phi1);
    os.writeEntry("Amplitude2", A2);
    os.writeEntry("Frequency2", f2);
    os.writeEntry("K2", K2);
    os.writeEntry("phase2", phi2);
    os.writeEntry("FD_type", FDtype);
    os.writeEntry("NfreeJoints", NfreeJoints);
    os.writeEntry("NcontrolledJoints", NcontrolledJoints);

}


// ************************************************************************* //
