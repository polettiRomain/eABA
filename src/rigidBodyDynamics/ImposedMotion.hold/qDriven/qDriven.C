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

#include "qDriven.H"
#include "rigidBodyModel.H"
#include "addToRunTimeSelectionTable.H"
#include "Time.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace RBD
{
namespace ImposedMotion
{
    defineTypeNameAndDebug(qDriven, 0);

    addToRunTimeSelectionTable
    (
        imposedmotion,
        qDriven,
        dictionary
    );
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::RBD::ImposedMotion::qDriven::qDriven
(
    const word& name,
    const dictionary& dict,
    const rigidBodyModel& model
)
:
    imposedmotion(name, dict, model)
{
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::RBD::ImposedMotion::qDriven::~qDriven()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::RBD::ImposedMotion::qDriven::imposedmotio
(
  //  scalarField& tau,
    Field<spatialVector>& q_doubleDot,
    Field<spatialVector>& FD_type,
    const rigidBodyModelState& state
) const
{
    /*point attachmentPt = bodyPoint(refAttachmentPt_);

    // Current axis of the spring
    vector r = attachmentPt - anchor_;
    scalar magR = mag(r);
    r /= (magR + VSMALL);

    // Velocity of the attached end of the spring
    vector v = bodyPointVelocity(refAttachmentPt_).l();

    // Force and moment on the master body including optional damping
    vector force
    (
        (-stiffness_*(magR - restLength_) - damping_*(r & v))*r
    );

    vector moment(attachmentPt ^ force);*/

if (motionType == "HM")
{
 // SET ACCELERATION

      scalar q_doubleDot_alpha = A2*4*M_PI*M_PI*f2*f2*sin(2*M_PI*f2*model_.time().value()); //TRUE
      //scalar q_doubleDot_alpha = M_PI*M_PI*M_PI*sin(2*M_PI*model_.time().value()); //TRY
      //scalar q_doubleDot_alpha = 31.0063*cos(2*M_PI*model_.time().value()); //TRY
      scalar q_doubleDot_phi = 4*M_PI*M_PI*A1*f1*f1*cos(2*M_PI*f1*model_.time().value()); //TRUE
      //scalar q_doubleDot_phi = 39.4784*A1*cos(2*M_PI*(model_.time().value()+0.25));
      //scalar q_doubleDot_phi = 20.76*cos(2*M_PI*model_.time().value()); //TRY
      //vector twoD_acceleration = vector(q_doubleDot_alpha,q_doubleDot_phi,0); //	IF ONLY ROTATION
	    vector twoD_acceleration = vector(q_doubleDot_phi,q_doubleDot_alpha,0);

 // SET POSITION

    // scalar q_doubleDot_alpha = M_PI*M_PI*M_PI*sin(2*M_PI*model_.time().value()); //TRUE
    scalar q_alpha = -A2*sin(2*M_PI*f2*model_.time().value());
    //scalar q_alpha = -A2*cos(2*M_PI*model_.time().value())+A2;
    scalar q_phi = -A1*cos(2*M_PI*f1*model_.time().value())+A1; //TRUE
    //scalar q_phi = -A1*cos(2*M_PI*(model_.time().value()+0.25))+A1;
    //  scalar q_doubleDot_phi = -M_PI*M_PI*4*A1*sin(2*M_PI*model_.time().value()); //TRY
    //vector twoD_acceleration = vector(q_doubleDot_alpha,q_doubleDot_phi,0); //	IF ONLY ROTATION
    vector twoD_position = vector(q_phi,q_alpha,0);

// SET Velocity

  scalar qDot_alpha = -A2*2*M_PI*f2*cos(2*M_PI*f2*model_.time().value());
  scalar qDot_phi = A1*2*M_PI*f1*sin(2*M_PI*f1*model_.time().value()); //TRUE
  vector twoD_velocity = vector(qDot_phi,qDot_alpha,0);
  twoD_acceleration[2]=qDot_phi;
  twoD_position[2]=qDot_alpha;

    q_doubleDot  += spatialVector(twoD_acceleration, twoD_position);
}

else if (motionType == "Wang")
{  
scalar ti;
ti = model_.time().value();
 // SET ACCELERATION
      scalar q_doubleDot_alpha = 4 * M_PI*M_PI * f2*f2 *(A2/tanh(factortanh))*factortanh* 1/sqr(cosh(factortanh*sin(2*M_PI*f2*ti))) * (2*factortanh*sqr(cos(2*M_PI*f2*ti)) * tanh(factortanh*sin(2*M_PI*f2*ti)) + sin(2*M_PI*f2*ti));
      scalar q_doubleDot_phi = -(4*M_PI*M_PI*A1*pow(factortanh_phi,3)*f1*f1*sin(2*M_PI*f1*(ti+0.25/f1))* sqr(cos(2*M_PI*f1*(ti+0.25/f1))))/(asin(factortanh_phi)*pow(1-factortanh_phi*factortanh_phi*sqr(sin(2*M_PI*f1*(ti+0.25/f1))),1.5)) -  (4*M_PI*M_PI*A1*factortanh_phi*f1*f1*sin(2*M_PI*f1*(ti+0.25/f1)))/(asin(factortanh_phi)*pow(1-factortanh_phi*factortanh_phi*sqr(sin(2*M_PI*f1*(ti+0.25/f1))),0.5));        
      vector twoD_acceleration = vector(q_doubleDot_phi,q_doubleDot_alpha,0);          
 // SET POSITION
    scalar q_alpha = - A2/tanh(factortanh)*tanh(factortanh*sin(2*M_PI*f2*ti));
    scalar q_phi = -A1/asin(factortanh_phi) * asin(factortanh_phi*sin(2*M_PI*f1*(ti+0.25/f1))) +A1; 
    vector twoD_position = vector(q_phi,q_alpha,0);
// SET Velocity
  scalar qDot_alpha = -2*M_PI*f2*(A2/tanh(factortanh)) *factortanh* cos(2*M_PI*f2*ti)/sqr(cosh(factortanh*sin(2*M_PI*f2*ti)));
  scalar qDot_phi =  -A1/asin(factortanh_phi) * 2*M_PI*f1*factortanh_phi * cos(2*M_PI*f1*(ti+0.25/f1)) / (sqrt(1-(factortanh_phi*factortanh_phi)*(sqr(sin(2*M_PI*f1*(ti+0.25/f1))))));
  vector twoD_velocity = vector(qDot_phi,qDot_alpha,0);

  twoD_acceleration[2]=qDot_phi;
  twoD_position[2]=qDot_alpha;
  
  q_doubleDot  += spatialVector(twoD_acceleration, twoD_position);
  
   Info << "t = " << model_.time().value() << endl;
   Info << "qphi test " << q_phi << "--" << q_alpha << endl;
   Info << "qphi dot  test " << qDot_phi << "--" << qDot_alpha  << endl;
   Info << "qphi  ddot test " << q_doubleDot_phi << "--" << q_doubleDot_alpha << endl;

}
else if (motionType == "Wang_midStart")
{  
scalar ti;
ti = model_.time().value();
 // SET ACCELERATION
      scalar q_doubleDot_alpha = 4 * M_PI*M_PI * f2*f2 *(A2/tanh(factortanh))*factortanh* 1/sqr(cosh(factortanh*cos(2*M_PI*f2*ti))) * (2*factortanh*sqr(sin(2*M_PI*f2*ti)) * tanh(factortanh*cos(2*M_PI*f2*ti)) + cos(2*M_PI*f2*ti));
      //scalar q_doubleDot_phi = -(4*M_PI*M_PI*A1*pow(factortanh_phi,3)*f1*f1*sin(2*M_PI*f1*ti)* sqr(cos(2*M_PI*f1*ti)))/(asin(factortanh_phi)*pow(1-factortanh_phi*factortanh_phi*sqr(sin(2*M_PI*f1*ti)),1.5)) -  (4*M_PI*M_PI*A1*factortanh_phi*f1*f1*sin(2*M_PI*f1*ti))/(asin(factortanh_phi)*pow(1-factortanh_phi*factortanh_phi*sqr(sin(2*M_PI*f1*ti)),0.5));  

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
/*
Info << "2D acc" << twoD_acceleration << endl;
Info << "2D pos" << twoD_position << endl;
Info << "spat vect aft " << q_doubleDot << endl;
 */    
}



    //q_doubleDot  += spatialVector(Zero, Zero); // I could put += instead
    // of = if I make many qDrivens ad for restraints. In that case the force is a sum of
    //all the restraints, here I just have one qDriven so += means =, it is something to
    // keep in mind if I make mani qDrivens
    Info << "q_doubleDot qDriven" << q_doubleDot << endl;
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


bool Foam::RBD::ImposedMotion::qDriven::read
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


void Foam::RBD::ImposedMotion::qDriven::write
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
