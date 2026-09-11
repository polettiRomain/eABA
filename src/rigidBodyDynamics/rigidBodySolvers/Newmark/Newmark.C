/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016 OpenFOAM Foundation
    Copyright (C) 2020 OpenCFD Ltd.
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

#include "Newmark.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace RBD
{
namespace rigidBodySolvers
{
    defineTypeNameAndDebug(Newmark, 0);
    addToRunTimeSelectionTable(rigidBodySolver, Newmark, dictionary);
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::RBD::rigidBodySolvers::Newmark::Newmark
(
    rigidBodyMotion& body,
    const dictionary& dict
)
:
    rigidBodySolver(body),
    gamma_(dict.getOrDefault<scalar>("gamma", 0.5)),
    beta_
    (
        max
        (
            0.25*sqr(gamma_ + 0.5),
            dict.getOrDefault<scalar>("beta", 0.25)
        )
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::RBD::rigidBodySolvers::Newmark::~Newmark()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::RBD::rigidBodySolvers::Newmark::solve
(
    const scalarField& tau,
    const Field<spatialVector>& fx
)
{
    // Accumulate the restraint forces
    scalarField rtau(tau);
    Field<spatialVector> rfx(fx);
    model_.applyRestraints(rtau, rfx, state());

    // Load the imposed q, qDot, qDdot for the joints for which the motion is user-defined 
    Field<label>  indexImposedJoints(model_.nBodies(), Zero); // dummy init, parameter updated inside applyImposedMotion
    Field<scalar> imposedJoints(model_.nBodies(), Zero); // dummy init, parameter updated inside applyImposedMotion
    model_.applyImposedMotion(indexImposedJoints,imposedJoints);
    
    // 2. Set the joint motion to the "imposed" joints only
    int j = 0;
    forAll(indexImposedJoints, i)
    {	
         label qi = indexImposedJoints[i];
         //Info << qi << endl;
	 q()[qi] = imposedJoints[j];
	 qDot()[qi] = imposedJoints[j+1];
	 qDdot()[qi] = imposedJoints[j+2];
	 j+= 3;
    }
    
   //3. Freeze the passsive dof until a user-defined time is reached 
   scalar t = model_.time().value();
   if (t < 0.0125) 
	{
	    bool isImposed;
	    for (label i=0; i<model_.nDoF(); i++)
	    {
		    isImposed = false;
		    forAll(indexImposedJoints, j)
		    {
			if (indexImposedJoints[j] == i)
			{
				isImposed = true;
				break;
			}
		    }	
		    if (!isImposed) 
	 	    { 
	    		q()[i] =0;
	    		qDot()[i] =0;
	    		qDdot()[i] =0;
		    }
	    }
	}



    //4. Calculate the accelerations for the given state and forces
    model_.forwardDynamics(state(), rtau, rfx,indexImposedJoints,imposedJoints);
   
   //5. Integrate the passive dof after the user-defined time is reached
   if (t >= 0.0125)
   {
	    bool isImposed;
	    for (label i=0; i<model_.nDoF(); i++)
	    {
		    isImposed = false;
		    forAll(indexImposedJoints, j)
		    {
			if (indexImposedJoints[j] == i)
			{
				isImposed = true;
				break;
			}
		    }	
		    if (!isImposed) 
	 	    { 

	    		qDot()[i] = qDot0()[i] + deltaT()*(gamma_*qDdot()[i] + (1 - gamma_)*qDdot0()[i]);
	    		q()[i]    = q0()[i] + deltaT()*qDot0()[i] + sqr(deltaT())*(beta_*qDdot()[i] + (0.5 - beta_)*qDdot0()[i]);
		    }
	    }
    }



      Info << "qDdot_newmark= " << qDdot() << endl;
      Info << "qDot_newmark= " << qDot() << endl;
      Info << "q_newmark= " << q() << endl;

    correctQuaternionJoints();
}


// ************************************************************************* //
