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

#include "symplectic.H"
#include "addToRunTimeSelectionTable.H"
#include "rigidBodyimposedmotion.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace RBD
{
namespace rigidBodySolvers
{
    defineTypeNameAndDebug(symplectic, 0);
    addToRunTimeSelectionTable(rigidBodySolver, symplectic, dictionary);
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::RBD::rigidBodySolvers::symplectic::symplectic
(
    rigidBodyMotion& body,
    const dictionary& dict
)
:
    rigidBodySolver(body)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::RBD::rigidBodySolvers::symplectic::~symplectic()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::RBD::rigidBodySolvers::symplectic::solve
(
    const scalarField& tau,
    const Field<spatialVector>& fx
)
{
    // First simplectic step:
    //     Half-step for linear and angular velocities
    //     Update position and orientation
    qDot() = qDot0() + 0.5*deltaT0()*qDdot();
    q() = q0() + deltaT()*qDot();

    // 1. Load the imposed q, qDot, qDdot for the active joints
    Field<label>  indexImposedJoints(model_.nBodies(), Zero); // dummy init, parameter updated inside applyImposedMotion
    Field<scalar> imposedJoints(model_.nBodies(), Zero);     // dummy init, parameter updated inside applyImposedMotion
    model_.applyImposedMotion(indexImposedJoints,imposedJoints);

    // 2. Set the joint motion to the "imposed" joints only
    int j = 0;
    forAll(indexImposedJoints, i)
    {	
         label qi = indexImposedJoints[i];
         Info << qi << endl;
	 q()[qi] = imposedJoints[j];
	 qDot()[qi] = imposedJoints[j+1];
	 qDdot()[qi] = imposedJoints[j+2];
	 j+= 3;
    }

   scalar t = model_.time().value();
   // freeze the passsive dof until a user-defined time is reached
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

    /*Info << "qDdot_sympl bef fwd = " << qDdot() << endl;
    Info << "qDot_sympl  = " << qDot() << endl;
    Info << "q_sympl  = " << q() << endl;*/

    correctQuaternionJoints();

    // Update the body-state prior to the evaluation of the restraints
    model_.forwardDynamicsCorrection(state());

    // Accumulate the restraint forces
    scalarField rtau(tau);
    Field<spatialVector> rfx(fx);
    model_.applyRestraints(rtau, rfx, state());

    // Calculate the body acceleration for the given state
    // and restraint forces
    model_.forwardDynamics(state(), rtau, rfx, indexImposedJoints, imposedJoints);

    // Second simplectic step:
    // Complete update of linear and angular velocities
    // Ensure that is applied only to the free joints
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
	    		qDot()[i] += 0.5*deltaT()*qDdot()[i];
		    }
	    }
    }

    Info << "qDdot_sympl = " << qDdot() << endl;
    Info << "qDot_sympl  = " << qDot() << endl;
    Info << "q_sympl  = " << q() << endl;

}


// ************************************************************************* //
