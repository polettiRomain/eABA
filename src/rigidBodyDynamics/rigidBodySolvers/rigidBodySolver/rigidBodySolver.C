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

#include "rigidBodySolver.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace RBD
{

    defineTypeNameAndDebug(rigidBodySolver, 0);
    defineRunTimeSelectionTable(rigidBodySolver, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::RBD::rigidBodySolver::rigidBodySolver(rigidBodyMotion& body)
:
    model_(body)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::RBD::rigidBodySolver::~rigidBodySolver()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::RBD::rigidBodySolver::correctQuaternionJoints()
{

    if (model_.unitQuaternions())
    {
        forAll(model_.joints(), i)
        {
            const label qi = model_.joints()[i].qIndex();

            if (model_.joints()[i].unitQuaternion())
            {
                // Calculate the change in the unit quaternion
                vector dv((q().block<vector>(qi) - q0().block<vector>(qi)));
                scalar magDv = mag(dv);

                if (magDv > SMALL)
                {
                    // Calculate the unit quaternion corresponding to the change
                    quaternion dQuat(dv/magDv, cos(magDv), true);
                    Info << "ENTRA" << endl;
                    // Transform the previous time unit quaternion
                    quaternion quat
                    (
                        model_.joints()[i].unitQuaternion(q0())*dQuat
                    );
                    quat.normalise();

                    // Update the joint unit quaternion
                    model_.joints()[i].unitQuaternion(quat, q());
                }
            }
        }
    }
}


// ************************************************************************* //
