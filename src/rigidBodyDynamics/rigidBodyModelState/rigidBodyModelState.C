/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2017 OpenFOAM Foundation
    Copyright (C) 2018-2020 OpenCFD Ltd.
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

#include "rigidBodyModelState.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::RBD::rigidBodyModelState::rigidBodyModelState
(
    const rigidBodyModel& model
)
:
    q_(model.nDoF(), Zero),
    qDot_(model.nDoF(), Zero),
    qDdot_(model.nDoF(), Zero),
    t_(-1),
    deltaT_(0)
{}


Foam::RBD::rigidBodyModelState::rigidBodyModelState
(
    const rigidBodyModel& model,
    const dictionary& dict
)
:
    q_(dict.getOrDefault("q", scalarField(model.nDoF(), Zero))),
    qDot_(dict.getOrDefault("qDot", scalarField(model.nDoF(), Zero))),
    qDdot_(dict.getOrDefault("qDdot", scalarField(model.nDoF(), Zero))),
    t_(dict.getOrDefault<scalar>("t", -1)),
    deltaT_(dict.getOrDefault<scalar>("deltaT", 0))

{

    if
    (
        q_.size() != model.nDoF()
     || qDot_.size() != model.nDoF()
     || qDdot_.size() != model.nDoF()
    )
    {
        FatalErrorInFunction << "State parameters 'q', 'qDot', 'qDdot'"
            << " do not have the same size as the number of DoF "
            << model.nDoF()
            << ". Is your \"rigidBodyMotionState\" state file consistent?"
            << exit(FatalError);
    }
}


// ************************************************************************* //
