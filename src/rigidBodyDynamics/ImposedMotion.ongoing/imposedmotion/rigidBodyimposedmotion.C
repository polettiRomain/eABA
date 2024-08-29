/*---------------------------------------------------------------------------* \
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

#include "rigidBodyimposedmotion.H"
#include "rigidBodyModel.H"
#include "joint.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace RBD
{
    defineTypeNameAndDebug(imposedmotion, 0);
    defineRunTimeSelectionTable(imposedmotion, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::RBD::imposedmotion::imposedmotion
(
    const word& name,
    const dictionary& dict,
    //const rigidBodyModel& model
    const joint& joint
)
:
    name_(name),
    coeffs_(dict),
    joint_(joint)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::RBD::imposedmotion::~imposedmotion()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

const Foam::dictionary& Foam::RBD::imposedmotion::coeffDict() const
{
    return coeffs_;
}


bool Foam::RBD::imposedmotion::read(const dictionary& dict)
{
    coeffs_ = dict;
    return true;
}


void Foam::RBD::imposedmotion::write(Ostream& os) const
{
    os.writeEntry("type", type());
    //os.writeEntry("body", model_.name(bodyID_));
}


// ************************************************************************* //
