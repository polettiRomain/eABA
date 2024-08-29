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

Class
    Foam::RBD::ImposedMotion::qDriven

Group
    grpRigidBodyDynamicsImposedMotion

Description
    Linear spring imposedmotion.

SourceFiles
    qDriven.C

\*---------------------------------------------------------------------------*/

#ifndef RBD_ImposedMotion_qDriven_H
#define RBD_ImposedMotion_qDriven_H

#include "rigidBodyimposedmotion.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RBD
{
namespace ImposedMotion
{

/*---------------------------------------------------------------------------*\
                          Class qDriven Declaration
\*---------------------------------------------------------------------------*/

class qDriven
:
    public imposedmotion
{
    // Private data

        //- Anchor point, where the spring is attached to an immovable
        //  object
        point anchor_;

        //- Reference point of attachment to the solid body
        point refAttachmentPt_;

        //- Spring stiffness coefficient [N/m]
        scalar stiffness_;

        //- Damping coefficient [Ns/m]
        scalar damping_;

        //- Rest length - length of spring when no forces are applied to it
        scalar restLength_;


public:

    //- Runtime type information
    TypeName("qDriven");


    // Constructors

        //- Construct from components
        qDriven
        (
            const word& name,
            const dictionary& dict,
            const rigidBodyModel& model
        );

        //- Construct and return a clone
        virtual autoPtr<imposedmotion> clone() const
        {
            return autoPtr<imposedmotion>
            (
                new qDriven(*this)
            );
        }


    //- Destructor
    virtual ~qDriven();


    // Member Functions

        //- Accumulate the imposedmotion internal joint forces into the tau field and
        //  external forces into the q_doubleDot field
        virtual void imposedmotio
        (
          //  scalarField& tau,
            Field<spatialVector>& q_doubleDot,
            const rigidBodyModelState& state
        ) const;

        //- Update properties from given dictionary
        virtual bool read(const dictionary& dict);

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace ImposedMotion
} // End namespace RBD
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
