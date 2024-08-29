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

Namespace
    Foam::RBD::imposedmotion

Description
    Namespace for rigid-body dynamics imposedmotion

Class
    Foam::RBD::imposedmotion

Description
    Base class for defining imposedmotion for rigid-body dynamics

SourceFiles
    rigidBodyimposedmotion.C
    rigidBodyimposedmotionNew.C

\*---------------------------------------------------------------------------*/

#ifndef RBD_rigidBodyimposedmotion_H
#define RBD_rigidBodyimposedmotion_H

#include "dictionary.H"
#include "autoPtr.H"
#include "spatialVector.H"
#include "point.H"
#include "scalarField.H"
#include "runTimeSelectionTables.H"
#include "rigidBodyModelState.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RBD
{

// Forward declaration of classes
class rigidBodyModel;

/*---------------------------------------------------------------------------*\
                          Class imposedmotion Declaration
\*---------------------------------------------------------------------------*/

class imposedmotion
{

protected:

    // Protected data

        //- Name of the imposedmotion
        word name_;

        //- ID of the body the imposedmotion is applied to
        label bodyID_;

        //- Index of the body the force is applied to
        label bodyIndex_;

        //- imposedmotion model specific coefficient dictionary
        dictionary coeffs_;

        //- Reference to the model
        const rigidBodyModel& model_;

        //- Transform the given point on the restrained body to the global frame
        inline point bodyPoint(const point& p) const;

        //- Transform the velocity of the given point on the restrained body
        //  to the global frame
        inline spatialVector bodyPointVelocity(const point& p) const;


public:

    //- Runtime type information
    TypeName("imposedmotion");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            imposedmotion,
            dictionary,
            (
                const word& name,
                const dictionary& dict,
                const rigidBodyModel& model
            ),
            (name, dict, model)
        );


    // Constructors

        //- Construct from the dict dictionary and Time
        imposedmotion
        (
            const word& name,
            const dictionary& dict,
            const rigidBodyModel& model
        );

        //- Construct and return a clone
        virtual autoPtr<imposedmotion> clone() const = 0;


    // Selectors

        //- Select constructed from the dict dictionary and Time
        static autoPtr<imposedmotion> New
        (
            const word& name,
            const dictionary& dict,
            const rigidBodyModel& model
        );


    //- Destructor
    virtual ~imposedmotion();


    // Member Functions

        //- Return the name
        const word& name() const
        {
            return name_;
        }

        label bodyID() const
        {
            return bodyID_;
        }

        //- Accumulate the imposedmotion internal joint forces into the tau field and
        //  external forces into the q_doubleDot field
        virtual void imposedmotio
        (
          //  scalarField& tau,
            Field<spatialVector>& q_doubleDot,
              Field<spatialVector>& FD_type,
            const rigidBodyModelState& state
        ) const = 0;

        //- Update properties from given dictionary
        virtual bool read(const dictionary& dict);

        //- Return access to coeffs
        const dictionary& coeffDict() const;

        //- Write
        virtual void write(Ostream&) const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RBD
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "rigidBodyimposedmotionI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
