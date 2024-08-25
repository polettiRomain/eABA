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
    Foam::RBD::ImposedMotion::wangParametrization

Group
    grpRigidBodyDynamicsImposedMotion

Description
    Linear spring imposedmotion.

SourceFiles
    wangParametrization.C

\*---------------------------------------------------------------------------*/

#ifndef RBD_ImposedMotion_wangParametrization_H
#define RBD_ImposedMotion_wangParametrization_H

#include "rigidBodyimposedmotion.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RBD
{
namespace ImposedMotion
{

/*---------------------------------------------------------------------------*\
                          Class wangParametrization Declaration
\*---------------------------------------------------------------------------*/
scalar thold_=0;
scalar Aphi_PID=-1.096;
scalar beta_PID=0;
scalar A_phi_off=0;
int flag_init_=0;
class wangParametrization
:
    public imposedmotion
{
    // Private data

	Field<label> jointList_;
	string motionType_;

        scalar Aphi_;
        scalar fphi_;
        scalar Kphi_;

        scalar Aalpha_;
        scalar falpha_;
        scalar Kalpha_;

        scalar Kpz_;
        scalar Kdz_;
        scalar Kiz_;

	scalar Kpxz_;
        scalar Kdxz_;
        scalar Kpthz_;
        scalar Kdthz_;

	scalar biasz_;
	scalar biasx_;
	scalar biasth_;
	
	scalar Aphi_max_;	
	scalar Aphi_min_;
	scalar beta_max_;	
	scalar beta_min_;
	scalar theta_max_;	
	scalar theta_min_;
	
        scalar Kpx_;
        scalar Kdx_;
        scalar Kix_;
        
	scalar Kpzx_;
        scalar Kdzx_;
        scalar Kpthx_;
        scalar Kdthx_;
        
        scalar Kptheta_;
        scalar Kdtheta_;
        scalar Kitheta_;
        
        scalar Kpzth_;
        scalar Kdzth_;
        scalar Kpxth_;
        scalar Kdxth_;

	scalar href_;
	scalar xref_;
	scalar hdref_;
	scalar xdref_;
	scalar thetaref_;
	scalar thetadref_;

public:

    //- Runtime type information
    TypeName("wangParametrization");

    // Constructors

        //- Construct from components
        wangParametrization
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
                new wangParametrization(*this)
            );
        }


    //- Destructor
    virtual ~wangParametrization();


    // Member Functions

        // Compute the q, qDot and qDdot of the joints with indices in jointIndex
	// Save it in imposedJoints [q0,qd0,qdd0,q1,qd1,qdd1,...]
        virtual void loadImposedMotion
        (
	Field<label>& jointIndex,
	Field<scalar>& imposedJoints
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
