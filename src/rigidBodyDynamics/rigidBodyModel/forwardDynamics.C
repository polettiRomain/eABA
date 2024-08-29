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

#include "rigidBodyModel.H"
#include "rigidBodyModelState.H"
#include "rigidBodyRestraint.H"
#include "rigidBodyimposedmotion.H"
#include "addToRunTimeSelectionTable.H"
#include "Time.H"

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::RBD::rigidBodyModel::applyRestraints
(
    scalarField& tau,
    Field<spatialVector>& fx,
    const rigidBodyModelState& state
) const
{
    if (restraints_.empty())
    {
        return;
    }

    forAll(restraints_, ri)
    {
        // Accumulate the restraint forces
        restraints_[ri].restrain(tau, fx, state);
    }
}

void Foam::RBD::rigidBodyModel::applyImposedMotion
(
    Field<label>& jointIndex,
    Field<scalar>& imposedJoints
) const
{
    if (ImposedMotion_.empty())
    {
        return;
    }

    forAll(ImposedMotion_, mi) // tested with only one entry in ImposedMotion dict so far
    {
    DebugInfo << "imposedmotion " << ImposedMotion_[mi].name();
    ImposedMotion_[mi].loadImposedMotion(jointIndex,imposedJoints);
    }
}

void Foam::RBD::rigidBodyModel::forwardDynamics
(
    rigidBodyModelState& state,
    const scalarField& tau,
    const Field<spatialVector>& fx,
    Field<label>&  indexImposedJoints,
    Field<scalar>& imposedJoints
) const
{
    const scalarField& q = state.q();
    const scalarField& qDot = state.qDot();
    scalarField& qDdot = state.qDdot();

    DebugInFunction
        << "q = " << q << nl
        << "qDot = " << qDot << nl
        << "tau = " << tau << endl;

    Info << "*************** fwd*****************" << endl;	
    Info << q << endl; 
    Info << qDot << endl; 
    Info << qDdot << endl; 


    // Joint state returned by jcalc
    joint::XSvc J;

    v_[0] = Zero;


    for (label i=1; i<nBodies(); i++)
    {
        // Select joint i
        const joint& jnt = joints()[i];

        // Fill the J object with J.X, J.S1, J.v, J.v.wx/y/z(), J.c. Value depends on the type of joint (i.e. different implementation of jcalc for diff joints)
	/*
        Rev joint - x
        J.X  = Xt(S_[0].l()*q[qIndex_]); Translational spatial transformation tensor for translation r (Xt takes a vector)
        J.S1 = S_[0] where S_[0] = (0,0,0,1,0,0), i.e. unit spatial vector in the direction of the motion
        J.v  = S_[0]*qDot[qIndex_]; rel velocity between body i and i-1 due to the joint
        J.c  = 0
	*/
	/*
        Prism joint - x
        J.X =Xrx(q[qIndex_]); Rotational spatial transformation tensor about the x-axis by omega radians (i.e. transformation matrix from frame i-1 to i: i_X_(i-1)). Xt Xrx a scalar
        J.S1 = S_[0] where S_[0] = (1,0,0,0,0,0), i.e. unit spatial vector in the direction of the motion
        J.v = 0: init
	J.v.wx() = this->v_[WX] = qDot[qIndex_];
        J.c = 0
	*/
        jnt.jcalc(J, q, qDot);
        

        S_[i] = J.S; // for joints that have more than 1 dof
        S1_[i] = J.S1;
        
	// XT_: Transform from the parent body frame to the joint frame
        //  Entry in dMdict: transform (1 0 0 0 1 0 0 0 1)(0 0 0);
        // Not used in flapping test case 
        Xlambda_[i] = J.X & XT_[i];


	//lambda_: List of indices of the parent of each body (i-1 in a pure serie linkage)
        const label lambdai = lambda_[i];
	// Transformation matrix applied to the forces because they are defined in the inertial frame 
        // X0_ stacks the Xlambda
        if (lambdai != 0)
        {
            X0_[i] = Xlambda_[i] & X0_[lambdai];
        }
        else
        {
            X0_[i] = Xlambda_[i];
        }
	
	// velocity of link i: qdot + i_X_(i-1) v_(i-1)
        // this means than J.v is defined in the frame of ref from the link
        v_[i] = (Xlambda_[i] & v_[lambdai]) + J.v;
	
	// velocity product term
        c_[i] = J.c + (v_[i] ^ J.v);
	// isolated inertia, the rest is filled afterwards
        IA_[i] = I(i);

	// zero-acceleration force
        pA_[i] = v_[i] ^* (I(i) & v_[i]);

        // zero-acceleration force: External forces put in link i frame
        if (fx.size())
        {
            pA_[i] -= *X0_[i] & fx[i];
        }

    }


    // Tip to base loop: articulated inertia and z.a. force 
    bool isImposed;
    for (label i=nBodies()-1; i>0; i--)
    {
        const joint& jnt = joints()[i];
        const label qi = jnt.qIndex();

        if (jnt.nDoF() == 1)
        {
            U1_[i] = IA_[i] & S1_[i];
            Dinv_[i].xx() = 1/(S1_[i] && U1_[i]);
            u_[i].x() = tau[qi] - (S1_[i] && pA_[i]);

            const label lambdai = lambda_[i];
	    

            if (lambdai != 0)
            {

		// Check if active joint
		isImposed = false;
		forAll(indexImposedJoints, j)
		   {
			if (indexImposedJoints[j] == qi)
			{
				//Info << "Imposed Joint detected" << endl;
				//Info<< qi<< indexImposedJoints[j] << endl;
				isImposed = true;
				break;
			}
		    }
                // Free joint	
                if (!isImposed)
		{
			//Info << "Imposed Joint NOT detected" << endl;
		        const spatialTensor Ia
		        (
		            IA_[i] - (U1_[i]*(Dinv_[i].xx()*U1_[i]))
		        );


		        const spatialVector pa
		        (
		            pA_[i] + (Ia & c_[i]) + U1_[i]*(Dinv_[i].xx()*u_[i].x())
		        );


		        IA_[lambdai] +=
		            spatialTensor(Xlambda_[i].T())
		          & Ia
		          & spatialTensor(Xlambda_[i]);

		        pA_[lambdai] += Xlambda_[i].T() & pa;
		}
		// Active joint
		else
		{		
			IA_[lambdai] += 
			          spatialTensor(Xlambda_[i].T())
				  & IA_[i]
				  & spatialTensor(Xlambda_[i]);

			pA_[lambdai] += spatialTensor(Xlambda_[i].T()) & (pA_[i] + (IA_[i] & c_[i]) + U1_[i] * qDdot[qi]);	
		}
            }
        }
        else
        {
            U_[i] = IA_[i] & S_[i];
            Dinv_[i] = (S_[i].T() & U_[i]).inv();

            u_[i] = tau.block<vector>(qi) - (S_[i].T() & pA_[i]);

            const label lambdai = lambda_[i];

            if (lambdai != 0)
            {
                spatialTensor Ia
                (
                    IA_[i]
                  - (U_[i] & Dinv_[i] & U_[i].T())
                );

                spatialVector pa
                (
                    pA_[i]
                  + (Ia & c_[i])
                  + (U_[i] & Dinv_[i] & u_[i])
                );

                IA_[lambdai] +=
                    spatialTensor(Xlambda_[i].T())
                  & Ia
                  & spatialTensor(Xlambda_[i]);

                pA_[lambdai] += Xlambda_[i].T() & pa;
            }
        }
    }

    a_[0] = spatialVector(Zero, -g_);
    
    // Base to tip loop: acceleration
    for (label i=1; i<nBodies(); i++)
    {
        const joint& jnt = joints()[i];
        const label qi = jnt.qIndex();

        a_[i] = (Xlambda_[i] & a_[lambda_[i]]) + c_[i];

        if (jnt.nDoF() == 1)
        {
	    
	    // Check if active joint
	    isImposed = false;
	    forAll(indexImposedJoints, j)
	    {
		if (indexImposedJoints[j] == qi)
		{
			//Info << "2. Imposed Joint detected" << endl;
			//Info<< qi<< indexImposedJoints[j] <<endl;
		        isImposed = true;
			break;
		}
	    }	
	    //  if not compute qdDot
	    if (!isImposed) 
 	    { 
		    //Info << "2. Imposed Joint NOT detected" << endl;
		    qDdot[qi] = Dinv_[i].xx()*(u_[i].x() - (U1_[i] && a_[i]));
	    }
            // Compute a for active and free joints
	    a_[i] += S1_[i]*qDdot[qi];

        }
        else
        {
            vector qDdoti(Dinv_[i] & (u_[i] - (U_[i].T() & a_[i])));

            // Need to add mutable "block<vector>" to Field
            qDdot[qi] = qDdoti.x();
            qDdot[qi+1] = qDdoti.y();
            qDdot[qi+2] = qDdoti.z();

            a_[i] += (S_[i] & qDdoti);
        }
    }

    DebugInfo
        << "qDdot = " << qDdot << nl
        << "a = " << a_ << endl;


}


void Foam::RBD::rigidBodyModel::forwardDynamicsCorrection
(
    const rigidBodyModelState& state
) const
{

    
    DebugInFunction << endl;

    const scalarField& q = state.q();
    const scalarField& qDot = state.qDot();
    const scalarField& qDdot = state.qDdot();

    // Joint state returned by jcalc
    joint::XSvc J;

    v_[0] = Zero;
    a_[0] = spatialVector(Zero, -g_);

    for (label i=1; i<nBodies(); i++)
    {
        const joint& jnt = joints()[i];
        const label qi = jnt.qIndex();

        jnt.jcalc(J, q, qDot);

        S_[i] = J.S;
        S1_[i] = J.S1;

        Xlambda_[i] = J.X & XT_[i];

        const label lambdai = lambda_[i];
	
        if (lambdai != 0)
        {
            X0_[i] = Xlambda_[i] & X0_[lambdai];
        }
        else
        {
            X0_[i] = Xlambda_[i];
        }

        v_[i] = (Xlambda_[i] & v_[lambdai]) + J.v;
        c_[i] = J.c + (v_[i] ^ J.v);
        a_[i] = (Xlambda_[i] & a_[lambdai]) + c_[i];
        

        if (jnt.nDoF() == 1)
        {
            a_[i] += S1_[i]*qDdot[qi];
        }
        else
        {
            a_[i] += S_[i] & qDdot.block<vector>(qi);
        }
    }
    DebugInfo<< "a = " << a_ << endl;
}




// ************************************************************************* //
