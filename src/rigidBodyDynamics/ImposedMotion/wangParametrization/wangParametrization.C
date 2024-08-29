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

#include "wangParametrization.H"
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
    defineTypeNameAndDebug(wangParametrization, 0);

    addToRunTimeSelectionTable
    (
        imposedmotion,
        wangParametrization,
        dictionary
    );
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::RBD::ImposedMotion::wangParametrization::wangParametrization
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

Foam::RBD::ImposedMotion::wangParametrization::~wangParametrization()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::RBD::ImposedMotion::wangParametrization::loadImposedMotion
(
   Field<label>& jointIndex,
   Field<scalar>& imposedJoints
) const
{

	// time
	scalar ti;
	ti = model_.time().value();

	if (motionType_ == "midStart")  // wing motion constant during sim
	{


		if (jointList_.size() == 2) // by definition the Wang parametrization has only 2 dof
		{

			Field<scalar> varList(jointList_.size()*3); // contain the pos, vel, acc of the 2 joints
			
			// Berman, G. J., & Wang, Z. J. (2007). Energy-minimizing kinematics in hovering insect flight
			// phi - flapping motion
			scalar q_phi = -Aphi_/asin(Kphi_) * asin(Kphi_*sin(2*M_PI*fphi_*ti)); 
			scalar qDot_phi =  -Aphi_/asin(Kphi_) * 2*M_PI*fphi_*Kphi_ * cos(2*M_PI*fphi_*ti) / (sqrt(1-(Kphi_*Kphi_)*(sqr(sin(2*M_PI*fphi_*ti)))));
			scalar qDdot_phi = -(Aphi_ * Kphi_ * 4*M_PI*M_PI*fphi_*fphi_ * sin(2*M_PI*fphi_*ti)* (sqr(Kphi_)*sqr(sin(2*M_PI*fphi_*ti)) + sqr(Kphi_)*sqr(cos(2*M_PI*fphi_*ti)) - 1)) / (asin(Kphi_) * pow((1 - sqr(Kphi_) * sqr(sin(2*M_PI*fphi_*ti))),1.5) );
			varList[0]   = q_phi; varList[1] = qDot_phi; varList[2] = qDdot_phi; 

			// alpha - pitching motion
			scalar q_alpha = -Aalpha_/tanh(Kalpha_)*tanh(Kalpha_*cos(2*M_PI*falpha_*ti)) + Aalpha_;
			scalar qDot_alpha = 2*M_PI*falpha_*(Aalpha_/tanh(Kalpha_)) *Kalpha_* sin(2*M_PI*falpha_*ti)/sqr(cosh(Kalpha_*sin(2*M_PI*falpha_*ti)));
			scalar qDdot_alpha = 4 * M_PI*M_PI * falpha_*falpha_ *(Aalpha_/tanh(Kalpha_))*Kalpha_* 1/sqr(cosh(Kalpha_*cos(2*M_PI*falpha_*ti))) * (2*Kalpha_*sqr(sin(2*M_PI*falpha_*ti)) * tanh(Kalpha_*cos(2*M_PI*falpha_*ti)) + cos(2*M_PI*falpha_*ti));
			varList[3] = q_alpha; varList[4] = qDot_alpha; varList[5] = qDdot_alpha; 
			
			jointIndex    = jointList_;
			imposedJoints = varList;

			Info << "PitchingAngle(" << q_alpha << endl;
			Info << "FlappingAngle(" << q_phi << endl;

		}
		else
		{
			Info << "Execution continue but no imposed motion was loaded" << endl;
			Info << "Please provide a list containing the index of the 2 joints in the ImposedMotion sub dictionnary " << endl;
		}

	}
	// 3 dof body: x,z, theta
	// 3 control var: Aphi, beta, Aoff
	// Parametrization: 'PD' for each control var: 
	// Aphi = SCALED(tanh(Kp*(z-zref) +  Kd*(z-zref)))
	// Beta = SCALED(tanh(Kp*(x-xref) +  Kd*(x-xref))
	// Aoff = SCALED(tanh(Kp*(th-thref) +  Kd*(th-thref)))
	// where SCALED is mapping from [-1,1] to [control_var_min, control_var_max]
	if (motionType_ == "midStartControl_xzth_PD") 
	{
		
		// Done only once to recover the correct thold if the simulation is not start from 0
		if (flag_init_ == 0)
		{
		 	scalar tinit = ti;
			thold_ = 0;
			while ((tinit-1/fphi_) >= 0)
			{
				tinit  -= 1/fphi_;
				thold_ += 1/fphi_;
			}
			Info << "thold start is: "<< thold_ << endl;
			flag_init_ = 1;
		}


		if (jointList_.size() == 3) // 3 dof: beta, phi, alpha
		{

			Field<scalar> varList(jointList_.size()*3);
			

			// End of cycle control 
			if ((ti - thold_) >= (1/fphi_) )
			{	
			
				Info << "End of cycle --- " << endl;

				thold_ = thold_ + 1/fphi_;

				label sphereBodyID = 3;

				vector xyz = model_.X0(sphereBodyID).r();
				vector xyzd = model_.v(sphereBodyID, Zero).l();
				scalar theta = model_.X0(sphereBodyID).E()[2];
				scalar thetad = model_.v(sphereBodyID, Zero).w()[1];

				// Vertical motion
				scalar e_z = xyz[2] - href_;
				scalar ed_z = xyzd[2] - hdref_;
				scalar ascaled = tanh(Kpz_*e_z + Kdz_*ed_z + biasz_);
				Aphi_PID = - ((ascaled + 1) / 2 * (Aphi_max_- Aphi_min_) + Aphi_min_);
				Info << "CONTROL --- ascaled " << ascaled << endl;
				Info << "CONTROL --- Aphi_PID " << Aphi_PID << endl;

				// Longitudinal motion
				scalar e_x = xyz[0] - xref_;
				scalar ed_x = xyzd[0] - xdref_;
				scalar betascaled = tanh(Kpx_*e_x + Kdx_*ed_x);
				beta_PID = - ((betascaled + 1) / 2 * (beta_max_- beta_min_) + beta_min_);
		                Info << "CONTROL --- bscaled " << betascaled << endl;
		                Info << "CONTROL --- beta_PID " << beta_PID << endl;
				
				// Pitch motion
				Info << "theta" << theta << endl;
		                Info << "theta" << asin(theta) << endl;
					
		                scalar e_theta = asin(theta) - thetaref_;
		                scalar ed_theta = thetad - thetadref_;

				Info << e_theta << endl;
		                Info << ed_theta << endl;

				scalar thetascaled = tanh(Kptheta_*e_theta + Kdtheta_*ed_theta);
				A_phi_off =  ((thetascaled + 1) / 2 * (theta_max_-theta_min_) + theta_min_);
		                Info << "CONTROL --- thetascaled " << thetascaled << endl;
		                Info << "CONTROL --- A_phi_off " << A_phi_off << endl;

			}
			scalar a;
			double one = 1;
			double tstar = ti*fphi_;
			if (fmod(tstar, one) < 0.5) // downstroke
			{
				a = (Aphi_ + A_phi_off)/Aphi_;
		        }
		    	else // upstroke
		    	{
				a = (Aphi_ -  A_phi_off)/Aphi_;
		        }
			// beta
			scalar q_beta = beta_PID;
			scalar qDot_beta =  0;
			scalar qDdot_beta = 0;
			varList[0]   = q_beta; varList[1] = qDot_beta; varList[2] = qDdot_beta; 

			
			// phi
			scalar q_phi = -Aphi_PID/asin(Kphi_) * asin(Kphi_*sin(2*M_PI*fphi_*ti))*a; 
			scalar qDot_phi =  -Aphi_PID/asin(Kphi_) * 2*M_PI*fphi_*Kphi_ * cos(2*M_PI*fphi_*ti) / (sqrt(1-(Kphi_*Kphi_)*(sqr(sin(2*M_PI*fphi_*ti)))))*a;
			scalar qDdot_phi = -(Aphi_PID * Kphi_ * 4*M_PI*M_PI*fphi_*fphi_ * sin(2*M_PI*fphi_*ti)* (sqr(Kphi_)*sqr(sin(2*M_PI*fphi_*ti)) + sqr(Kphi_)*sqr(cos(2*M_PI*fphi_*ti)) - 1)) / (asin(Kphi_) * pow((1 - sqr(Kphi_) * sqr(sin(2*M_PI*fphi_*ti))),1.5) )*a;
			varList[3]   = q_phi; varList[4] = qDot_phi; varList[5] = qDdot_phi; 

			// alpha
			scalar q_alpha = -Aalpha_/tanh(Kalpha_)*tanh(Kalpha_*cos(2*M_PI*falpha_*ti)) + Aalpha_;
			scalar qDot_alpha = 2*M_PI*falpha_*(Aalpha_/tanh(Kalpha_)) *Kalpha_* sin(2*M_PI*falpha_*ti)/sqr(cosh(Kalpha_*sin(2*M_PI*falpha_*ti)));
			scalar qDdot_alpha = 4 * M_PI*M_PI * falpha_*falpha_ *(Aalpha_/tanh(Kalpha_))*Kalpha_* 1/sqr(cosh(Kalpha_*cos(2*M_PI*falpha_*ti))) * (2*Kalpha_*sqr(sin(2*M_PI*falpha_*ti)) * tanh(Kalpha_*cos(2*M_PI*falpha_*ti)) + cos(2*M_PI*falpha_*ti));
			varList[6] = q_alpha; varList[7] = qDot_alpha; varList[8] = qDdot_alpha; 
			
			jointIndex    = jointList_;
			imposedJoints = varList;

			Info << "PitchingAngle(" << q_alpha << endl;
			Info << "FlappingAngle(" << q_phi << endl;

		}
		else
		{
			Info << "Execution continue but no imposed motion was loaded" << endl;
			Info << "Please provide a list containing the index of the 2 joints in the ImposedMotion sub dictionnary " << endl;
		}

	}
	// 3 dof body: x,z, theta
	// 3 control var: Aphi, beta, Aoff
	// Parametrization: 'linear' = PD for each control var, function of x,z,theta and xd,zd, thetad
	if (motionType_ == "midStartControl_xzth_linearPD") 
	{
		
		// Done only once to recover the correct thold
		if (flag_init_ == 0)
		{
		 	scalar tinit = ti;
			thold_ = 0;
			while ((tinit-1/fphi_) >= 0)
			{
				tinit  -= 1/fphi_;
				thold_ += 1/fphi_;
			}
			Info << "thold start is: "<< thold_ << endl;
			flag_init_ = 1;
		}


		if (jointList_.size() == 3) // 3 dof: beta, phi, alpha
		{

			Field<scalar> varList(jointList_.size()*3);		

			// End of cycle control
			if ((ti - thold_) >= (1/fphi_) )
			{	
			
				Info << "End of cycle --- " << endl;

				thold_ = thold_ + 1/fphi_;

				label sphereBodyID = 3;

				vector xyz = model_.X0(sphereBodyID).r();
				vector xyzd = model_.v(sphereBodyID, Zero).l();
				scalar msintheta = model_.X0(sphereBodyID).E()[2]; 
				scalar thetad = model_.v(sphereBodyID, Zero).w()[1];

				// Error definition
				scalar e_z  = -(xyz[2] - href_);    // "-" to be consistent with python environment where z is pointing down
				scalar ed_z = -(xyzd[2] - hdref_);  // "-" to be consistent with python environment where z is pointing down
				scalar e_x  = xyz[0] - xref_;					
				scalar ed_x = xyzd[0] - xdref_;
				scalar e_theta  = -asin(msintheta) - thetaref_;
		                scalar ed_theta = thetad - thetadref_;  
				
				// Flapping amplitude
				scalar ascaled    = tanh(Kpz_*e_z + Kdz_*ed_z + Kpxz_*e_x + Kdxz_*ed_x + Kpthz_*e_theta + Kdthz_*ed_theta + biasz_);
				Aphi_PID       = - ((ascaled + 1) / 2 * (Aphi_max_- Aphi_min_) + Aphi_min_);	 // - due to the parametrization below
				// Stroke angle
				scalar betascaled = tanh(Kpx_*e_x + Kdx_*ed_x + Kpzx_*e_z + Kdzx_*ed_z + Kpthx_*e_theta + Kdthx_*ed_theta + biasx_);
				beta_PID          = - ((betascaled + 1) / 2 * (beta_max_- beta_min_) + beta_min_); // - because beta is defined in the other direction in python
				// Flapping amplitude offset
				scalar thetascaled = tanh(Kptheta_*e_theta + Kdtheta_*ed_theta +  Kpzth_*e_z + Kdzth_*ed_z  + Kpxth_*e_x + Kdxth_*ed_x +  biasth_);
				A_phi_off =  ((thetascaled + 1) / 2 * (theta_max_-theta_min_) + theta_min_);
				
				
				Info << "CONTROL --- Aphi_PID " << Aphi_PID << endl;
		                Info << "CONTROL --- beta_PID " << beta_PID << endl;
		                Info << "CONTROL --- A_phi_off_PID " << A_phi_off << endl;
		                
		                Info << e_theta << endl;
		                Info << ed_theta << endl;


			}
			scalar a;
			double one = 1;
			double tstar = ti*fphi_;
			if (fmod(tstar, one) < 0.5) // downstroke
			{
				a = (Aphi_ + A_phi_off)/Aphi_;
		        }
		    	else // upstroke
		    	{
				a = (Aphi_ -  A_phi_off)/Aphi_;
		        }
			// beta
			scalar q_beta = beta_PID;
			scalar qDot_beta =  0;
			scalar qDdot_beta = 0;
			varList[0]   = q_beta; varList[1] = qDot_beta; varList[2] = qDdot_beta; 

			
			// phi
			scalar q_phi = -Aphi_PID/asin(Kphi_) * asin(Kphi_*sin(2*M_PI*fphi_*ti))*a; 
			scalar qDot_phi =  -Aphi_PID/asin(Kphi_) * 2*M_PI*fphi_*Kphi_ * cos(2*M_PI*fphi_*ti) / (sqrt(1-(Kphi_*Kphi_)*(sqr(sin(2*M_PI*fphi_*ti)))))*a;
			scalar qDdot_phi = -(Aphi_PID * Kphi_ * 4*M_PI*M_PI*fphi_*fphi_ * sin(2*M_PI*fphi_*ti)* (sqr(Kphi_)*sqr(sin(2*M_PI*fphi_*ti)) + sqr(Kphi_)*sqr(cos(2*M_PI*fphi_*ti)) - 1)) / (asin(Kphi_) * pow((1 - sqr(Kphi_) * sqr(sin(2*M_PI*fphi_*ti))),1.5) )*a;
			varList[3]   = q_phi; varList[4] = qDot_phi; varList[5] = qDdot_phi; 

			// alpha
			scalar q_alpha = -Aalpha_/tanh(Kalpha_)*tanh(Kalpha_*cos(2*M_PI*falpha_*ti)) + Aalpha_;
			scalar qDot_alpha = 2*M_PI*falpha_*(Aalpha_/tanh(Kalpha_)) *Kalpha_* sin(2*M_PI*falpha_*ti)/sqr(cosh(Kalpha_*sin(2*M_PI*falpha_*ti)));
			scalar qDdot_alpha = 4 * M_PI*M_PI * falpha_*falpha_ *(Aalpha_/tanh(Kalpha_))*Kalpha_* 1/sqr(cosh(Kalpha_*cos(2*M_PI*falpha_*ti))) * (2*Kalpha_*sqr(sin(2*M_PI*falpha_*ti)) * tanh(Kalpha_*cos(2*M_PI*falpha_*ti)) + cos(2*M_PI*falpha_*ti));
			varList[6] = q_alpha; varList[7] = qDot_alpha; varList[8] = qDdot_alpha; 
			
			jointIndex    = jointList_;
			imposedJoints = varList;

			Info << "PitchingAngle(" << q_alpha << endl;
			Info << "FlappingAngle(" << q_phi << endl;

		}
		else
		{
			Info << "Execution continue but no imposed motion was loaded" << endl;
			Info << "Please provide a list containing the index of the 2 joints in the ImposedMotion sub dictionnary " << endl;
		}

	}
	// 1 controlled dof body: z
	// 1 control var: Aphi
	// Parametrization: PD function of z and zd
	else if (motionType_ == "midStartControl_z_PD") 
	{
		if (jointList_.size() == 2) 
		{
			Field<scalar> varList(jointList_.size()*3);

			if ((ti - thold_) >= (1/fphi_) )
			{	
				Info << "End of cycle --- " << endl;

				thold_ = thold_ + 1/fphi_;

				vector xyz = model_.X0(bodyID_).r();
				vector xyzd = model_.v(bodyID_, Zero).l();

				// Vertical motion
				scalar e_z = xyz[2] - href_;
				scalar ed_z = xyzd[2] - hdref_;
				scalar ascaled = tanh(Kpz_*e_z + Kdz_*ed_z + biasz_);
				Aphi_PID = - ((ascaled + 1) / 2 * (Aphi_max_- Aphi_min_) + Aphi_min_);
				Info << "CONTROL --- ascaled " << ascaled << endl;
				Info << "CONTROL --- Aphi_PID " << Aphi_PID << endl;

			}
			
			// phi
			scalar q_phi = -Aphi_PID/asin(Kphi_) * asin(Kphi_*sin(2*M_PI*fphi_*ti)); 
			scalar qDot_phi =  -Aphi_PID/asin(Kphi_) * 2*M_PI*fphi_*Kphi_ * cos(2*M_PI*fphi_*ti) / (sqrt(1-(Kphi_*Kphi_)*(sqr(sin(2*M_PI*fphi_*ti)))));
			scalar qDdot_phi = -(Aphi_PID * Kphi_ * 4*M_PI*M_PI*fphi_*fphi_ * sin(2*M_PI*fphi_*ti)* (sqr(Kphi_)*sqr(sin(2*M_PI*fphi_*ti)) + sqr(Kphi_)*sqr(cos(2*M_PI*fphi_*ti)) - 1)) / (asin(Kphi_) * pow((1 - sqr(Kphi_) * sqr(sin(2*M_PI*fphi_*ti))),1.5) );
			varList[0]   = q_phi; varList[1] = qDot_phi; varList[2] = qDdot_phi; 

			// alpha
			scalar q_alpha = -Aalpha_/tanh(Kalpha_)*tanh(Kalpha_*cos(2*M_PI*falpha_*ti)) + Aalpha_;
			scalar qDot_alpha = 2*M_PI*falpha_*(Aalpha_/tanh(Kalpha_)) *Kalpha_* sin(2*M_PI*falpha_*ti)/sqr(cosh(Kalpha_*sin(2*M_PI*falpha_*ti)));
			scalar qDdot_alpha = 4 * M_PI*M_PI * falpha_*falpha_ *(Aalpha_/tanh(Kalpha_))*Kalpha_* 1/sqr(cosh(Kalpha_*cos(2*M_PI*falpha_*ti))) * (2*Kalpha_*sqr(sin(2*M_PI*falpha_*ti)) * tanh(Kalpha_*cos(2*M_PI*falpha_*ti)) + cos(2*M_PI*falpha_*ti));
			varList[3] = q_alpha; varList[4] = qDot_alpha; varList[5] = qDdot_alpha; 
			
			jointIndex    = jointList_;
			imposedJoints = varList;

			Info << "PitchingAngle(" << q_alpha << endl;
			Info << "FlappingAngle(" << q_phi << endl;

		}
		else
		{
			Info << "Execution continue but no imposed motion was loaded" << endl;
			Info << "Please provide a list containing the index of the 2 joints in the ImposedMotion sub dictionnary " << endl;
		}

	}
	// no control
	// drone with 2 wings
	else if (motionType_ == "midStart_2wings") 
	{
	
		if (jointList_.size() == 4) 
		{

			Field<scalar> varList(jointList_.size()*3);
			
			// phi
			scalar q_phi = -Aphi_/asin(Kphi_) * asin(Kphi_*sin(2*M_PI*fphi_*ti)); 
			scalar qDot_phi =  -Aphi_/asin(Kphi_) * 2*M_PI*fphi_*Kphi_ * cos(2*M_PI*fphi_*ti) / (sqrt(1-(Kphi_*Kphi_)*(sqr(sin(2*M_PI*fphi_*ti)))));
			scalar qDdot_phi = -(Aphi_ * Kphi_ * 4*M_PI*M_PI*fphi_*fphi_ * sin(2*M_PI*fphi_*ti)* (sqr(Kphi_)*sqr(sin(2*M_PI*fphi_*ti)) + sqr(Kphi_)*sqr(cos(2*M_PI*fphi_*ti)) - 1)) / (asin(Kphi_) * pow((1 - sqr(Kphi_) * sqr(sin(2*M_PI*fphi_*ti))),1.5) );
			varList[0]   = q_phi; varList[1] = qDot_phi; varList[2] = qDdot_phi; 
			
			// alpha
			scalar q_alpha = -Aalpha_/tanh(Kalpha_)*tanh(Kalpha_*cos(2*M_PI*falpha_*ti)) + Aalpha_;
			scalar qDot_alpha = 2*M_PI*falpha_*(Aalpha_/tanh(Kalpha_)) *Kalpha_* sin(2*M_PI*falpha_*ti)/sqr(cosh(Kalpha_*sin(2*M_PI*falpha_*ti)));
			scalar qDdot_alpha = 4 * M_PI*M_PI * falpha_*falpha_ *(Aalpha_/tanh(Kalpha_))*Kalpha_* 1/sqr(cosh(Kalpha_*cos(2*M_PI*falpha_*ti))) * (2*Kalpha_*sqr(sin(2*M_PI*falpha_*ti)) * tanh(Kalpha_*cos(2*M_PI*falpha_*ti)) + cos(2*M_PI*falpha_*ti));
			varList[3] = q_alpha; varList[4] = qDot_alpha; varList[5] = qDdot_alpha; 

			varList[6]   = -q_phi; varList[7] = -qDot_phi; varList[8] = -qDdot_phi; 
			varList[9] = q_alpha; varList[10] = qDot_alpha; varList[11] = qDdot_alpha; 

			jointIndex    = jointList_;
			imposedJoints = varList;

			Info << "PitchingAngle(" << q_alpha << endl;
			Info << "FlappingAngle(" << q_phi << endl;

		}
		else
		{
			Info << "Execution continue but no imposed motion was loaded" << endl;
			Info << "Please provide a list containing the index of the 2 joints in the ImposedMotion sub dictionnary " << endl;
		}

	}
}


bool Foam::RBD::ImposedMotion::wangParametrization::read
(
    const dictionary& dict
)
{
    imposedmotion::read(dict);
    
    
    // joint index and type of joint's motion
    coeffs_.readEntry("jointList", jointList_);
    coeffs_.readEntry("motionType", motionType_);
    
    // Flapping motion
    coeffs_.readEntry("Aphi", Aphi_);
    coeffs_.readEntry("fphi", fphi_);
    coeffs_.readEntry("Kphi", Kphi_);
    coeffs_.readEntry("Aalpha", Aalpha_);
    coeffs_.readEntry("falpha", falpha_);
    coeffs_.readEntry("Kalpha", Kalpha_);
    
    // If controlled flight: load control parameters
    if ((motionType_ == "midStartControl_xzth_linearPD") || (motionType_ == "midStartControl_xzth_PD"))
    {

	    //coeffs_.readEntry("thold", thold_);
	    coeffs_.readEntry("Kpz", Kpz_);
	    coeffs_.readEntry("Kdz", Kdz_);
	    coeffs_.readEntry("Kiz", Kiz_);
	    
	    coeffs_.readEntry("Kpxz", Kpxz_);
	    coeffs_.readEntry("Kdxz", Kdxz_);
	    coeffs_.readEntry("Kpthz", Kpthz_);
	    coeffs_.readEntry("Kdthz", Kdthz_);
	    
	    coeffs_.readEntry("Kpx", Kpx_);
	    coeffs_.readEntry("Kdx", Kdx_);
	    coeffs_.readEntry("Kix", Kix_);
	    
	    coeffs_.readEntry("Kpzx", Kpzx_);
	    coeffs_.readEntry("Kdzx", Kdzx_);
	    coeffs_.readEntry("Kpthx", Kpthx_);
	    coeffs_.readEntry("Kdthx", Kdthx_);
	    
	    coeffs_.readEntry("Kptheta", Kptheta_);
	    coeffs_.readEntry("Kdtheta", Kdtheta_);
	    coeffs_.readEntry("Kitheta", Kitheta_);
	    
	    coeffs_.readEntry("Kpzth", Kpzth_);
	    coeffs_.readEntry("Kdzth", Kdzth_);
	    coeffs_.readEntry("Kpxth", Kpxth_);
	    coeffs_.readEntry("Kdxth", Kdxth_);
	    
	    
	    coeffs_.readEntry("biasz", biasz_);
	    coeffs_.readEntry("biasx", biasx_);
	    coeffs_.readEntry("biasth", biasth_);
	    
	    coeffs_.readEntry("Aphi_max", Aphi_max_);
	    coeffs_.readEntry("Aphi_min", Aphi_min_);
	    coeffs_.readEntry("beta_max", beta_max_);
	    coeffs_.readEntry("beta_min", beta_min_);
	    coeffs_.readEntry("theta_max", theta_max_);
	    coeffs_.readEntry("theta_min", theta_min_);
	    
	    coeffs_.readEntry("href", href_);
	    coeffs_.readEntry("hdref", hdref_);
	    coeffs_.readEntry("xref", xref_);
	    coeffs_.readEntry("xdref", xdref_);
	    coeffs_.readEntry("thetaref", thetaref_);
	    coeffs_.readEntry("thetadref", thetadref_);
    }
    else if (motionType_ == "midStartControl_z_PD")
    {

	    //coeffs_.readEntry("thold", thold_);
	    coeffs_.readEntry("Kpz", Kpz_);
	    coeffs_.readEntry("Kdz", Kdz_);
	    coeffs_.readEntry("Kiz", Kiz_);

	    coeffs_.readEntry("biasz", biasz_);
	    
	    coeffs_.readEntry("Aphi_max", Aphi_max_);
	    coeffs_.readEntry("Aphi_min", Aphi_min_);

	    coeffs_.readEntry("href", href_);
	    coeffs_.readEntry("hdref", hdref_);
    }
    return true;
}


void Foam::RBD::ImposedMotion::wangParametrization::write
(
    Ostream& os
) const
{
    imposedmotion::write(os);

    os.writeEntry("jointList", jointList_);
    os.writeEntry("motionType", motionType_);

    os.writeEntry("Aphi", Aphi_);
    os.writeEntry("fphi", fphi_);
    os.writeEntry("Kphi", Kphi_);
    os.writeEntry("Aalpha", Aalpha_);
    os.writeEntry("falpha", falpha_);
    os.writeEntry("Kalpha", Kalpha_);

    if ((motionType_ == "midStartControl_xzth_linearPD") || (motionType_ == "midStartControl_xzth_PD"))
    {
	    os.writeEntry("Kpz", Kpz_);
	    os.writeEntry("Kdz", Kdz_);
	    os.writeEntry("Kiz", Kiz_);
	    
	    os.writeEntry("Kpxz", Kpxz_);
	    os.writeEntry("Kdxz", Kdxz_);
	    os.writeEntry("Kpthz", Kpthz_);
	    os.writeEntry("Kdthz", Kdthz_);
	    
	    os.writeEntry("Kpx", Kpx_);
	    os.writeEntry("Kdx", Kdx_);
	    os.writeEntry("Kix", Kix_);
	    
	    os.writeEntry("Kpzx", Kpzx_);
	    os.writeEntry("Kdzx", Kdzx_);
	    os.writeEntry("Kpthx", Kpthx_);
	    os.writeEntry("Kdthx", Kdthx_);
	    
	    os.writeEntry("Kptheta", Kptheta_);
	    os.writeEntry("Kdtheta", Kdtheta_);
	    os.writeEntry("Kitheta", Kitheta_);
	    
	    os.writeEntry("Kpzth", Kpzth_);
	    os.writeEntry("Kdzth", Kdzth_);
	    os.writeEntry("Kpxth", Kpxth_);
	    os.writeEntry("Kdxth", Kdxth_);
	    
	    
	    os.writeEntry("biasz", biasz_);
	    os.writeEntry("biasx", biasx_);
	    os.writeEntry("biasth", biasth_);
	    
	    os.writeEntry("Aphi_max", Aphi_max_);
	    os.writeEntry("Aphi_min", Aphi_min_);
	    os.writeEntry("beta_max", beta_max_);
	    os.writeEntry("beta_min", beta_min_);
	    os.writeEntry("theta_max", theta_max_);
	    os.writeEntry("theta_min", theta_min_);
	    
	    os.writeEntry("href", href_);
	    os.writeEntry("hdref", hdref_);
	    os.writeEntry("xref", xref_);
	    os.writeEntry("xdref", xdref_);
	    os.writeEntry("thetaref", thetaref_);
	    os.writeEntry("thetadref", thetadref_);
	}
     else if ((motionType_ == "midStartControl_xzth_linearPD") || (motionType_ == "midStartControl_xzth_PD"))
    {
	    os.writeEntry("Kpz", Kpz_);
	    os.writeEntry("Kdz", Kdz_);
	    os.writeEntry("Kiz", Kiz_);
	   
	    os.writeEntry("biasz", biasz_);
	    
	    os.writeEntry("Aphi_max", Aphi_max_);
	    os.writeEntry("Aphi_min", Aphi_min_);

	    
	    os.writeEntry("href", href_);
	    os.writeEntry("hdref", hdref_);
	    }


}


// ************************************************************************* //
