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

#include "rigidBodyModel.H"
#include "masslessBody.H"
#include "compositeBody.H"
#include "jointBody.H"
#include "nullJoint.H"
#include "rigidBodyRestraint.H"
#include "rigidBodyimposedmotion.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace RBD
{
    defineTypeNameAndDebug(rigidBodyModel, 0);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::RBD::rigidBodyModel::initializeRootBody()
{
    bodies_.append(new masslessBody("root"));
    lambda_.append(0);
    bodyIDs_.insert("root", 0);
    joints_.append(new joints::null());
    XT_.append(spatialTransform());

    nDoF_ = 0;
    unitQuaternions_ = false;

    resizeState();
}


void Foam::RBD::rigidBodyModel::resizeState()
{
    Xlambda_.append(spatialTransform());
    X0_.append(spatialTransform());

    v_.append(Zero);
    a_.append(Zero);
    c_.append(Zero);

    IA_.append(spatialTensor::I);
    pA_.append(Zero);

    S_.append(Zero);
    S1_.append(Zero);
    U_.append(Zero);
    U1_.append(Zero);
    Dinv_.append(Zero);
    u_.append(Zero);
}


void Foam::RBD::rigidBodyModel::addRestraints
(
    const dictionary& dict
)
{
    if (dict.found("restraints"))
    {
        const dictionary& restraintDict = dict.subDict("restraints");

        label i = 0;

        restraints_.setSize(restraintDict.size());

        for (const entry& dEntry : restraintDict)
        {
            if (dEntry.isDict())
            {
                restraints_.set
                (
                    i++,
                    restraint::New
                    (
                        dEntry.keyword(),
                        dEntry.dict(),
                        *this
                    )
                );
            }
        }

        restraints_.setSize(i);
    }
}



void Foam::RBD::rigidBodyModel::addImposedMotion
(
    const dictionary& dict
)
{
    if (dict.found("ImposedMotion"))
    {
        const dictionary& imposedmotionDict = dict.subDict("ImposedMotion");
        //ImposedMotion_ = imposedmotion::New(imposedmotionDict,*this);

        label i = 0;

        ImposedMotion_.setSize(imposedmotionDict.size());

        for (const entry& dEntry : imposedmotionDict)
        {
            if (dEntry.isDict())
            {
                ImposedMotion_.set
                (
                    i++,
                    imposedmotion::New
                    (
                        dEntry.keyword(),
                        dEntry.dict(),
                        *this
                    )
                );
            }
        }

        ImposedMotion_.setSize(i);
        
    }
}



// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

Foam::label Foam::RBD::rigidBodyModel::join_
(
    const label parentID,
    const spatialTransform& XT,
    autoPtr<joint> jointPtr,
    autoPtr<rigidBody> bodyPtr
)
{
    // Append the body
    const rigidBody& body = bodyPtr();
    bodies_.append(bodyPtr);
    const label bodyID = nBodies()-1;
    bodyIDs_.insert(body.name(), bodyID);

    // If the parentID refers to a merged body find the parent into which it has
    // been merged and set lambda and XT accordingly
    if (merged(parentID))
    {
        const subBody& sBody = mergedBody(parentID);
        lambda_.append(sBody.masterID());
        XT_.append(XT & sBody.masterXT());
    }
    else
    {
        lambda_.append(parentID);
        XT_.append(XT);
    }

    // Append the joint
    const joint& prevJoint = joints_[joints_.size() - 1];
    joints_.append(jointPtr);
    joint& curJoint = joints_[joints_.size() - 1];
    curJoint.index() = joints_.size() - 1;
    curJoint.qIndex() = prevJoint.qIndex() + prevJoint.nDoF();

    // Increment the degrees of freedom
    nDoF_ += curJoint.nDoF();
    unitQuaternions_ = unitQuaternions_ || curJoint.unitQuaternion();

    resizeState();

    return bodyID;
}


// * * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * * //

Foam::RBD::rigidBodyModel::rigidBodyModel(const Time& time)
:
    time_(time),
    g_(Zero)
{
    initializeRootBody();
}


Foam::RBD::rigidBodyModel::rigidBodyModel
(
    const Time& time,
    const dictionary& dict
)
:
    time_(time),
    g_(Zero)
{
    initializeRootBody();

    const dictionary& bodiesDict = dict.subDict("bodies");
    Info << bodiesDict << endl;

    for (const entry& dEntry : bodiesDict)
    {
        const keyType& key = dEntry.keyword();
        const dictionary& bodyDict = dEntry.dict();

        Info << key << endl;
        Info << bodyDict << endl;

        if (bodyDict.found("mergeWith"))
        {
            merge
            (
                bodyID(bodyDict.get<word>("mergeWith")),
                bodyDict.get<spatialTransform>("transform"),
                rigidBody::New(key, bodyDict)
            );
        }
        else
        {
            Info << "---" << bodyDict.subDict("joint") << endl;
            join
            (
                bodyID(bodyDict.get<word>("parent")),
                bodyDict.get<spatialTransform>("transform"),
                joint::New(bodyDict.subDict("joint")),
                rigidBody::New(key, bodyDict)
            );
        }
    }

    // Read the restraints and any other re-readable settings.
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::RBD::rigidBodyModel::~rigidBodyModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::label Foam::RBD::rigidBodyModel::join
(
    const label parentID,
    const spatialTransform& XT,
    autoPtr<joint> jointPtr,
    autoPtr<rigidBody> bodyPtr
)
{
    if (isA<joints::composite>(jointPtr()))
    {
        return join
        (
            parentID,
            XT,
            autoPtr<joints::composite>
            (
                dynamic_cast<joints::composite*>(jointPtr.ptr())
            ),
            bodyPtr
        );
    }
    else
    {
        return join_
        (
            parentID,
            XT,
            jointPtr,
            bodyPtr
        );
    }
}


Foam::label Foam::RBD::rigidBodyModel::join
(
    const label parentID,
    const spatialTransform& XT,
    autoPtr<joints::composite> cJointPtr,
    autoPtr<rigidBody> bodyPtr
)
{
    label parent = parentID;
    joints::composite& cJoint = cJointPtr();

    // For all but the final joint in the set add a jointBody with the
    // joint and transform
    for (label j=0; j<cJoint.size()-1; j++)
    {
        parent = join_
        (
            parent,
            j == 0 ? XT : spatialTransform(),
            cJoint[j].clone(),
            autoPtr<rigidBody>(new jointBody)
        );
    }

    // For the final joint in the set add the real body
    parent = join_
    (
        parent,
        cJoint.size() == 1 ? XT : spatialTransform(),
        autoPtr<joint>(cJointPtr.ptr()),
        bodyPtr
    );

    // Set the properties of the last joint in the list to those set
    // by rigidBodyModel
    cJoint.setLastJoint();

    return parent;
}


void Foam::RBD::rigidBodyModel::makeComposite(const label bodyID)
{
    if (!isA<compositeBody>(bodies_[bodyID]))
    {
        // Retrieve the un-merged body
        autoPtr<rigidBody> bodyPtr = bodies_.release(bodyID);

        // Insert the compositeBody containing the original body
        bodies_.set
        (
            bodyID,
            new compositeBody(bodyPtr)
        );
    }
}


Foam::label Foam::RBD::rigidBodyModel::merge
(
    const label parentID,
    const spatialTransform& XT,
    autoPtr<rigidBody> bodyPtr
)
{
    autoPtr<subBody> sBodyPtr;

    // If the parentID refers to a merged body find the parent into which it has
    // been merged and merge this on into the same parent with the appropriate
    // transform
    if (merged(parentID))
    {
        const subBody& sBody = mergedBody(parentID);

        makeComposite(sBody.masterID());

        sBodyPtr.reset
        (
            new subBody
            (
                bodyPtr,
                bodies_[sBody.masterID()].name(),
                sBody.masterID(),
                XT & sBody.masterXT()
            )
        );
    }
    else
    {
        makeComposite(parentID);

        sBodyPtr.reset
        (
            new subBody
            (
                bodyPtr,
                bodies_[parentID].name(),
                parentID,
                XT
            )
        );
    }

    const subBody& sBody = sBodyPtr();
    mergedBodies_.append(sBodyPtr);

    // Merge the sub-body with the parent
    bodies_[sBody.masterID()].merge(sBody);

    const label sBodyID = mergedBodyID(mergedBodies_.size() - 1);
    bodyIDs_.insert(sBody.name(), sBodyID);

    return sBodyID;
}


Foam::spatialTransform Foam::RBD::rigidBodyModel::X0
(
    const label bodyId
) const
{
    if (merged(bodyId))
    {
        const subBody& mBody = mergedBody(bodyId);
        return mBody.masterXT() & X0_[mBody.masterID()];
    }

    return X0_[bodyId];
}


void Foam::RBD::rigidBodyModel::write(Ostream& os) const
{
    os.beginBlock("bodies");

    // Write the moving bodies
    for (label i=1; i<nBodies(); i++)
    {
        // Do not write joint-bodies created automatically to support elements
        // of composite joints
        if (!isType<jointBody>(bodies_[i]))
        {
            os.beginBlock(bodies_[i].name());

            bodies_[i].write(os);

            os.writeEntry("parent", bodies_[lambda_[i]].name());
            os.writeEntry("transform", XT_[i]);

            os  << indent << "joint" << nl
                << joints_[i] << endl;

            os.endBlock();
        }
    }

    // Write the bodies merged into the parent bodies for efficiency
    forAll(mergedBodies_, i)
    {
        os.beginBlock(mergedBodies_[i].name());

        mergedBodies_[i].body().write(os);

        os.writeEntry("transform", mergedBodies_[i].masterXT());
        os.writeEntry("mergeWith", mergedBodies_[i].masterName());

        os.endBlock();
    }

    os.endBlock();


    if (!restraints_.empty())
    {
        os.beginBlock("restraints");

        forAll(restraints_, ri)
        {
            // const word& restraintType(restraints_[ri].type());

            os.beginBlock(restraints_[ri].name());

            restraints_[ri].write(os);

            os.endBlock();
        }

        os.endBlock();
    }
    
    // Similarly to restraints, load the ImposedMotion library
    if (!ImposedMotion_.empty())
    {
        os.beginBlock("ImposedMotion");

        forAll(ImposedMotion_, mi)
        {

            os.beginBlock(ImposedMotion_[mi].name());

            ImposedMotion_[mi].write(os);

            os.endBlock();
        }

        os.endBlock();
    }
}


bool Foam::RBD::rigidBodyModel::read(const dictionary& dict)
{
     restraints_.clear();
     addRestraints(dict);


    ImposedMotion_.clear();
    addImposedMotion(dict);

    return true;
}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

Foam::Ostream& Foam::RBD::operator<<(Ostream& os, const rigidBodyModel& rbm)
{
    rbm.write(os);
    return os;
}


// ************************************************************************* //
