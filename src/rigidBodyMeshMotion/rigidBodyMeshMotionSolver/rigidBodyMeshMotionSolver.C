/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2017 OpenFOAM Foundation
    Copyright (C) 2016-2020 OpenCFD Ltd.
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

#include "rigidBodyMeshMotionSolver.H"
#include "addToRunTimeSelectionTable.H"
#include "polyMesh.H"
#include "pointPatchDist.H"
#include "pointConstraints.H"
#include "uniformDimensionedFields.H"
#include "forces.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(rigidBodyMeshMotionSolver, 0);

    addToRunTimeSelectionTable
    (
        motionSolver,
        rigidBodyMeshMotionSolver,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::rigidBodyMeshMotionSolver::bodyMesh::bodyMesh
(
    const polyMesh& mesh,
    const word& name,
    const label bodyID,
    const dictionary& dict
)
:
    name_(name),
    bodyID_(bodyID),
    patches_(dict.get<wordRes>("patches")),
    patchSet_(mesh.boundaryMesh().patchSet(patches_))
{}


Foam::rigidBodyMeshMotionSolver::rigidBodyMeshMotionSolver
(
    const polyMesh& mesh,
    const IOdictionary& dict
)
:
    motionSolver(mesh, dict, typeName),
    model_
    (
        mesh.time(),
        coeffDict(),
        IOobject
        (
            "rigidBodyMotionState",
            mesh.time().timeName(),
            "uniform",
            mesh
        ).typeHeaderOk<IOdictionary>(true)
      ? IOdictionary
        (
            IOobject
            (
                "rigidBodyMotionState",
                mesh.time().timeName(),
                "uniform",
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE,
                false
            )
        )
      : coeffDict()
    ),
    test_(coeffDict().getOrDefault("test", false)),
    rhoInf_(1.0),
    rhoName_(coeffDict().getOrDefault<word>("rho", "rho")),
    curTimeIndex_(-1),
    meshSolverPtr_
    (
        motionSolver::New
        (
            mesh,
            IOdictionary
            (
                IOobject
                (
                    "rigidBodyMotionSolver:meshSolver",
                    mesh.time().constant(),
                    mesh
                ),
                coeffDict().subDict("meshSolver")
            )
        )
    ),
    meshSolver_(refCast<displacementMotionSolver>(meshSolverPtr_()))
{
    if (rhoName_ == "rhoInf")
    {
        coeffDict().readEntry("rhoInf", rhoInf_);
    }

    const dictionary& bodiesDict = coeffDict().subDict("bodies");

    for (const entry& dEntry : bodiesDict)
    {
        const keyType& bodyName = dEntry.keyword();
        const dictionary& bodyDict = dEntry.dict();

        if (bodyDict.found("patches"))
        {
            const label bodyID = model_.bodyID(bodyName);

            if (bodyID == -1)
            {
                FatalErrorInFunction
                    << "Body " << bodyName
                    << " has been merged with another body"
                       " and cannot be assigned a set of patches"
                    << exit(FatalError);
            }

            bodyMeshes_.append
            (
                new bodyMesh
                (
                    mesh,
                    bodyName,
                    bodyID,
                    bodyDict
                )
            );
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::pointField>
Foam::rigidBodyMeshMotionSolver::curPoints() const
{
    return meshSolverPtr_->curPoints();
}


void Foam::rigidBodyMeshMotionSolver::solve()
{
    const Time& t = mesh().time();

    if (mesh().nPoints() != meshSolver_.points0().size())
    {
        FatalErrorInFunction
            << "The number of points in the mesh seems to have changed." << endl
            << "In constant/polyMesh there are " << meshSolver_.points0().size()
            << " points; in the current mesh there are " << mesh().nPoints()
            << " points." << exit(FatalError);
    }

    // Store the motion state at the beginning of the time-step
    if (curTimeIndex_ != this->db().time().timeIndex())
    {
        model_.newTime();
        curTimeIndex_ = this->db().time().timeIndex();
    }

    if (t.foundObject<uniformDimensionedVectorField>("g"))
    {
        model_.g() = t.lookupObject<uniformDimensionedVectorField>("g").value();
    }

    if (test_)
    {
        const label nIter(coeffDict().get<label>("nIter"));

        for (label i=0; i<nIter; i++)
        {
            model_.solve
            (
                t.value(),
                t.deltaTValue(),
                scalarField(model_.nDoF(), Zero),
                Field<spatialVector>(model_.nBodies(), Zero)
            );
        }
    }
    else
    {
        Field<spatialVector> fx(model_.nBodies(), Zero);


        forAll(bodyMeshes_, bi)
        {
            const label bodyID = bodyMeshes_[bi].bodyID_;

            dictionary forcesDict;
            forcesDict.add("type", functionObjects::forces::typeName);
            forcesDict.add("patches", bodyMeshes_[bi].patches_);
            forcesDict.add("rhoInf", rhoInf_);
            forcesDict.add("rho", rhoName_);
            forcesDict.add("CofR", vector::zero);

            functionObjects::forces f("forces", db(), forcesDict);
            f.calcForcesMoments();

            fx[bodyID] = spatialVector(f.momentEff(), f.forceEff());

        }

        model_.solve
        (
            t.value(),
            t.deltaTValue(),
            scalarField(model_.nDoF(), Zero),
            fx
        );

    }

    if (Pstream::master() && model_.report())
    {
        forAll(bodyMeshes_, bi)
        {
            model_.status(bodyMeshes_[bi].bodyID_);
        }
    }

    // Update the displacements
    forAll(bodyMeshes_, bi)
    {
        for (const label patchi : bodyMeshes_[bi].patchSet_)
        {
            pointField patchPoints0
            (
                meshSolver_.pointDisplacement().boundaryField()[patchi]
               .patchInternalField(meshSolver_.points0())
            );

            meshSolver_.pointDisplacement().boundaryFieldRef()[patchi] ==
            (
                model_.transformPoints
                (
                    bodyMeshes_[bi].bodyID_,
                    patchPoints0
                ) - patchPoints0
            )();
        }
    }

    meshSolverPtr_->solve();
}


bool Foam::rigidBodyMeshMotionSolver::writeObject
(
    IOstreamOption streamOpt,
    const bool valid
) const
{
    // Force ASCII writing
    streamOpt.format(IOstream::ASCII);

    IOdictionary dict
    (
        IOobject
        (
            "rigidBodyMotionState",
            mesh().time().timeName(),
            "uniform",
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        )
    );

    model_.state().write(dict);
    return dict.regIOobject::writeObject(streamOpt, valid);
}


bool Foam::rigidBodyMeshMotionSolver::read()
{
    if (motionSolver::read())
    {
        model_.read(coeffDict());

        return true;
    }

    return false;
}


void Foam::rigidBodyMeshMotionSolver::movePoints(const pointField& points)
{
    meshSolverPtr_->movePoints(points);
}


void Foam::rigidBodyMeshMotionSolver::updateMesh(const mapPolyMesh& mpm)
{
    meshSolverPtr_->updateMesh(mpm);
}


// ************************************************************************* //
