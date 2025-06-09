/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2018 OpenFOAM Foundation
    Copyright (C) 2025 OpenCFD Ltd.
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

#include "smoothTriSurfaceMesh.H"
#include "addToRunTimeSelectionTable.H"
#include "unitConversion.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(smoothTriSurfaceMesh, 0);
addToRunTimeSelectionTable(searchableSurface, smoothTriSurfaceMesh, dict);

}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::smoothTriSurfaceMesh::calcFeatureEdges(const scalar featureAngle)
{
    if (featureAngle <= -180)
    {
        return;
    }

    const triSurface& s = *this;

    scalar cosAngle = Foam::cos(degToRad(featureAngle));

    const labelListList& edgeFaces = s.edgeFaces();
    const edgeList& edges = s.edges();

    forAll(edgeFaces, edgei)
    {
        const labelList& eFaces = edgeFaces[edgei];

        if (eFaces.size() > 2)
        {
            isBorderEdge_[edgei] = true;
        }
        else if (eFaces.size() == 2)
        {
            const vector n0(s[eFaces[0]].unitNormal(points()));
            const vector n1(s[eFaces[1]].unitNormal(points()));
            if ((n0&n1) < cosAngle)
            {
                isBorderEdge_[edgei] = true;
            }
        }
    }

    forAll(isBorderEdge_, edgei)
    {
        if (isBorderEdge_[edgei])
        {
            const edge& e = edges[edgei];
            isPointOnBorderEdge_[e[0]] = true;
            isPointOnBorderEdge_[e[1]] = true;
        }
    }
}


Foam::vector Foam::smoothTriSurfaceMesh::pointNormal
(
    const label startFacei,
    const label localPointi
) const
{
    const triSurface& s = *this;

    if (!isPointOnBorderEdge_[localPointi])
    {
        const labelList& pFaces = pointFaces()[localPointi];
        vector pn(Zero);
        forAll(pFaces, i)
        {
            label facei = pFaces[i];
            pn += s[facei].unitNormal(points());
        }
        return pn.normalise();
    }


    // Calculate the local point normal on the face. This routine
    // only gets called if the point is on a border edge so we can
    // walk and always hit a border edge.

    const edgeList& edges = triSurface::edges();
    const labelList& fEdges = faceEdges()[startFacei];

    // Get the two edges on the point
    label e0 = -1;
    label e1 = -1;
    forAll(fEdges, i)
    {
        const edge& e = edges[fEdges[i]];
        if (e.otherVertex(localPointi) == -1)
        {
            e0 = fEdges[fEdges.rcIndex(i)];
            e1 = fEdges[fEdges.fcIndex(i)];
            break;
        }
    }

    label facei = startFacei;
    label edgei = e0;

    // Initialise normal
    vector n(s[facei].unitNormal(points()));

    while (!isBorderEdge_[edgei])
    {
        // Cross edge to next face
        const labelList& eFaces = edgeFaces()[edgei];
        if (eFaces.size() != 2)
        {
            break;
        }

        forAll(eFaces, i)
        {
            if (eFaces[i] != facei)
            {
                facei = eFaces[i];
                break;
            }
        }

        if (facei == startFacei)
        {
            break;
        }

        n += s[facei].unitNormal(points());

        // Cross face to next edge
        const labelList& fEdges = faceEdges()[facei];

        forAll(fEdges, i)
        {
            label ei = fEdges[i];
            if (ei != edgei && edges[ei].otherVertex(localPointi) != -1)
            {
                edgei = ei;
                break;
            }
        }
    }


    // Walk in other direction
    {

        facei = startFacei;
        edgei = e1;

        while (!isBorderEdge_[edgei])
        {
            // Cross edge to next face
            const labelList& eFaces = edgeFaces()[edgei];
            if (eFaces.size() != 2)
            {
                break;
            }

            forAll(eFaces, i)
            {
                if (eFaces[i] != facei)
                {
                    facei = eFaces[i];
                    break;
                }
            }

            if (facei == startFacei)
            {
                break;
            }

            n += s[facei].unitNormal(points());

            // Cross face to next edge
            const labelList& fEdges = faceEdges()[facei];

            forAll(fEdges, i)
            {
                label ei = fEdges[i];
                if (ei != edgei && edges[ei].otherVertex(localPointi) != -1)
                {
                    edgei = ei;
                    break;
                }
            }
        }
    }

    return n.normalise();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::smoothTriSurfaceMesh::smoothTriSurfaceMesh
(
    const IOobject& io,
    const scalar featureAngle
)
:
    triSurfaceMesh(io),
    isBorderEdge_(nEdges()),
    isPointOnBorderEdge_(nPoints())
{
    calcFeatureEdges(featureAngle);
}


Foam::smoothTriSurfaceMesh::smoothTriSurfaceMesh
(
    const IOobject& io,
    const dictionary& dict
)
:
    triSurfaceMesh(io, dict),
    isBorderEdge_(nEdges()),
    isPointOnBorderEdge_(nPoints())
{
    if (dict.found("featureAngle"))
    {
        calcFeatureEdges(readScalar(dict.lookup("featureAngle")));
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::smoothTriSurfaceMesh::~smoothTriSurfaceMesh()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::smoothTriSurfaceMesh::getNormal
(
    const List<pointIndexHit>& info,
    vectorField& normal
) const
{
    const triSurface& s = static_cast<const triSurface&>(*this);
    const pointField& pts = s.points();

    DebugPout<< "smoothTriSurfaceMesh::getNormal :"
        << " surface:" << this->searchableSurface::name()
        << " tris:" << s.size()
        << " feat edges:" << isBorderEdge_.count()
        << endl;

    normal.setSize(info.size());

    forAll(info, i)
    {
        if (info[i].hit())
        {
            label facei = info[i].index();

            // Get local coordinates in triangle
            barycentric2D coordinates
            (
                s[facei].tri(pts).pointToBarycentric(info[i].hitPoint())
            );

            // Average point normals
            const triFace& localTri = s.localFaces()[facei];

            vector n(Zero);
            forAll(localTri, fp)
            {
                n += coordinates[fp]*pointNormal(facei, localTri[fp]);
            }
            normal[i] = n.normalise();
        }
        else
        {
            // Set to what?
            normal[i] = Zero;
        }
    }
}


// ************************************************************************* //
