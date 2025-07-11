/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     5.0
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

Global
    CourantNo

Description
    Calculates and outputs the maximum Courant Number.

\*---------------------------------------------------------------------------*/

scalar CoNum = 0.0;

surfaceScalarField magPhi = mag(phi)/fvc::interpolate(fluid.alphaVoid());

forAll(phases, phasei)
{
    if (phases[phasei].stationary()) continue;

    magPhi = max
    (
        magPhi,
        mag(phases[phasei].phi())
    );
}

surfaceScalarField SfUfbyDelta =
    mesh.surfaceInterpolation::deltaCoeffs()*magPhi;

const scalar deltaT = runTime.deltaT().value();

CoNum = max(SfUfbyDelta/mesh.magSf()).value()*deltaT;

Info<< "Courant Number max: " << CoNum << endl;

// ************************************************************************* //
