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

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline const Foam::fvMesh& Foam::phaseSystem::mesh() const
{
    return mesh_;
}


inline const Foam::volScalarField& Foam::phaseSystem::p() const
{
    return p_;
}


inline Foam::volScalarField& Foam::phaseSystem::p()
{
    return p_;
}


inline const Foam::volScalarField& Foam::phaseSystem::alphaVoid() const
{
    return alphaVoid_;
}


inline Foam::volScalarField& Foam::phaseSystem::alphaVoid()
{
    return alphaVoid_;
}


inline const Foam::volScalarField& Foam::phaseSystem::porosity() const
{
    return porosity_;
}


inline Foam::volScalarField& Foam::phaseSystem::porosity()
{
    return porosity_;
}


inline Foam::tmp<Foam::surfaceScalarField> Foam::phaseSystem::snGradPSf() const
{
    return snGradPSf_.valid() ? tmp<surfaceScalarField>(snGradPSf_()) : tmp<surfaceScalarField>();
}


inline void Foam::phaseSystem::snGradPSf(tmp<surfaceScalarField> snGradPSf)
{
    if (!snGradPSf_.valid())
    {
        snGradPSf_ = snGradPSf;
        snGradPSf_().rename("snGradPSf");
        snGradPSf_().checkIn();
    }
    else
    {
        snGradPSf_() = snGradPSf;
    }
}


inline const Foam::phaseSystem::phaseModelList&
Foam::phaseSystem::phases() const
{
    return phaseModels_;
}


inline Foam::phaseSystem::phaseModelList&
Foam::phaseSystem::phases()
{
    return phaseModels_;
}


inline const Foam::phaseSystem::phaseModelPartialList&
Foam::phaseSystem::movingPhases() const
{
    return movingPhaseModels_;
}


inline Foam::phaseSystem::phaseModelPartialList&
Foam::phaseSystem::movingPhases()
{
    return movingPhaseModels_;
}


inline const Foam::phaseSystem::phaseModelPartialList&
Foam::phaseSystem::stationaryPhases() const
{
    return stationaryPhaseModels_;
}


inline Foam::phaseSystem::phaseModelPartialList&
Foam::phaseSystem::stationaryPhases()
{
    return stationaryPhaseModels_;
}


inline const Foam::phaseSystem::phaseModelPartialList&
Foam::phaseSystem::porousPhases() const
{
    return porousPhaseModels_;
}


inline Foam::phaseSystem::phaseModelPartialList&
Foam::phaseSystem::porousPhases()
{
    return porousPhaseModels_;
}


inline const Foam::phaseSystem::phaseModelPartialList&
Foam::phaseSystem::anisothermalPhases() const
{
    return anisothermalPhaseModels_;
}


inline Foam::phaseSystem::phaseModelPartialList&
Foam::phaseSystem::anisothermalPhases()
{
    return anisothermalPhaseModels_;
}


inline const Foam::phaseSystem::phaseModelPartialList&
Foam::phaseSystem::multiComponentPhases() const
{
    return multiComponentPhaseModels_;
}


inline Foam::phaseSystem::phaseModelPartialList&
Foam::phaseSystem::multiComponentPhases()
{
    return multiComponentPhaseModels_;
}


inline const Foam::surfaceScalarField& Foam::phaseSystem::phi() const
{
    return phi_;
}


inline Foam::surfaceScalarField& Foam::phaseSystem::phi()
{
    return phi_;
}


inline const Foam::surfaceScalarField& Foam::phaseSystem::phip() const
{
    return phip_;
}


inline Foam::surfaceScalarField& Foam::phaseSystem::phip()
{
    return phip_;
}


inline const Foam::surfaceScalarField& Foam::phaseSystem::phig() const
{
    return phig_;
}


inline Foam::surfaceScalarField& Foam::phaseSystem::phig()
{
    return phig_;
}


inline Foam::tmp<fv::convectionScheme<scalar>> Foam::phaseSystem::mvConvP() const
{
    return tmp<fv::convectionScheme<scalar>>(mvConvP_());
}


inline void Foam::phaseSystem::mvConvP(tmp<fv::convectionScheme<scalar>> mvConv)
{
    if (!mvConvP_.valid())
    {
        mvConvP_ = mvConv;
    }
    else
    {
        mvConvP_() = mvConv;
    }
}


inline Foam::tmp<fv::convectionScheme<scalar>> Foam::phaseSystem::mvConvG() const
{
    return tmp<fv::convectionScheme<scalar>>(mvConvG_());
}


inline void Foam::phaseSystem::mvConvG(tmp<fv::convectionScheme<scalar>> mvConv)
{
    if (!mvConvG_.valid())
    {
        mvConvG_ = mvConv;
    }
    else
    {
        mvConvG_() = mvConv;
    }
}


inline bool Foam::phaseSystem::hasCapillarity() const
{
    return hasCapillarity_;
}


inline void Foam::phaseSystem::setCapillarity(bool capillarity)
{
    hasCapillarity_ = capillarity;
}


// ************************************************************************* //
