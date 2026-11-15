/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017 OpenFOAM Foundation
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

#include "singleComponentMixture.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ThermoType>
Foam::singleComponentMixture<ThermoType>::singleComponentMixture
(
    const dictionary& thermoDict,
    const fvMesh& mesh,
    const word& phaseName
)
:
    basicSpecieMixture(thermoDict, wordList::null(), mesh, phaseName),
    thermo_(thermoDict.subDict("mixture"))
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class ThermoType>
void Foam::singleComponentMixture<ThermoType>::read
(
    const dictionary& thermoDict
)
{
    thermo_ = ThermoType(thermoDict.subDict("mixture"));
}


//////////////////////////////////////////////////////////////////////////
//added 20240604
template<class ThermoType>
Foam::scalar Foam::singleComponentMixture<ThermoType>::muCellMixture
(
    const scalar p,
    const scalar T,
    const label celli
) const
{
    return 1.0;
}


template<class ThermoType>
Foam::scalar Foam::singleComponentMixture<ThermoType>::alphahCellMixture
(
    const scalar p,
    const scalar T,
    const label celli
) const
{
    return 1.0;
}


template<class ThermoType>
Foam::scalar Foam::singleComponentMixture<ThermoType>::muPatchFaceMixture
(
    const scalar p,
    const scalar T,    
    const label patchi,
    const label facei
) const
{
    return 1.0;
}


template<class ThermoType>
Foam::scalar Foam::singleComponentMixture<ThermoType>::alphahPatchFaceMixture
(
    const scalar p,
    const scalar T,
    const label patchi,
    const label facei
) const
{
    return 1.0;
}
//////////////////////////////////////////////////////////////////////////


// ************************************************************************* //
