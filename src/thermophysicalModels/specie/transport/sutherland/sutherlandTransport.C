/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
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

#include "sutherlandTransport.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Thermo>
Foam::scalar Foam::sutherlandTransport<Thermo>::readCoeff
(
    const word& coeffName,
    const dictionary& dict
)
{
    return dict.subDict("transport").get<scalar>(coeffName);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Thermo>
Foam::sutherlandTransport<Thermo>::sutherlandTransport(const dictionary& dict)
:
    Thermo(dict),
    //As_(readCoeff("As", dict)),//commentted 20240604
    //Ts_(readCoeff("Ts", dict)) //commentted 20240604
    geometry_(readCoeff("geometry", dict)),   //added 20240604
    welldepth_(readCoeff("welldepth", dict)), //added 20240604
    diameter_(readCoeff("diameter", dict)),   //added 20240604
    dipole_(readCoeff("dipole", dict)),       //added 20240604
    polar_(readCoeff("polar", dict)),         //added 20240604
    relaxation_(readCoeff("relaxation", dict))//added 20240604
{}


template<class Thermo>
Foam::sutherlandTransport<Thermo>::sutherlandTransport
(
    const Thermo& t,
    const dictionary& dict
)
:
    Thermo(t),
    //As_(readCoeff("As", dict)),//commentted 20240604
    //Ts_(readCoeff("Ts", dict)) //commentted 20240604
    geometry_(readCoeff("geometry", dict)),   //added 20240604
    welldepth_(readCoeff("welldepth", dict)), //added 20240604
    diameter_(readCoeff("diameter", dict)),   //added 20240604
    dipole_(readCoeff("dipole", dict)),       //added 20240604
    polar_(readCoeff("polar", dict)),         //added 20240604
    relaxation_(readCoeff("relaxation", dict))//added 20240604
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Thermo>
void Foam::sutherlandTransport<Thermo>::write(Ostream& os) const
{
    os.beginBlock(this->specie::name());

    Thermo::write(os);

    // Entries in dictionary format
    {
        os.beginBlock("transport");
        //os.writeEntry("As", As_);//commentted 20240604
        //os.writeEntry("Ts", Ts_);//commentted 20240604
        os.writeEntry("geometry", geometry_);    //added 20240604
        os.writeEntry("welldepth", welldepth_);  //added 20240604
        os.writeEntry("diameter", diameter_);    //added 20240604
        os.writeEntry("dipole", dipole_);        //added 20240604
        os.writeEntry("polar", polar_);          //added 20240604
        os.writeEntry("relaxation", relaxation_);//added 20240604
        os.endBlock();
    }

    os.endBlock();
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class Thermo>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const sutherlandTransport<Thermo>& st
)
{
    st.write(os);
    return os;
}


// ************************************************************************* //
