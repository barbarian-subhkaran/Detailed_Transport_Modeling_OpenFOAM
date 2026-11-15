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

#include "specie.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //
//commentted 20240604
/*
template<class Thermo>
inline void Foam::sutherlandTransport<Thermo>::calcCoeffs
(
    const scalar mu1, const scalar T1,
    const scalar mu2, const scalar T2
)
{
    scalar rootT1 = sqrt(T1);
    scalar mu1rootT2 = mu1*sqrt(T2);
    scalar mu2rootT1 = mu2*rootT1;

    Ts_ = (mu2rootT1 - mu1rootT2)/(mu1rootT2/T1 - mu2rootT1/T2);

    As_ = mu1*(1.0 + Ts_/T1)/rootT1;
}
*/

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Thermo>
inline Foam::sutherlandTransport<Thermo>::sutherlandTransport
(
    const Thermo& t,
    //const scalar As,//commentted 20240604
    //const scalar Ts //commentted 20240604
    const scalar geometry, //added 20240604
    const scalar welldepth,//added 20240604
    const scalar diameter, //added 20240604
    const scalar dipole,   //added 20240604
    const scalar polar,    //added 20240604
    const scalar relaxation//added 20240604
)
:
    Thermo(t),
    //As_(As),//commentted 20240604
    //Ts_(Ts) //commentted 20240604
    geometry_(geometry),   //added 20240604
    welldepth_(welldepth), //added 20240604
    diameter_(diameter),   //added 20240604
    dipole_(dipole),       //added 20240604
    polar_(polar),         //added 20240604
    relaxation_(relaxation)//added 20240604
{}

//commentted 20240604
/*
template<class Thermo>
inline Foam::sutherlandTransport<Thermo>::sutherlandTransport
(
    const Thermo& t,
    const scalar mu1, const scalar T1,
    const scalar mu2, const scalar T2
)
:
    Thermo(t)
{
    calcCoeffs(mu1, T1, mu2, T2);
}
*/

template<class Thermo>
inline Foam::sutherlandTransport<Thermo>::sutherlandTransport
(
    const word& name,
    const sutherlandTransport& st
)
:
    Thermo(name, st),
    //As_(st.As_),//commentted 20240604
    //Ts_(st.Ts_) //commentted 20240604
    geometry_(st.geometry_),   //added 20240604
    welldepth_(st.welldepth_), //added 20240604
    diameter_(st.diameter_),   //added 20240604
    dipole_(st.dipole_),       //added 20240604
    polar_(st.polar_),         //added 20240604
    relaxation_(st.relaxation_)//added 20240604
{}


template<class Thermo>
inline Foam::autoPtr<Foam::sutherlandTransport<Thermo>>
Foam::sutherlandTransport<Thermo>::clone() const
{
    return autoPtr<sutherlandTransport<Thermo>>::New(*this);
}


template<class Thermo>
inline Foam::autoPtr<Foam::sutherlandTransport<Thermo>>
Foam::sutherlandTransport<Thermo>::New
(
    const dictionary& dict
)
{
    return autoPtr<sutherlandTransport<Thermo>>::New(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Thermo>
inline Foam::scalar Foam::sutherlandTransport<Thermo>::mu
(
    const scalar p,
    const scalar T
) const
{
    //return As_*::sqrt(T)/(1.0 + Ts_/T);//commentted 20240604

    scalar sigma = 1.0E-10 * diameter_;//collision diameter [m]

    //- Equation (12.8) "Chemically reacting flow" (2003)
    scalar T_star =  T / welldepth_;


    scalar epsilon = kB * welldepth_;// [J]
    scalar miu = 1.0E-21 * pow(10,-7.0/2.0) * dipole_;// [m^(3/2) J^(1/2)]


    //- Equation (12.18) "Chemically reacting flow" (2003)
    scalar delta_star = 1.0/2.0*sqr(miu)/(epsilon*pow(sigma,3.0));

    //- Equation (12.100) "Chemically reacting flow" (2003)
    return 5.0/16.0*sqrt((this->W()/1000.0)/NA*kB*T/PI)/(sqr(sigma)*Omega(2,T_star,delta_star));
}


template<class Thermo>
inline Foam::scalar Foam::sutherlandTransport<Thermo>::Omega
(
    const label index,
    const scalar T_star,
    const scalar delta_star
) const
{
    scalar omega = 0.0;

    if (index == 1)
    {
        //- Table 12.1 "Chemically reacting flow" (2003) 
        scalar a[6] = {1.0548, 0.15504, 0.55909, 2.1705, 0.093193, 1.5};

        //- Equation (12.19) "Chemically reacting flow" (2003)
        omega = a[0]*pow(T_star,-a[1])+pow(T_star+a[2],-a[3]);//Equation (12.6)
        scalar f = 1.0+(exp(a[4]/T_star)-exp(-a[5]/T_star))*sqr(delta_star)/(2.0+2.5*delta_star);

        omega *= f;
    }
    else if (index == 2)
    {
        //- Table 12.1 "Chemically reacting flow" (2003) 
        scalar b[6] = {1.0413, 0.11930, 0.43628, 1.6041, 0.095661, 2.0};

        //- Equation (12.20) "Chemically reacting flow" (2003)
        omega = b[0]*pow(T_star,-b[1])+pow(T_star+b[2],-b[3]);//Equation (12.7)
        scalar f = 1.0+(exp(b[4]/T_star)-exp(-b[5]/T_star))*sqr(delta_star)/(2.0+2.5*delta_star);

        omega *= f;
    }
    else
    {
        FatalErrorInFunction
            << "Error: collision integrals" << nl
            << exit(FatalError);
    }

    return omega;
}


template<class Thermo>
inline Foam::scalar Foam::sutherlandTransport<Thermo>::Dkk
(
    const scalar p,
    const scalar T
) const
{
    scalar sigma = 1.0E-10 * diameter_;//collision diameter [m]

    //- Equation (12.8) "Chemically reacting flow" (2003)
    scalar T_star =  T / welldepth_;


    scalar epsilon = kB * welldepth_;// [J]
    scalar miu = 1.0E-21 * pow(10,-7.0/2.0) * dipole_;// [m^(3/2) J^(1/2)]


    //- Equation (12.18) "Chemically reacting flow" (2003)
    scalar delta_star = 1.0/2.0*sqr(miu)/(epsilon*pow(sigma,3.0));

    //- Equation (12.110) "Chemically reacting flow" (2003)
    return 3.0/8.0*sqrt(NA*pow(kB*T,3.0)/(PI*(this->W()/1000.0)))/(p*sqr(sigma)*Omega(1,T_star,delta_star));
}


template<class Thermo>
inline Foam::scalar Foam::sutherlandTransport<Thermo>::F(const scalar T) const
{
    //- Equation (12.112) "Chemically reacting flow" (2003)
    return 1.0+1.0/2.0*pow(PI,3.0/2.0)*sqrt(welldepth_/T)+(sqr(PI)/4.0+2.0)*(welldepth_/T)+pow(PI,3.0/2.0)*pow(welldepth_/T,3.0/2.0);
}


template<class Thermo>
inline Foam::scalar Foam::sutherlandTransport<Thermo>::Zrot(const scalar T) const
{
    //- Equation (12.111) "Chemically reacting flow" (2003)
    return relaxation_*F(298)/F(T);
}


template<class Thermo>
inline Foam::scalar Foam::sutherlandTransport<Thermo>::kappa
(
    const scalar p, const scalar T
) const
{
    //scalar Cv_ = this->Cv(p, T);                     //commentted 20240604
    //return mu(p, T)*Cv_*(1.32 + 1.77*this->R()/Cv_); //commentted 20240604
    
    //- Transformation of the unit from [J/(kg K)] to [J/(mol K)]
    scalar Cv_total = (this->W()/1000.0)*this->Cv(p,T);

    //- Equation (8.124) "Chemically reacting flow" (2003)
    scalar Cv_trans = 3.0/2.0*R;  // translational contribution

    scalar Cv_rot = 0.0;          // rotational contribution
    scalar Cv_vib = 0.0;          // vibrational contribution

    if (geometry_ == 0)
    {
        Cv_rot = 0.0;
        Cv_vib =0.0;
    }
    else if (geometry_ == 1)
    {
        Cv_rot = R; //- Equation (8.127) "Chemically reacting flow" (2003)
        Cv_vib = Cv_total-Cv_trans-Cv_rot;//- Equation (12.108) "Chemically reacting flow" (2003)
    }
    else if (geometry_ == 2)
    {
        Cv_rot = 3.0/2.0*R;//- Equation (8.128) "Chemically reacting flow" (2003)
        Cv_vib = Cv_total-Cv_trans-Cv_rot;//- Equation (12.108) "Chemically reacting flow" (2003)
    }
    else
    {
        FatalErrorInFunction
            << "Error: geometry" << nl
            << exit(FatalError);
    }

    //- Equation (12.106) "Chemically reacting flow" (2003)
    scalar A = 5.0/2.0-this->rho(p,T)*Dkk(p,T)/mu(p,T);
    //- Equation (12.107) "Chemically reacting flow" (2003)
    scalar B = Zrot(T)+2.0/PI*(5.0/3.0*Cv_rot/R+this->rho(p,T)*Dkk(p,T)/mu(p,T));

    //- Equation (12.103) "Chemically reacting flow" (2003)
    scalar f_trans = 5.0/2.0*(1.0-2.0/PI*Cv_rot/Cv_trans*A/B);
    //- Equation (12.104) "Chemically reacting flow" (2003)
    scalar f_rot = this->rho(p,T)*Dkk(p,T)/mu(p,T)*(1.0+2.0/PI*A/B);
    //- Equation (12.105) "Chemically reacting flow" (2003)
    scalar f_vib = this->rho(p,T)*Dkk(p,T)/mu(p,T);


    //- Equation (12.102) "Chemically reacting flow" (2003)
    return mu(p,T)/(this->W()/1000.0)*(f_trans*Cv_trans+f_rot*Cv_rot+f_vib*Cv_vib);
}


template<class Thermo>
inline Foam::scalar Foam::sutherlandTransport<Thermo>::alphah
(
    const scalar p,
    const scalar T
) const
{

    return kappa(p, T)/this->Cp(p, T);
}



template<class Thermo>
inline Foam::scalar Foam::sutherlandTransport<Thermo>::geometry() const
{
    return geometry_;
}


template<class Thermo>
inline Foam::scalar Foam::sutherlandTransport<Thermo>::welldepth() const
{
    return welldepth_;
}


template<class Thermo>
inline Foam::scalar Foam::sutherlandTransport<Thermo>::diameter() const
{
    return diameter_;
}


template<class Thermo>
inline Foam::scalar Foam::sutherlandTransport<Thermo>::dipole() const
{
    return dipole_;
}


template<class Thermo>
inline Foam::scalar Foam::sutherlandTransport<Thermo>::polar() const
{
    return polar_;
}


template<class Thermo>
inline Foam::scalar Foam::sutherlandTransport<Thermo>::relaxation() const
{
    return relaxation_;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //
//added 20240604
template<class Thermo>
inline void Foam::sutherlandTransport<Thermo>::operator=
(
    const sutherlandTransport<Thermo>& st
)
{
    Thermo::operator=(st);

    geometry_ = st.geometry_;
    welldepth_ = st.welldepth_;
    diameter_ = st.diameter_;
    dipole_ = st.dipole_;
    polar_ = st.polar_;
    relaxation_ = st.relaxation_;
}


template<class Thermo>
inline void Foam::sutherlandTransport<Thermo>::operator+=
(
    const sutherlandTransport<Thermo>& st
)
{
    scalar Y1 = this->Y();

    Thermo::operator+=(st);

    if (mag(this->Y()) > SMALL)
    {
        Y1 /= this->Y();
        //scalar Y2 = st.Y()/this->Y();

        //As_ = Y1*As_ + Y2*st.As_;//commentted 20240604
        //Ts_ = Y1*Ts_ + Y2*st.Ts_;//commentted 20240604
        geometry_ = st.geometry_;    //added 20240604
        welldepth_ = st.welldepth_;  //added 20240604
        diameter_ = st.diameter_;    //added 20240604
        dipole_ = st.dipole_;        //added 20240604
        polar_ = st.polar_;          //added 20240604
        relaxation_ = st.relaxation_;//added 20240604
    }
}


template<class Thermo>
inline void Foam::sutherlandTransport<Thermo>::operator*=
(
    const scalar s
)
{
    Thermo::operator*=(s);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class Thermo>
inline Foam::sutherlandTransport<Thermo> Foam::operator+
(
    const sutherlandTransport<Thermo>& st1,
    const sutherlandTransport<Thermo>& st2
)
{
    Thermo t
    (
        static_cast<const Thermo&>(st1) + static_cast<const Thermo&>(st2)
    );

    if (mag(t.Y()) < SMALL)
    {
        return sutherlandTransport<Thermo>
        (
            t,
            //st1.As_,//commentted 20240604
            //st1.Ts_ //commentted 20240604
            st1.geometry_,  //added 20240604
            st1.welldepth_, //added 20240604
            st1.diameter_,  //added 20240604
            st1.dipole_,    //added 20240604
            st1.polar_,     //added 20240604
            st1.relaxation_ //added 20240604
        );
    }
    else
    {
        scalar Y1 = st1.Y()/t.Y();
        scalar Y2 = st2.Y()/t.Y();

        return sutherlandTransport<Thermo>
        (
            t,
            //Y1*st1.As_ + Y2*st2.As_,//commentted 20240604
            //Y1*st1.Ts_ + Y2*st2.Ts_ //commentted 20240604
            st1.geometry_,  //added 20240604
            st1.welldepth_, //added 20240604
            st1.diameter_,  //added 20240604
            st1.dipole_,    //added 20240604
            st1.polar_,     //added 20240604
            st1.relaxation_ //added 20240604
        );
    }
}


template<class Thermo>
inline Foam::sutherlandTransport<Thermo> Foam::operator*
(
    const scalar s,
    const sutherlandTransport<Thermo>& st
)
{
    return sutherlandTransport<Thermo>
    (
        s*static_cast<const Thermo&>(st),
        //st.As_,//commentted 20240604
        //st.Ts_ //commentted 20240604
        st.geometry_,  //added 20240604
        st.welldepth_, //added 20240604
        st.diameter_,  //added 20240604
        st.dipole_,    //added 20240604
        st.polar_,     //added 20240604
        st.relaxation_ //added 20240604
    );
}


// ************************************************************************* //
