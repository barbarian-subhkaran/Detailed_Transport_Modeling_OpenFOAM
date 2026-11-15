/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2021-2022 OpenCFD Ltd.
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

#include "multiComponentMixture.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class ThermoType>
const ThermoType& Foam::multiComponentMixture<ThermoType>::constructSpeciesData
(
    const dictionary& thermoDict
)
{
    forAll(species_, i)
    {
        speciesData_.set
        (
            i,
            new ThermoType(thermoDict.subDict(species_[i]))
        );
    }

    return speciesData_[0];
}


template<class ThermoType>
void Foam::multiComponentMixture<ThermoType>::correctMassFractions()
{
    // Multiplication by 1.0 changes Yt patches to "calculated"
    volScalarField Yt("Yt", 1.0*Y_[0]);

    for (label n=1; n<Y_.size(); n++)
    {
        Yt += Y_[n];
    }

    if (mag(min(Yt).value()) < ROOTVSMALL)
    {
        FatalErrorInFunction
            << "Sum of mass fractions is zero for species " << this->species()
            << nl
            << incrIndent << indent
            << "Min of mass fraction sum " << min(Yt).value()
            << decrIndent
            << exit(FatalError);
    }

    const scalar diff(mag(max(Yt).value()) - scalar(1));

    if (diff > ROOTVSMALL)
    {
        WarningInFunction
            << "Sum of mass fractions is different from one for species "
            << this->species()
            << nl << incrIndent << indent
            << "Max of mass fraction sum differs from 1 by " << diff
            << decrIndent << nl;
    }

    forAll(Y_, n)
    {
        Y_[n] /= Yt;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ThermoType>
Foam::multiComponentMixture<ThermoType>::multiComponentMixture
(
    const dictionary& thermoDict,
    const wordList& specieNames,
    const ReactionTable<ThermoType>& thermoData,
    const fvMesh& mesh,
    const word& phaseName
)
:
    basicSpecieMixture(thermoDict, specieNames, mesh, phaseName),
    speciesData_(species_.size()),
    mixture_("mixture", *thermoData[specieNames[0]]),
    mixtureVol_("volMixture", *thermoData[specieNames[0]])
{
    forAll(species_, i)
    {
        speciesData_.set
        (
            i,
            new ThermoType(*thermoData[species_[i]])
        );
    }

    //added 20240604
    Info << "\nSelecting mixture-averaged transport model" << nl << endl;

    correctMassFractions();
}


template<class ThermoType>
Foam::multiComponentMixture<ThermoType>::multiComponentMixture
(
    const dictionary& thermoDict,
    const fvMesh& mesh,
    const word& phaseName
)
:
    basicSpecieMixture
    (
        thermoDict,
        thermoDict.lookup("species"),
        mesh,
        phaseName
    ),
    speciesData_(species_.size()),
    mixture_("mixture", constructSpeciesData(thermoDict)),
    mixtureVol_("volMixture", speciesData_[0])
{
    //added 20240604
    Info << "\nSelecting mixture-averaged transport model" << nl << endl;

    correctMassFractions();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ThermoType>
const ThermoType& Foam::multiComponentMixture<ThermoType>::cellMixture
(
    const label celli
) const
{
    mixture_ = Y_[0][celli]*speciesData_[0];

    for (label n=1; n<Y_.size(); n++)
    {
        mixture_ += Y_[n][celli]*speciesData_[n];
    }

    return mixture_;
}


template<class ThermoType>
const ThermoType& Foam::multiComponentMixture<ThermoType>::patchFaceMixture
(
    const label patchi,
    const label facei
) const
{
    mixture_ = Y_[0].boundaryField()[patchi][facei]*speciesData_[0];

    for (label n=1; n<Y_.size(); n++)
    {
        mixture_ += Y_[n].boundaryField()[patchi][facei]*speciesData_[n];
    }

    return mixture_;
}


template<class ThermoType>
const ThermoType& Foam::multiComponentMixture<ThermoType>::cellVolMixture
(
    const scalar p,
    const scalar T,
    const label celli
) const
{
    scalar rhoInv = 0.0;
    forAll(speciesData_, i)
    {
        rhoInv += Y_[i][celli]/speciesData_[i].rho(p, T);
    }

    mixtureVol_ =
        Y_[0][celli]/speciesData_[0].rho(p, T)/rhoInv*speciesData_[0];

    for (label n=1; n<Y_.size(); n++)
    {
        mixtureVol_ +=
            Y_[n][celli]/speciesData_[n].rho(p, T)/rhoInv*speciesData_[n];
    }

    return mixtureVol_;
}


template<class ThermoType>
const ThermoType& Foam::multiComponentMixture<ThermoType>::
patchFaceVolMixture
(
    const scalar p,
    const scalar T,
    const label patchi,
    const label facei
) const
{
    scalar rhoInv = 0.0;
    forAll(speciesData_, i)
    {
        rhoInv +=
            Y_[i].boundaryField()[patchi][facei]/speciesData_[i].rho(p, T);
    }

    mixtureVol_ =
        Y_[0].boundaryField()[patchi][facei]/speciesData_[0].rho(p, T)/rhoInv
      * speciesData_[0];

    for (label n=1; n<Y_.size(); n++)
    {
        mixtureVol_ +=
            Y_[n].boundaryField()[patchi][facei]/speciesData_[n].rho(p,T)
          / rhoInv*speciesData_[n];
    }

    return mixtureVol_;
}


template<class ThermoType>
void Foam::multiComponentMixture<ThermoType>::read
(
    const dictionary& thermoDict
)
{
    forAll(species_, i)
    {
        speciesData_[i] = ThermoType(thermoDict.subDict(species_[i]));
    }
}


//////////////////////////////////////////////////////////////////////////
//added 20240604
template<class ThermoType>
Foam::scalar Foam::multiComponentMixture<ThermoType>::Dkj
(
    const label k,
    const label j,
    const scalar p,
    const scalar T
) const
{
    scalar Wk = speciesData_[k].W();// [kg/kmol]
    scalar Wj = speciesData_[j].W();// [kg/kmol]
    //- Equation (12.92) "Chemically reacting flow" (2003)
    scalar Wkj = Wk*Wj/(Wk+Wj);

    scalar epsilonk = kB * speciesData_[k].welldepth();// [J]
    scalar epsilonj = kB * speciesData_[j].welldepth();// [J]
    //- Equation (12.36) "Chemically reacting flow" (2003)
    scalar epsilonkj = sqrt(epsilonk*epsilonj);

    scalar sigmak = 1.0E-10 * speciesData_[k].diameter();// [m]
    scalar sigmaj = 1.0E-10 * speciesData_[j].diameter();// [m]
    //- Equation (12.35) "Chemically reacting flow" (2003)
    scalar sigmakj = 1.0/2.0*(sigmak+sigmaj);

    scalar miuk = 1.0E-21 * pow(10,-7.0/2.0) * speciesData_[k].dipole();// [m^(3/2) J^(1/2)]
    scalar miuj = 1.0E-21 * pow(10,-7.0/2.0) * speciesData_[j].dipole();// [m^(3/2) J^(1/2)]

    scalar polark = 1.0E-30 * speciesData_[k].polar();// [m]
    scalar polarj = 1.0E-30 * speciesData_[j].polar();// [m]


    scalar delta_star = 0.0;

    if (polark == 0.0)//k-polar
    {
        if (polarj == 0.0)//j-polar
        {
            //- Equation (12.37) "Chemically reacting flow" (2003)
            delta_star = 1.0/2.0*miuk*miuj/(epsilonkj*pow(sigmakj,3.0));
        }
        else//j-nonpolar
        {
            //- Equation (12.41) "Chemically reacting flow" (2003)
            scalar alphajstar = polarj/pow(sigmaj,3.0);//non-polar
            //- Equation (12.42) "Chemically reacting flow" (2003)
            scalar miukstar = miuk/sqrt(epsilonk*pow(sigmak,3.0));//polar
            //- Equation (12.40) "Chemically reacting flow" (2003)
            scalar ksi = 1.0+1.0/4.0*alphajstar*miukstar*sqrt(epsilonk/epsilonj);
            //- Equation (12.38) "Chemically reacting flow" (2003)
            epsilonkj *= sqr(ksi);
            //- Equation (12.39) "Chemically reacting flow" (2003)
            sigmakj *= pow(ksi, -1.0/6.0);
            //- Equation (12.37) "Chemically reacting flow" (2003)
            delta_star = 1.0/2.0*miuk*miuj/(epsilonkj*pow(sigmakj,3.0));
        }
    }
    else//k-nonpolar
    {
        if (polarj == 0.0)//j-polar
        {
            //- Equation (12.41) "Chemically reacting flow" (2003)
            scalar alphakstar = polark/pow(sigmak,3.0);//non-polar
            //- Equation (12.42) "Chemically reacting flow" (2003)
            scalar miujstar = miuj/sqrt(epsilonj*pow(sigmaj,3.0));//polar
            //- Equation (12.40) "Chemically reacting flow" (2003)
            scalar ksi = 1.0+1.0/4.0*alphakstar*miujstar*sqrt(epsilonj/epsilonk);
            //- Equation (12.38) "Chemically reacting flow" (2003)
            epsilonkj *= sqr(ksi);
            //- Equation (12.39) "Chemically reacting flow" (2003)
            sigmakj *= pow(ksi, -1.0/6.0);
            //- Equation (12.37) "Chemically reacting flow" (2003)
            delta_star = 1.0/2.0*miuk*miuj/(epsilonkj*pow(sigmakj,3.0));
        }
        else//j-nonpolar
        {
            delta_star = 0.0;
        }
    }

    //- Equation (12.8) "Chemically reacting flow" (2003)
    scalar T_star =  T * kB / epsilonkj;
    

    //- Table 12.1 "Chemically reacting flow" (2003) 
    scalar a[6] = {1.0548, 0.15504, 0.55909, 2.1705, 0.093193, 1.5};

    //- Equation (12.19) "Chemically reacting flow" (2003)
    scalar omega = a[0]*pow(T_star,-a[1])+pow(T_star+a[2],-a[3]);//Equation (12.6)
    scalar f = 1.0+(exp(a[4]/T_star)-exp(-a[5]/T_star))*sqr(delta_star)/(2.0+2.5*delta_star);

    omega *= f;


    //- Equation (12.113) "Chemically reacting flow" (2003)
    return 3.0/16.0*sqrt(2.0*NA*pow(kB*T,3.0)/(PI*Wkj/1000.0))/(p*sqr(sigmakj)*omega);
}


template<class ThermoType>
Foam::scalar Foam::multiComponentMixture<ThermoType>::Phikj
(
    const label k,
    const label j,
    const scalar p,
    const scalar T
) const
{
    scalar Wk = speciesData_[k].W();// [kg/kmol]
    scalar Wj = speciesData_[j].W();// [kg/kmol]

    //- Equation (12.100) "Chemically reacting flow" (2003)
    scalar muk = speciesData_[k].mu(p,T);
    scalar muj = speciesData_[j].mu(p,T);

    //- Equation (12.118) "Chemically reacting flow" (2003)
    return 1.0/sqrt(8.0)/sqrt(1.0+Wk/Wj) * sqr(1.0+sqrt(muk/muj)*pow(Wj/Wk,1.0/4.0));
}


template<class ThermoType>
Foam::scalar Foam::multiComponentMixture<ThermoType>::muCellMixture
(
    const scalar p,
    const scalar T,
    const label celli
) const
{
    scalar muMixture = 0.0;

    List<scalar> X(species_.size(),0.0);//mole fraction
    List<scalar> Y(species_.size(),0.0);//mass fraction
    scalar inverseMoleWeight = 0.0; 

    forAll(species_, i)
    {
        Y[i] = Y_[i][celli];
        inverseMoleWeight += Y[i]/speciesData_[i].W();
    }

    if (inverseMoleWeight == 0.0)
    {
        return 0.0;
    }

    scalar meanMoleWeight = 1.0/inverseMoleWeight;// [kg/kmol]

    forAll(species_, i)
    {
        //- Equation (3.84) "Chemically reacting flow" (2003)
        X[i] = Y[i]*meanMoleWeight/speciesData_[i].W();
    }

    List<scalar> sum(species_.size(),0.0);

    forAll(species_, k)
    {
        forAll(species_, j)
        {
            sum[k] += X[j]*Phikj(k,j,p,T);
        }
    }


    forAll(species_, k)
    {
        //- Equation (12.117) "Chemically reacting flow" (2003)
        muMixture += X[k]*speciesData_[k].mu(p,T)/sum[k];
    }
    
    return muMixture;
}


template<class ThermoType>
Foam::scalar Foam::multiComponentMixture<ThermoType>::alphahCellMixture
(
    const scalar p,
    const scalar T,
    const label celli
) const
{
    scalar alphahMixture = 0.0;

    List<scalar> X(species_.size(),0.0);//mole fraction
    List<scalar> Y(species_.size(),0.0);//mass fraction
    scalar inverseMoleWeight = 0.0; 

    forAll(species_, i)
    {
        Y[i] = Y_[i][celli];
        inverseMoleWeight += Y[i]/speciesData_[i].W();
    }

    if (inverseMoleWeight == 0.0)
    {
        return 0.0;
    }

    scalar meanMoleWeight = 1.0/inverseMoleWeight;

    forAll(species_, i)
    {
        //- Equation (3.84) "Chemically reacting flow" (2003)
        X[i] = Y[i]*meanMoleWeight/speciesData_[i].W();
    }


    scalar sum1 = 0.0;
    scalar sum2 = 0.0;
    scalar Cpmix = 0.0;// [J/kg/K]

    forAll(species_, k)
    {
        sum1 += X[k]*speciesData_[k].kappa(p,T);
        sum2 += X[k]/speciesData_[k].kappa(p,T);
        //- Equation (1.12) "Theoretical and numerical combustion" (2012)        
        Cpmix += Y[k]*speciesData_[k].Cp(p,T);
    }

    //- Equation (12.119) "Chemically reacting flow" (2003)
    scalar kappa = 1.0/2.0*(sum1+1.0/sum2);

    alphahMixture = kappa/Cpmix;
    
    return alphahMixture;
}


template<class ThermoType>
Foam::scalar Foam::multiComponentMixture<ThermoType>::muPatchFaceMixture
(
    const scalar p,
    const scalar T,    
    const label patchi,
    const label facei
) const
{
    scalar muMixture = 0.0;

    List<scalar> X(species_.size(),0.0);//mole fraction
    List<scalar> Y(species_.size(),0.0);//mass fraction
    scalar inverseMoleWeight = 0.0; 

    forAll(species_, i)
    {
        Y[i] = Y_[i].boundaryField()[patchi][facei];
        inverseMoleWeight += Y[i]/speciesData_[i].W();
    }

    if (inverseMoleWeight == 0.0)
    {
        return 0.0;
    }

    scalar meanMoleWeight = 1.0/inverseMoleWeight;// [kg/kmol]

    forAll(species_, i)
    {
        //- Equation (3.84) "Chemically reacting flow" (2003)
        X[i] = Y[i]*meanMoleWeight/speciesData_[i].W();
    }

    List<scalar> sum(species_.size(),0.0);

    forAll(species_, k)
    {
        forAll(species_, j)
        {
            sum[k] += X[j]*Phikj(k,j,p,T);
        }
    }


    forAll(species_, k)
    {
        //- Equation (12.117) "Chemically reacting flow" (2003)
        muMixture += X[k]*speciesData_[k].mu(p,T)/sum[k];
    }
    
    return muMixture;
}


template<class ThermoType>
Foam::scalar Foam::multiComponentMixture<ThermoType>::alphahPatchFaceMixture
(
    const scalar p,
    const scalar T,
    const label patchi,
    const label facei
) const
{
    scalar alphahMixture = 0.0;

    List<scalar> X(species_.size(),0.0);//mole fraction
    List<scalar> Y(species_.size(),0.0);//mass fraction
    scalar inverseMoleWeight = 0.0; 

    forAll(species_, i)
    {
        Y[i] = Y_[i].boundaryField()[patchi][facei];
        inverseMoleWeight += Y[i]/speciesData_[i].W();
    }

    if (inverseMoleWeight == 0.0)
    {
        return 0.0;
    }

    scalar meanMoleWeight = 1.0/inverseMoleWeight;

    forAll(species_, i)
    {
        //- Equation (3.84) "Chemically reacting flow" (2003)
        X[i] = Y[i]*meanMoleWeight/speciesData_[i].W();
    }


    scalar sum1 = 0.0;
    scalar sum2 = 0.0;
    scalar Cpmix = 0.0;// [J/kg/K]

    forAll(species_, k)
    {
        sum1 += X[k]*speciesData_[k].kappa(p,T);
        sum2 += X[k]/speciesData_[k].kappa(p,T);
        //- Equation (1.12) "Theoretical and numerical combustion" (2012) 
        Cpmix += Y[k]*speciesData_[k].Cp(p,T);
    }

    //- Equation (12.119) "Chemically reacting flow" (2003)
    scalar kappa = 1.0/2.0*(sum1+1.0/sum2);

    alphahMixture = kappa/Cpmix;

    return alphahMixture;
}
//////////////////////////////////////////////////////////////////////////


// ************************************************************************* //
