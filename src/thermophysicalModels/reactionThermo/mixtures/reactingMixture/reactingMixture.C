/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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

#include "reactingMixture.H"
#include "fvMesh.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class ThermoType>
Foam::autoPtr<Foam::chemistryReader<ThermoType>>&
Foam::reactingMixture<ThermoType>::reader()
{
    return static_cast<autoPtr<chemistryReader<ThermoType>>&>(*this);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ThermoType>
Foam::reactingMixture<ThermoType>::reactingMixture
(
    const dictionary& thermoDict,
    const fvMesh& mesh,
    const word& phaseName
)
:
    speciesTable(),
    autoPtr<chemistryReader<ThermoType>>
    (
        chemistryReader<ThermoType>::New(thermoDict, *this)
    ),
    multiComponentMixture<ThermoType>
    (
        thermoDict,
        *this,
        this->reader()->speciesThermo(),
        mesh,
        phaseName
    ),
    PtrList<Reaction<ThermoType>>
    (
        this->reader()->reactions()
    ),
    speciesComposition_
    (
        this->reader()->specieComposition()
    )
{
    // Done with reader
    reader().clear();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ThermoType>
void Foam::reactingMixture<ThermoType>::read(const dictionary& thermoDict)
{}


//////////////////////////////////////////////////////////////////////////
//- Equation (12.180) "Chemically reacting flow" (2003)
template<class ThermoType>
inline Foam::PtrList<Foam::volScalarField> Foam::reactingMixture<ThermoType>::Dkmix
(
    const volScalarField& p, 
    const volScalarField& T, 
    const PtrList<Foam::volScalarField>& Y
) const
{       
    label nsp_ = speciesComposition_.size();

    PtrList<volScalarField> Dimix(nsp_);
    
    forAll(Y, i)
    {
        Dimix.set
        (
            i,        
            new volScalarField
            (
                IOobject
                (
                    "Dimix",
                    p.mesh().time().timeName(),
                    p.mesh()
                ),
                p.mesh(),
                dimensionedScalar
                (
                    "Dimix", 
                    dimensionSet(0,2,-1,0,0,0,0), 
                    Zero
                )
            )        
        );              
    }
    
    List<scalar> XNew(nsp_, 0.0);
    List<scalar> YNew(nsp_, 0.0);

    // for internal field
    forAll(p.mesh().cells(), celli)
    {
        List<scalar> mixDiffCoeffs(nsp_, 0.0);

        scalar sumY = 0.0;
    
        forAll (Y, i)
        {
            YNew[i] = Y[i][celli];
            sumY += YNew[i];
        }

        if (sumY != 0.0)
        {
            scalar inverseMoleWeight = 0.0; 

            forAll (Y, i)
            {
                YNew[i] /= sumY;

                inverseMoleWeight += YNew[i]/this->W(i);
            }

            scalar meanMoleWeight = 1.0/inverseMoleWeight;// [kg/kmol]

            forAll (Y, i)
            {
                //- Equation (3.84) "Chemically reacting flow" (2003)
                XNew[i] = YNew[i]*meanMoleWeight/this->W(i);
            }


            List<scalar> sum(nsp_, 0.0);
            label flag = -1;

            forAll (Y, k)
            {
                forAll (Y, j)
                {
                    if (j != k)
                    {
                        sum[k] += XNew[j]/this->Dkj(k,j,p[celli],T[celli]);
                    }
                }

                if (YNew[k] == 1.0)
                {
                    flag = k;
                }
            }

            if (flag != -1)
            {
                mixDiffCoeffs[flag] = this->speciesData()[flag].Dkk(p[celli],T[celli]);
            }
            else
            {
                forAll (Y, k)
                {
                    //- Equation (12.180) "Chemically reacting flow" (2003)
                    mixDiffCoeffs[k] = (1.0-YNew[k])/sum[k];
                }
            }
        }


        forAll (Y, i)
        {
            Dimix[i][celli] = mixDiffCoeffs[i];
        }
    }


    // for boundary fields         
    forAll(p.boundaryField(), patchi)
    {
        forAll(p.boundaryField()[patchi], facei)
        {
            List<scalar> mixDiffCoeffs(nsp_, 0.0);

            scalar sumY = 0.0;
        
            forAll (Y, i)
            {
                YNew[i] = Y[i].boundaryField()[patchi][facei];
                sumY += YNew[i];
            }

            if (sumY != 0.0)
            {
                scalar inverseMoleWeight = 0.0; 

                forAll (Y, i)
                {
                    YNew[i] /= sumY;

                    inverseMoleWeight += YNew[i]/this->W(i);
                }

                scalar meanMoleWeight = 1.0/inverseMoleWeight;// [kg/kmol]

                forAll (Y, i)
                {
                    //- Equation (3.84) "Chemically reacting flow" (2003)
                    XNew[i] = YNew[i]*meanMoleWeight/this->W(i);
                }


                List<scalar> sum(nsp_, 0.0);
                label flag = -1;

                forAll (Y, k)
                {
                    forAll (Y, j)
                    {
                        if (j != k)
                        {
                            sum[k] += XNew[j]/this->Dkj(k,j,p.boundaryField()[patchi][facei],T.boundaryField()[patchi][facei]);
                        }
                    }

                    if (YNew[k] == 1.0)
                    {
                        flag = k;
                    }
                }


                if (flag != -1)
                {
                    mixDiffCoeffs[flag] = this->speciesData()[flag].Dkk(p.boundaryField()[patchi][facei],T.boundaryField()[patchi][facei]);
                }
                else
                {
                    forAll (Y, k)
                    {
                        //- Equation (12.180) "Chemically reacting flow" (2003)
                        mixDiffCoeffs[k] = (1.0-YNew[k])/sum[k];
                    }
                }
            }


            forAll (Y, i)
            {
                Dimix[i].boundaryFieldRef()[patchi][facei] = mixDiffCoeffs[i];
            }
        }
    }

    return Dimix;
}
//////////////////////////////////////////////////////////////////////////


// ************************************************************************* //
