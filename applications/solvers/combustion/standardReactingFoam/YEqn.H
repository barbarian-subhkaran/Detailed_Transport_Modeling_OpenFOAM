tmp<fv::convectionScheme<scalar>> mvConvection
(
    fv::convectionScheme<scalar>::New
    (
        mesh,
        fields,
        phi,
        mesh.divScheme("div(phi,Yi_h)")
    )
);



//////////////////////////////////////////////////////////////////////////
//commentted 20240604
/*
{
    reaction->correct();
    Qdot = reaction->Qdot();
    volScalarField Yt(0.0*Y[0]);

    forAll(Y, i)
    {
        if (i != inertIndex && composition.active(i))
        {
            volScalarField& Yi = Y[i];

            fvScalarMatrix YiEqn
            (
                fvm::ddt(rho, Yi)
              + mvConvection->fvmDiv(phi, Yi)
              - fvm::laplacian(turbulence->muEff(), Yi)
             ==
                reaction->R(Yi)
              + fvOptions(rho, Yi)
            );

            YiEqn.relax();

            fvOptions.constrain(YiEqn);

            YiEqn.solve(mesh.solver("Yi"));

            fvOptions.correct(Yi);

            Yi.clamp_min(0);
            Yt += Yi;
        }
    }

    Y[inertIndex] = scalar(1) - Yt;
    Y[inertIndex].clamp_min(0);
}
*/


//////////////////////////////////////////////////////////////////////////
//added 20240604
// Mixture-averaged diffusion coefficient (m2/s)
PtrList<volScalarField> Dk = dynamic_cast<const reactingMixture<gasHThermoPhysics>&>(composition).Dkmix(p,T,Y);

PtrList<volScalarField> hsi(Y.size());
PtrList<surfaceScalarField> J(Y.size());
PtrList<volVectorField> YV(Y.size());


forAll(Y, i)
{
    hsi.set
    (
        i,        
        new volScalarField
        (
            IOobject
            (
                "hsi",
                mesh.time().timeName(),
                mesh
            ),
            mesh,
            dimensionedScalar
            (
                "hsi", 
                dimEnergy/dimMass, 
                Zero
            )
        )
    );

    // Product of mass fraction and diffusion velocity of each species
    YV.set
    (
        i,
        new volVectorField
        (
            IOobject
            (
                "YV",
                mesh.time().timeName(),
                mesh
            ),
            mesh,
            dimensionedVector("YV", dimVelocity, vector(0.0, 0.0, 0.0))
        )
    );
}


forAll(p.mesh().cells(), celli)
{
    forAll(Y, i)
    {
        hsi[i][celli] = composition.Hs(i, p[celli], T[celli]);
    }
}

forAll(p.boundaryField(), patchi)
{
    forAll(p.boundaryField()[patchi], facei)
    {
        forAll(Y, i)
        {
            hsi[i].boundaryFieldRef()[patchi][facei] = composition.Hs(i, p.boundaryField()[patchi][facei], T.boundaryField()[patchi][facei]);
        }
    }
}



{
    reaction->correct();
    Qdot = reaction->Qdot();
    volScalarField Yt(0.0*Y[0]);
    


    // Compute correction velocity
    Vc *= 0.0;

    forAll(Y, i)
    {
        Vc += Dk[i]*fvc::grad(Y[i]) + (Dk[i]*Y[i]/Wmix)*fvc::grad(Wmix);
    }


    phiVc = fvc::interpolate(rho*Vc) & mesh.Sf();


    forAll(Y, i)
    {
        if (i != inertIndex && composition.active(i))
        {
            volScalarField& Yi = Y[i];

            fvScalarMatrix YiEqn
            (
                fvm::ddt(rho, Yi)
              + mvConvection->fvmDiv(phi, Yi)
              //- fvm::laplacian(turbulence->muEff(), Yi)     //commentted 20240604
              - fvm::laplacian(rho*Dk[i], Yi)                 //added 20240604
              - fvc::laplacian(rho*Dk[i]*Yi/Wmix, Wmix)       //added 20240604
              + fvm::div(phiVc, Yi, "div(phiV,Yi_h)")         //added 20240604
             ==
                reaction->R(Yi)
              + fvOptions(rho, Yi)
            );

            YiEqn.relax();

            fvOptions.constrain(YiEqn);

            YiEqn.solve(mesh.solver("Yi"));

            fvOptions.correct(Yi);

            Yi.max(0.0);
            Yt += Yi;
        }
    }

    Y[inertIndex] = scalar(1) - Yt;
    Y[inertIndex].max(0.0);
    
    
    
    volVectorField YV_balance(0.0*YV[0]);

    forAll(Y, i)
    {
        // Calculate diffusivity velocity for each species (including inert species) for use in EEqn.H
        YV[i] = (-Dk[i]*fvc::grad(Y[i])) - (Dk[i]*Y[i]/Wmix)*fvc::grad(Wmix);
    
        if (i != inertIndex && composition.active(i))
        {
            YV_balance += YV[i]; 
        }
    }
    
    // since sum(Y[i]*V[i]) over all species = 0 for global mass conservation
    YV[inertIndex] = -YV_balance;

    forAll(Y, i)
    {
        // Species diffusion flux term used in EEqn.H
        J.set
        (
            i,      
            fvc::interpolate(YV[i]*rho) & mesh.Sf() 
        );  
    } 
}
//////////////////////////////////////////////////////////////////////////
