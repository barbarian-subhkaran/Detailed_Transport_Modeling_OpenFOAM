{
    volScalarField& he = thermo.he();



    //////////////////////////////////////////////////////////////////////
    //added 20240604
    //-Reset and then calculate Heat by diffusion terms for each time step
    sum1 *= 0.0;
    sum2 *= 0.0;

    forAll(Y, k)
    {
        sum1 += fvc::laplacian((kappa/Cp)*hsi[k], Y[k]); 
        sum2 += fvc::div(J[k], hsi[k], "div(Jk,hsi)");
    }
    //////////////////////////////////////////////////////////////////////



    fvScalarMatrix EEqn
    (
        fvm::ddt(rho, he) + mvConvection->fvmDiv(phi, he)
      + fvc::ddt(rho, K) + fvc::div(phi, K)
      + (
            he.name() == "e"
          ? fvc::div
            (
                fvc::absolute(phi/fvc::interpolate(rho), U),
                p,
                "div(phiv,p)"
            )
          : -dpdt
        )
      //- fvm::laplacian(turbulence->alphaEff(), he)     //commentted 20240604
      - fvm::laplacian(kappa/Cp, he)                     //added 20240604
      + sum1                                             //added 20240604
      + sum2                                             //added 20240604
     ==
        Qdot
      + fvOptions(rho, he)
    );

    EEqn.relax();

    fvOptions.constrain(EEqn);

    EEqn.solve();

    fvOptions.correct(he);

    thermo.correct();

    Info<< "min/max(T) = "
        << min(T).value() << ", " << max(T).value() << endl;
}
