    runTime.write();

    Info<< "Qdot = " << Qdot
        << ", T = " << thermo.T()[0]
        << ", p = " << thermo.p()[0]
        << ", " << Y[0].name() << " = " << Y[0][0]
        << endl;

    post<< runTime.value() << token::TAB << thermo.T()[0] << token::TAB
        << thermo.p()[0] << endl;


//*************************************************************************//
//added 20240604
    ofs << runTime.value() << token::TAB;
    forAll(Y, i)
    {
        ofs << Y[i][0] << token::TAB;
    }
    ofs << endl; 


    Q<< runTime.value() << token::TAB << Qdot << endl;
//*************************************************************************//
