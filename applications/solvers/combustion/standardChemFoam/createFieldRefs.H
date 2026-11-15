BasicChemistryModel<psiReactionThermo>& chemistry = pChemistry();
scalar dtChem = refCast<const BasicChemistryModel<psiReactionThermo>>
(
    chemistry
).deltaTChem()[0];
basicSpecieMixture& composition = thermo.composition();
PtrList<volScalarField>& Y = composition.Y();
volScalarField& p = thermo.p();


//*************************************************************************//
//added 20240604
OFstream ofs(args.path()/"standardChemFoamYi.out");

ofs << "# Time" << token::TAB;
forAll(Y, i)
{
    ofs << Y[i].name() << token::TAB;
}
ofs << endl; 


OFstream Q(args.path()/"standardChemFoamQdot.out");
Q<< "# Time" << token::TAB << "Qdot" << endl;
//*************************************************************************//