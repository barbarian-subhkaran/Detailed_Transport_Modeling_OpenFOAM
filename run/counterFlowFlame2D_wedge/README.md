cd constant
standardChemkinToFoam ../chemkin/chem.inp ../chemkin/therm.dat ../chemkin/trans.dat reactions.dat thermo.dat
cd ..
blockMesh -merge-points
checkMesh
standardReactingFoam

