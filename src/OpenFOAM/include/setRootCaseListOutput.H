// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Process some "standard" list options

#ifndef setRootCaseListOutput_H
#define setRootCaseListOutput_H

{
    bool listOptions = false;

    if (args.found("listSwitches"))
    {
        debug::listSwitches(args.found("listUnsetSwitches"));
        listOptions = true;
    }

    if (args.found("listRegisteredSwitches"))
    {
        debug::listRegisteredSwitches(args.found("listUnsetSwitches"));
        listOptions = true;
    }

    #ifdef fvPatchField_H
    if (args.found("listScalarBCs"))
    {
        Info<< "scalarBCs"
            << fvPatchField<Foam::scalar>::
               dictionaryConstructorTable().sortedToc()
            << endl;
        listOptions = true;
    }

    if (args.found("listVectorBCs"))
    {
        Info<< "vectorBCs"
            <<  fvPatchField<Foam::vector>::
                dictionaryConstructorTable().sortedToc()
            << endl;
        listOptions = true;
    }
    #endif

    #ifdef functionObject_H
    if (args.found("listFunctionObjects"))
    {
        Info<< "functionObjects"
            << functionObject::dictionaryConstructorTable().sortedToc()
            << endl;
        listOptions = true;
    }
    #endif

    #ifdef fvOption_H
    if (args.found("listFvOptions"))
    {
        Info<< "fvOptions"
            << fv::option::dictionaryConstructorTable().sortedToc()
            << endl;
        listOptions = true;
    }
    #endif

    #ifdef turbulentTransportModel_H
    if (args.found("listTurbulenceModels"))
    {
        Info<< "Turbulence models"
            << incompressible::turbulenceModel::
               dictionaryConstructorTable().sortedToc()
            << endl;

        Info<< "RAS models"
            << incompressible::RASModel::
               dictionaryConstructorTable().sortedToc()
            << endl;

        Info<< "LES models"
            << incompressible::LESModel::
               dictionaryConstructorTable().sortedToc()
            << endl;
        listOptions = true;
    }
    #elif defined(turbulentFluidThermoModel_H)
    if (args.found("listTurbulenceModels"))
    {
        Info<< "Turbulence models"
            << compressible::turbulenceModel::
               dictionaryConstructorTable().sortedToc()
            << endl;

        Info<< "RAS models"
            << compressible::RASModel::
               dictionaryConstructorTable().sortedToc()
            << endl;

        Info<< "LES models"
            << compressible::LESModel::
               dictionaryConstructorTable().sortedToc()
            << endl;
        listOptions = true;
    }
    #endif

    if (listOptions)
    {
        exit(0);
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
