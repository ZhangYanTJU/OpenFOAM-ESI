    scalar CoNum = -GREAT;

    forAll(fluidRegions, regioni)
    {
        CoNum = Foam::max
        (
            compressibleCourantNo
            (
                fluidRegions[regioni],
                runTime,
                rhoFluid[regioni],
                phiFluid[regioni]
            ),
            CoNum
        );
    }
/*
    forAll(porousFluidRegions, porousi)
    {
        CoNum = Foam::max
        (
            compressibleCourantNo
            (
                porousFluidRegions[porousi],
                runTime,
                rhoPorous[porousi],
                phiPorous[porousi]
            ),
            CoNum
        );
    }
*/
