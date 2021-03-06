/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
     \\/     M anipulation  |
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

#include "BinghamPapanastasiou.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace viscosityModels
{
    defineTypeNameAndDebug(BinghamPapanastasiou, 0);

    addToRunTimeSelectionTable
    (
        viscosityModel,
        BinghamPapanastasiou,
        dictionary
    );
}
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::viscosityModels::BinghamPapanastasiou::calcNu() const
{
    tmp<volScalarField> sr(strainRate());

    return
    (
	 	k_ + tau0_ * (1.0 - exp(-sr() / (sqrt(2.0) * eps_)))
		/ (max(sr(), dimensionedScalar ("VSMALL", dimless/dimTime, VSMALL)))
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::viscosityModels::BinghamPapanastasiou::BinghamPapanastasiou
(
    const word& name,
    const dictionary& viscosityProperties,
    const volVectorField& U,
    const surfaceScalarField& phi
)
:
    viscosityModel(name, viscosityProperties, U, phi),
    BinghamPapanastasiouCoeffs_
    (
        viscosityProperties.optionalSubDict(typeName + "Coeffs")
    ),
    k_("k", dimViscosity, BinghamPapanastasiouCoeffs_),
    tau0_("tau0", dimViscosity/dimTime, BinghamPapanastasiouCoeffs_),
    eps_("eps", dimless/dimTime, BinghamPapanastasiouCoeffs_),
    nu_
    (
        IOobject
        (
            name,
            U_.time().timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        calcNu()
    )
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::viscosityModels::BinghamPapanastasiou::read
(
    const dictionary& viscosityProperties
)
{
    viscosityModel::read(viscosityProperties);

    BinghamPapanastasiouCoeffs_ =
        viscosityProperties.optionalSubDict(typeName + "Coeffs");

    BinghamPapanastasiouCoeffs_.lookup("k") >> k_;
    BinghamPapanastasiouCoeffs_.lookup("tau0") >> tau0_;
    BinghamPapanastasiouCoeffs_.lookup("eps") >> eps_;

    return true;
}


// ************************************************************************* //
