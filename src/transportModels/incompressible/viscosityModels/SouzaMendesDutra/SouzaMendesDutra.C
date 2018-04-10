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

#include "SouzaMendesDutra.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace viscosityModels
{
    defineTypeNameAndDebug(SouzaMendesDutra, 0);

    addToRunTimeSelectionTable
    (
        viscosityModel,
        SouzaMendesDutra,
        dictionary
    );
}
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::viscosityModels::SouzaMendesDutra::calcNu() const
{
    dimensionedScalar tone("tone", dimTime, 1.0);
    dimensionedScalar rtone("rtone", dimless/dimTime, 1.0);

    tmp<volScalarField> sr(strainRate());

    return
    (
		(tau0_ + k_ * rtone * pow(tone * sr(), n_)) * (1.0 - exp(-nu0_ * sr() / tau0_))
		/ (max(sr(), dimensionedScalar ("VSMALL", dimless/dimTime, VSMALL)))
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::viscosityModels::SouzaMendesDutra::SouzaMendesDutra
(
    const word& name,
    const dictionary& viscosityProperties,
    const volVectorField& U,
    const surfaceScalarField& phi
)
:
    viscosityModel(name, viscosityProperties, U, phi),
    SouzaMendesDutraCoeffs_
    (
        viscosityProperties.optionalSubDict(typeName + "Coeffs")
    ),
    k_("k", dimViscosity, SouzaMendesDutraCoeffs_),
    n_("n", dimless, SouzaMendesDutraCoeffs_),
    tau0_("tau0", dimViscosity/dimTime, SouzaMendesDutraCoeffs_),
    nu0_("nu0", dimViscosity, SouzaMendesDutraCoeffs_),
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

bool Foam::viscosityModels::SouzaMendesDutra::read
(
    const dictionary& viscosityProperties
)
{
    viscosityModel::read(viscosityProperties);

    SouzaMendesDutraCoeffs_ =
        viscosityProperties.optionalSubDict(typeName + "Coeffs");

    SouzaMendesDutraCoeffs_.lookup("k") >> k_;
    SouzaMendesDutraCoeffs_.lookup("n") >> n_;
    SouzaMendesDutraCoeffs_.lookup("tau0") >> tau0_;
    SouzaMendesDutraCoeffs_.lookup("nu0") >> nu0_;

    return true;
}


// ************************************************************************* //
