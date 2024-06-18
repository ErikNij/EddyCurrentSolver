/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.0
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline label Foam::pimpleApp::Control::pRefCell() const
{
    return pRefCell_;
}

inline const scalar& Foam::pimpleApp::Control::pRefValue() const
{
    return pRefValue_;
}

inline void Foam::pimpleApp::Control::setpRefCell(const volScalarField& p)
{
    setRefCell(p, dict(), pRefCell_, pRefValue_);
}

inline void Foam::pimpleApp::Control::setpRefCell(label i, const scalar& v)
{
    pRefCell_ = i;
    pRefValue_ = v;
}

inline bool Foam::pimpleApp::Control::firstInnerIter() const
{
    return
       corrPISO_ == 1
    && corrNonOrtho_ == 1;
}

// ************************************************************************* //
