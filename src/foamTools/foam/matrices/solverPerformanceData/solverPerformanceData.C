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

#include "solverPerformanceData.H"
#include "token.H"

// * * * * * * * * * * * * * * * Member operators  * * * * * * * * * * * * * //

bool Foam::solverPerformanceData::operator!=
(
    const solverPerformanceData& spd
) const
{
    return
    (
        solverName()      != spd.solverName()
     || fieldName()       != spd.fieldName()
     || initialResidual() != spd.initialResidual()
     || finalResidual()   != spd.finalResidual()
     || nIterations()     != spd.nIterations()
     || converged()       != spd.converged()
     || singular()        != spd.singular()
    );
}


// * * * * * * * * * * * * * * * Istream Operator  * * * * * * * * * * * * * //

Foam::Istream& Foam::operator>>
(
    Istream& is,
    solverPerformanceData& spd
)
{
    is.readBeginList("solverPerformanceData");
    is  >> spd.solverName_
        >> spd.fieldName_;

    List<scalarField*> residual(2);
    residual[0] = &spd.initialResidual_;
    residual[1] = &spd.finalResidual_;

    forAll (residual, i)
    {
        residual[i]->clear();

        token irt(is);
        if (irt != token::BEGIN_LIST)
        {
            scalar val = 0;

            if (irt.isLabel())
            {
                val = irt.labelToken();
            }
            else if (irt.isScalar())
            {
                val = irt.scalarToken();
            }
            else
            {
                // TODO: FatalError
            }

            residual[i]->setSize(1, val);
        }
        else
        {
            while (!is.eof())
            {
                is >> irt;

                if (irt != token::END_LIST)
                {
                    scalar val = 0;

                    if (irt.isLabel())
                    {
                        val = irt.labelToken();
                    }
                    else if (irt.isScalar())
                    {
                        val = irt.scalarToken();
                    }
                    else
                    {
                        // TODO: FatalError
                    }

                    residual[i]->setSize
                    (
                        residual[i]->size() + 1,
                        val
                    );
                }
                else
                {
                    break;
                }
            }
        }
    }

    is  >> spd.nIterations_
        >> spd.converged_
        >> spd.singular_;
    is.readEndList("solverPerformanceData");

    return is;
}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const solverPerformanceData& spd
)
{
    os  << token::BEGIN_LIST
        << spd.solverName_ << token::SPACE
        << spd.fieldName_ << token::SPACE;

    os  << token::BEGIN_LIST << token::SPACE;
    forAll (spd.initialResidual_, cmpt)
    {
        os << spd.initialResidual_[cmpt] << token::SPACE;
    }
    os  << token::END_LIST << token::SPACE;

    os  << token::BEGIN_LIST << token::SPACE;
    forAll (spd.finalResidual_, cmpt)
    {
        os << spd.finalResidual_[cmpt] << token::SPACE;
    }
    os  << token::END_LIST << token::SPACE;

    os  << spd.nIterations_ << token::SPACE
        << spd.converged_ << token::SPACE
        << spd.singular_ << token::SPACE
        << token::END_LIST;

    return os;
}


// ************************************************************************* //
