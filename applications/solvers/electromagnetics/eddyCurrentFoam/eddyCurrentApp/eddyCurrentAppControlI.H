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

inline bool
Foam::eddyCurrentApp::Control::meshIs3D() const
{
    return meshIs3D_;
};


inline const Foam::Vector<Foam::label>&
Foam::eddyCurrentApp::Control::solutionDir() const
{
    return solutionDir_;
};


inline const Foam::word&
Foam::eddyCurrentApp::Control::interfaceName() const
{
    return interfacePatchName_;
};


inline Foam::label
Foam::eddyCurrentApp::Control::interfaceLabel() const
{
    return interfacePatchLabel_;
};


inline bool
Foam::eddyCurrentApp::Control::switchV()
{
    switchV_ = !switchV_;

    return switchV_;
};


inline Foam::label
Foam::eddyCurrentApp::Control::subCorr() const
{
    return subCorr_;
}


#ifdef eddyCurrentAppLink_H

inline bool
Foam::eddyCurrentApp::Control::updateZeroCounter() const
{
    return (emUpdateData_.counter == 0);
}


inline bool
Foam::eddyCurrentApp::Control::updateOutputTimeIndex() const
{
    return (emUpdateData_.outputTimeIndex > 0)
        && (emUpdateData_.outputTimeIndex % emSettings_.outputTimeIndexCycle == 0);
}


inline bool
Foam::eddyCurrentApp::Control::updateTimeIndex() const
{
    return (mesh_.time().timeIndex() > 0)
        && (mesh_.time().timeIndex() % emSettings_.timeIndexCycle == 0);
}


inline bool
Foam::eddyCurrentApp::Control::updateTime() const
{
    return (mesh_.time().value() - emUpdateData_.lastTime > emSettings_.timeCycle);
}

#endif


// ************************************************************************* //
