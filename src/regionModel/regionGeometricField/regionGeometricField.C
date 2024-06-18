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

#include "regionGeometricField.H"
#include "fvcExtrapolate.H"
#include "calculatedFvPatchFields.H"
#include "demandDrivenData.H"
#include "surfaceInterpolate.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// check mesh for two fields
#define checkField(gf1, gf2, op)                                    \
if ((gf1).mesh() != (gf2).mesh())                                   \
{                                                                   \
    FatalErrorIn("checkField(gf1, gf2, op)")                        \
        << "different mesh for fields "                             \
        << (gf1).name() << " and " << (gf2).name()                  \
        << " during operatrion " <<  op                             \
        << abort(FatalError);                                       \
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template
<
    class Type, template<class> class PatchField, class GeoMesh,
    class RegionGeoMesh
>
const Foam::IOobject
Foam::regionGeometricField<Type, PatchField, GeoMesh, RegionGeoMesh>::
regionIOobject
(
    label regionI,
    const IOobject& IOo
) const
{
    const objectRegistry& db = regionMesh_[regionI].thisDb();

    return IOobject
    (
        IOo.name(),
        IOo.time().timeName(),
        db,
        IOo.readOpt(),
        IOo.writeOpt()
    );
}


template
<
    class Type, template<class> class PatchField, class GeoMesh,
    class RegionGeoMesh
>
bool
Foam::regionGeometricField<Type, PatchField, GeoMesh, RegionGeoMesh>::
linkFieldPtr
(
    label regionI,
    const IOobject& IOo
) const
{
    const objectRegistry& dbI = regionMesh_[regionI].thisDb();

    if
    (
        dbI.foundObject<GeometricFieldType>(IOo.name())
    )
    {
// TODO: Add member function to objectRegistry to get write access!
        GeometricFieldType& fieldRefI =
            const_cast<GeometricFieldType&>
            (
                dbI.lookupObject<GeometricFieldType>(IOo.name())
            );

        fieldPtrs_[regionI] = &fieldRefI;

        fieldLinked_[regionI] = true;

        return true;
    }
    else
    {
        return false;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template
<
    class Type, template<class> class PatchField, class GeoMesh,
    class RegionGeoMesh
>
Foam::regionGeometricField<Type, PatchField, GeoMesh, RegionGeoMesh>::
regionGeometricField
(
    const IOobject& IOo,
    const RegionMesh& regionMesh,
    const dimensioned<Type>& dim,
    const HashTable<IOobject> IOoOverride
)
:
    regIOobject
    (
        IOobject
        (
            IOo.name(),
            IOo.time().timeName(),
            regionMesh.thisDb(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            true
        )
    ),
    regionMesh_(regionMesh),
    fieldPtrs_(regionMesh.size(), NULL),
    fieldActive_(regionMesh.size(), true),
    fieldLinked_(regionMesh.size(), false)
{
    if (debug)
    {
        Info<< "Foam::regionGeometricField::regionGeometricField(...) : "
            << "Create region field " << IOo.name()
            << endl;
    }

    forAll (*this, regionI)
    {
        IOobject regionIOoI = regionIOobject(regionI, IOo);

        if (IOoOverride.found(regions()[regionI]))
        {
            regionIOoI = IOoOverride[regions()[regionI]];
        }

// TODO: Add member functions linkOrNew for GeometricField
        if (linkFieldPtr(regionI, regionIOoI))
        {
            if (debug)
            {
                Info<< "Foam::regionGeometricField::regionGeometricField(...) : "
                    << "Linked existing field " << regionIOoI.name() << " for region "
                    << regionMesh_.regions()[regionI]
                    << endl;
            }
        }
        else
        {
            if (debug)
            {
                Info<< "Foam::regionGeometricField::regionGeometricField(...) : "
                    << "Create new field " << regionIOoI.name() << " for region "
                    << regionMesh_.regions()[regionI]
                    << endl;
            }

            if (regionIOoI.readOpt() == IOobject::MUST_READ)
            {
                fieldPtrs_[regionI] =
                    new GeometricFieldType
                    (
                        regionIOoI,
                        regionMesh_[regionI]
                    );
            }
            else
            {
                fieldPtrs_[regionI] =
                    new GeometricFieldType
                    (
                        regionIOoI,
                        regionMesh_[regionI],
                        dim
                    );
            }
        }
    }
}


template
<
    class Type, template<class> class PatchField, class GeoMesh,
    class RegionGeoMesh
>
Foam::regionGeometricField<Type, PatchField, GeoMesh, RegionGeoMesh>::
regionGeometricField
(
    const IOobject& IOo,
    const RegionMesh& regionMesh,
    const dimensioned<Type>& dim,
    const word& patchFieldType,
    const HashTable<IOobject> IOoOverride
)
:
    regIOobject
    (
        IOobject
        (
            IOo.name(),
            IOo.time().timeName(),
            regionMesh.thisDb(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            true
        )
    ),
    regionMesh_(regionMesh),
    fieldPtrs_(regionMesh.size(), NULL),
    fieldActive_(regionMesh.size(), true),
    fieldLinked_(regionMesh.size(), false)
{
    if (debug)
    {
        Info<< "Foam::regionGeometricField::regionGeometricField(...) : "
            << "Create region field " << IOo.name()
            << endl;
    }

    forAll (*this, regionI)
    {
        IOobject regionIOoI = regionIOobject(regionI, IOo);

        if (IOoOverride.found(regions()[regionI]))
        {
            regionIOoI = IOoOverride[regions()[regionI]];
        }

// TODO: Add member functions linkOrNew for GeometricField
        if (linkFieldPtr(regionI, regionIOoI))
        {
            if (debug)
            {
                Info<< "Foam::regionGeometricField::regionGeometricField(...) : "
                    << "Linked existing field " << regionIOoI.name() << " for region "
                    << regionMesh_.regions()[regionI]
                    << endl;
            }
        }
        else
        {
            if (debug)
            {
                Info<< "Foam::regionGeometricField::regionGeometricField(...) : "
                    << "Create new field " << regionIOoI.name() << " for region "
                    << regionMesh_.regions()[regionI]
                    << endl;
            }

            if (regionIOoI.readOpt() == IOobject::MUST_READ)
            {
                fieldPtrs_[regionI] =
                    new GeometricFieldType
                    (
                        regionIOoI,
                        regionMesh_[regionI]
                    );
            }
            else
            {
                fieldPtrs_[regionI] =
                    new GeometricFieldType
                    (
                        regionIOoI,
                        regionMesh_[regionI],
                        dim,
                        patchFieldType
                    );
            }
        }
    }
}


template
<
    class Type, template<class> class PatchField, class GeoMesh,
    class RegionGeoMesh
>
Foam::regionGeometricField<Type, PatchField, GeoMesh, RegionGeoMesh>::
regionGeometricField
(
    const IOobject& IOo,
    const RegionMesh& regionMesh,
    const dimensioned<Type>& dim,
    const wordList& patchFieldTypes,
    const HashTable<IOobject> IOoOverride
)
:
    regIOobject
    (
        IOobject
        (
            IOo.name(),
            IOo.time().timeName(),
            regionMesh.thisDb(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            true
        )
    ),
    regionMesh_(regionMesh),
    fieldPtrs_(regionMesh.size(), NULL),
    fieldActive_(regionMesh.size(), true),
    fieldLinked_(regionMesh.size(), false)
{
    if (debug)
    {
        Info<< "Foam::regionGeometricField::regionGeometricField(...) : "
            << "Create region field " << IOo.name()
            << endl;
    }

    forAll (*this, regionI)
    {
        IOobject regionIOoI = regionIOobject(regionI, IOo);

        if (IOoOverride.found(regions()[regionI]))
        {
            regionIOoI = IOoOverride[regions()[regionI]];
        }

// TODO: Add member functions linkOrNew for GeometricField
        if (linkFieldPtr(regionI, regionIOoI))
        {
            if (debug)
            {
                Info<< "Foam::regionGeometricField::regionGeometricField(...) : "
                    << "Linked existing field " << regionIOoI.name() << " for region "
                    << regionMesh_.regions()[regionI]
                    << endl;
            }
        }
        else
        {
            if (debug)
            {
                Info<< "Foam::regionGeometricField::regionGeometricField(...) : "
                    << "Create field " << regionIOoI.name() << " for region "
                    << regionMesh_.regions()[regionI]
                    << endl;
            }

            if (regionIOoI.readOpt() == IOobject::MUST_READ)
            {
                fieldPtrs_[regionI] =
                    new GeometricFieldType
                    (
                        regionIOoI,
                        regionMesh_[regionI]
                    );
            }
            else
            {
                fieldPtrs_[regionI] =
                    new GeometricFieldType
                    (
                        regionIOoI,
                        regionMesh_[regionI],
                        dim,
                        patchFieldTypes
                    );
            }
        }
    }
}


template
<
    class Type, template<class> class PatchField, class GeoMesh,
    class RegionGeoMesh
>
Foam::regionGeometricField<Type, PatchField, GeoMesh, RegionGeoMesh>::
regionGeometricField
(
    const IOobject& IOo,
    const RegionGeometricFieldType& rgf,
    const HashTable<IOobject> IOoOverride
)
:
    regIOobject
    (
        IOobject
        (
            IOo.name(),
            IOo.time().timeName(),
            rgf.mesh().thisDb(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            true
        )
    ),
    regionMesh_(rgf.mesh()),
    fieldPtrs_(rgf.mesh().size(), NULL),
    fieldActive_(rgf.mesh().size(), true),
    fieldLinked_(rgf.mesh().size(), false)
{
    if (debug)
    {
        Info<< "Foam::regionGeometricField::regionGeometricField(...) : "
            << "Create region field " << IOo.name()
            << endl;
    }

    forAll (*this, regionI)
    {
        IOobject regionIOoI = regionIOobject(regionI, IOo);

        if (IOoOverride.found(regions()[regionI]))
        {
            regionIOoI = IOoOverride[regions()[regionI]];
        }

// TODO: Add member functions linkOrNew for GeometricField
        if (linkFieldPtr(regionI, regionIOoI))
        {
            if (debug)
            {
                Info<< "Foam::regionGeometricField::regionGeometricField(...) : "
                    << "Linked existing field " << regionIOoI.name() << " for region "
                    << regionMesh_.regions()[regionI]
                    << endl;
            }
        }
        else
        {
            if (debug)
            {
                Info<< "Foam::regionGeometricField::regionGeometricField(...) : "
                    << "Create field " << regionIOoI.name() << " for region "
                    << regionMesh_.regions()[regionI]
                    << endl;
            }

            if (regionIOoI.readOpt() == IOobject::MUST_READ)
            {
                fieldPtrs_[regionI] =
                    new GeometricFieldType
                    (
                        regionIOoI,
                        regionMesh_[regionI]
                    );
            }
            else
            {
                fieldPtrs_[regionI] =
                    new GeometricFieldType
                    (
                        regionIOoI,
                        rgf.field(regionI)
                    );
            }
        }
    }
}

// * * * * * * * * * * * * * Constructor link-wrapper  * * * * * * * * * * * //

template
<
    class Type, template<class> class PatchField, class GeoMesh,
    class RegionGeoMesh
>
Foam::regionGeometricField<Type, PatchField, GeoMesh, RegionGeoMesh>*
Foam::regionGeometricField<Type, PatchField, GeoMesh, RegionGeoMesh>::LinkOrNew
(
    const IOobject& IOo,
    const RegionMesh& regionMesh,
    const dimensioned<Type>& dim,
    const HashTable<IOobject> IOoOverride
)
{
    const objectRegistry& db = regionMesh.thisDb();

    RegionGeometricFieldType* ptr = NULL;

    if
    (
        db.foundObject<RegionGeometricFieldType>(IOo.name())
    )
    {
// TODO: Add member function to objectRegistry to get write access!
        RegionGeometricFieldType& ref =
            const_cast<RegionGeometricFieldType&>
            (
                db.lookupObject<RegionGeometricFieldType>(IOo.name())
            );

        ptr = &ref;
    }
    else
    {
        ptr = new RegionGeometricFieldType
            (
                IOo,
                regionMesh,
                dim,
                IOoOverride
            );
    }

    return ptr;
}


template
<
    class Type, template<class> class PatchField, class GeoMesh,
    class RegionGeoMesh
>
Foam::regionGeometricField<Type, PatchField, GeoMesh, RegionGeoMesh>*
Foam::regionGeometricField<Type, PatchField, GeoMesh, RegionGeoMesh>::LinkOrNew
(
    const IOobject& IOo,
    const RegionMesh& regionMesh,
    const dimensioned<Type>& dim,
    const word& patchFieldType,
    const HashTable<IOobject> IOoOverride
)
{
    const objectRegistry& db = regionMesh.thisDb();

    RegionGeometricFieldType* ptr = NULL;

    if
    (
        db.foundObject<RegionGeometricFieldType>(IOo.name())
    )
    {
// TODO: Add member function to objectRegistry to get write access!
        RegionGeometricFieldType& ref =
            const_cast<RegionGeometricFieldType&>
            (
                db.lookupObject<RegionGeometricFieldType>(IOo.name())
            );

        ptr = &ref;
    }
    else
    {
        ptr = new RegionGeometricFieldType
            (
                IOo,
                regionMesh,
                dim,
                patchFieldType,
                IOoOverride
            );
    }

    return ptr;
}


template
<
    class Type, template<class> class PatchField, class GeoMesh,
    class RegionGeoMesh
>
Foam::regionGeometricField<Type, PatchField, GeoMesh, RegionGeoMesh>*
Foam::regionGeometricField<Type, PatchField, GeoMesh, RegionGeoMesh>::LinkOrNew
(
    const IOobject& IOo,
    const RegionMesh& regionMesh,
    const dimensioned<Type>& dim,
    const wordList& patchFieldTypes,
    const HashTable<IOobject> IOoOverride
)
{
    const objectRegistry& db = regionMesh.thisDb();

    RegionGeometricFieldType* ptr = NULL;

    if
    (
        db.foundObject<RegionGeometricFieldType>(IOo.name())
    )
    {
// TODO: Add member function to objectRegistry to get write access!
        RegionGeometricFieldType& ref =
            const_cast<RegionGeometricFieldType&>
            (
                db.lookupObject<RegionGeometricFieldType>(IOo.name())
            );

        ptr = &ref;
    }
    else
    {
        ptr = new RegionGeometricFieldType
            (
                IOo,
                regionMesh,
                dim,
                patchFieldTypes,
                IOoOverride
            );
    }

    return ptr;
}


template
<
    class Type, template<class> class PatchField, class GeoMesh,
    class RegionGeoMesh
>
Foam::regionGeometricField<Type, PatchField, GeoMesh, RegionGeoMesh>*
Foam::regionGeometricField<Type, PatchField, GeoMesh, RegionGeoMesh>::LinkOrNew
(
    const IOobject& IOo,
    const RegionGeometricFieldType& rgf,
    const HashTable<IOobject> IOoOverride
)
{
    const objectRegistry& db = rgf.mesh().thisDb();

    RegionGeometricFieldType* ptr = NULL;

    if
    (
        db.foundObject<RegionGeometricFieldType>(IOo.name())
    )
    {
// TODO: Add member function to objectRegistry to get write access!
        RegionGeometricFieldType& ref =
            const_cast<RegionGeometricFieldType&>
            (
                db.lookupObject<RegionGeometricFieldType>(IOo.name())
            );

        ptr = &ref;
    }
    else
    {
        ptr = new RegionGeometricFieldType
            (
                IOo,
                rgf,
                IOoOverride
            );
    }

    return ptr;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template
<
    class Type, template<class> class PatchField, class GeoMesh,
    class RegionGeoMesh
>
Foam::regionGeometricField<Type, PatchField, GeoMesh, RegionGeoMesh>::
~regionGeometricField()
{
    forAll (*this, regionI)
    {
        if (fieldActive_[regionI] && !fieldLinked_[regionI])
        {
            deleteDemandDrivenData(fieldPtrs_[regionI]);
        }
    }
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template
<
    class Type, template<class> class PatchField, class GeoMesh,
    class RegionGeoMesh
>
void
Foam::regionGeometricField<Type, PatchField, GeoMesh, RegionGeoMesh>::operator=
(
    const RegionGeometricFieldType& rgf
)
{
    if (this == &rgf)
    {
        FatalErrorIn
        (
            "Foam::regionGeometricField<>::operator=(...)"
        )   << "Attempted assignment to self"
            << abort(FatalError);
    }

    checkField(*this, rgf, "=");

    forAll (*this, regionI)
    {
        checkField(operator[](regionI), rgf.field(regionI), "=");

        // only equate field contents not ID

        operator[](regionI).dimensionedInternalField()
            = rgf.field(regionI).dimensionedInternalField();

        operator[](regionI).boundaryField()
            = rgf.field(regionI).boundaryField();
    }
}


template
<
    class Type, template<class> class PatchField, class GeoMesh,
    class RegionGeoMesh
>
void
Foam::regionGeometricField<Type, PatchField, GeoMesh, RegionGeoMesh>::operator=
(
    const tmp <RegionGeometricFieldType>& trgf
)
{
    if (this == &(trgf()))
    {
        FatalErrorIn
        (
            "Foam::regionGeometricField<>::operator=(...)"
        )   << "Attempted assignment to self"
            << abort(FatalError);
    }

    const RegionGeometricFieldType& rgf = trgf();

    checkField(*this, rgf, "=");

    forAll (*this, regionI)
    {
        checkField(operator[](regionI), rgf.field(regionI), "=");

        // only equate field contents not ID

        operator[](regionI).dimensions()
            = rgf.field(regionI).dimensions();

        operator[](regionI).internalField().transfer
        (
            const_cast<Field<Type>&>(rgf.field(regionI).internalField())
        );

        operator[](regionI).boundaryField()
            = rgf.field(regionI).boundaryField();
    }

    trgf.clear();
}


// * * * * * * * * * * * * * * * * Mapping * * * * * * * * * * * * * * * * * //

template
<
    class Type, template<class> class PatchField, class GeoMesh,
    class RegionGeoMesh
>
void
Foam::regionGeometricField<Type, PatchField, GeoMesh, RegionGeoMesh>::
mapBoundaryField
(
    label regionI
) const
{
    label regionI0 = regions()[polyMesh::defaultRegion];

    const GeometricFieldType& vf0 = operator[](regionI0);

    GeometricFieldType& vf = operator[](regionI);

    const polyBoundaryMesh& pbm0 = vf0.mesh().boundaryMesh();
    const polyBoundaryMesh& pbm = vf.mesh().boundaryMesh();

    forAll (pbm, patchI)
    {
        if
        (
            vf.boundaryField()[patchI].type()
         == calculatedFvPatchField<Type>::typeName
        )
        {
            const polyPatch& patch = pbm[patchI];

            label patchI0 = pbm0.findPatchID(patch.name());

            // Patch is present in regionI AND also in
            // default region. Mapping is done based on
            // direct addressing.
            if (patchI0 > -1)
            {
                const polyPatch& patch0 = pbm0[patchI0];
                const Field<Type>& patchField0 = vf0.boundaryField()[patchI0];
                label patchStart0 = patch0.start();
                label patchSize0 = patchField0.size();

                Field<Type>& patchField = vf.boundaryField()[patchI];
                label patchStart = patch.start();

                const labelList& fmap =
                    mesh().typeMap(addressingTypes::FACE, regionI, regionI0);

                forAll (patchField, facei)
                {
                    label faceI = patchStart + facei;
                    label faceI0 = fmap[faceI];

                    label facei0 = faceI0 - patchStart0;

                    if (facei0 > -1 && facei0 < patchSize0)
                    {
                        patchField[facei] = patchField0[facei0];
                    }
                    else
                    {
                        // TODO: Mapped from internal face. Do what?
                        //       one side of an intersecting patch?
                        // NOTE: This should already be avoided during
                        //       construction of regionPolyMesh
                        //       (or the first time faceMap is used) by
                        //       checking all boundary faces?
                    }
                }
            }
            // Patch is only present in regionI but NOT in
            // default region. Thus, all faces of this patch
            // correspond to interal faces of default region.
            // Interpolation is applied to get values for the
            // corresponding patch field.
            else if (patchI0 == -1)
            {
                const GeometricField<Type, fvsPatchField, surfaceMesh>
                    sf0 = fvc::interpolate(vf0);
                const Field<Type>& sf0In = sf0.internalField();

                Field<Type>& patchField = vf.boundaryField()[patchI];
                label patchStart = patch.start();

                const labelList& fmap =
                    mesh().typeMap(addressingTypes::FACE, regionI, regionI0);

                forAll (patchField, facei)
                {
                    label faceI = patchStart + facei;
                    label faceI0 = fmap[faceI];

                    patchField[facei] = sf0In[faceI0];
                }
            }
            else
            {
                // TODO: What if someone changes the name of
                //       one side of an intersecting patch?
                // NOTE: This should already be avoided during
                //       construction of regionPolyMesh
                //       (or the first time faceMap is used) by
                //       checking all boundary faces?
            }
        }
    }

    vf.correctBoundaryConditions();
};


// ************************************************************************* //
