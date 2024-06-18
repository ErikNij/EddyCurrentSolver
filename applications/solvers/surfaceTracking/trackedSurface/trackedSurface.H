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

Class
    trackedSurface

Description
    Implementation of Surface (interface) tracking method
    based on moving mesh.

SourceFiles
    trackedSurface.C
    trackedSurfaceDemandDrivenData.C
    trackedSurfaceMoveFixedPatches.C
    trackedSurfacePointDisplacement.C
    trackedSurfaceWrite.C

Authors
    Zejlko Tukovic, Pascal Beckstein

\*---------------------------------------------------------------------------*/

#ifndef trackedSurface_H
#define trackedSurface_H


#include "IOdictionary.H"
#include "fvCFD.H"
#include "wallFvPatch.H"
#include "faCFD.H"
#include "emptyFaPatch.H"
#include "wedgeFaPatch.H"
#include "wedgeFaPatchFields.H"
#include "slipFaPatchFields.H"
#include "fixedValueFaPatchFields.H"
#include "zeroGradientFaPatchFields.H"
#include "faSubMesh.H"
#include "IOpatchToPatchInterpolation.H"
#include "dynamicFvMesh.H"
#include "twoPhaseMixture.H"
#include "turbulenceModel.H"
#include "surfactantProperties.H"
#include "demandDrivenData.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class trackedSurface Declaration
\*---------------------------------------------------------------------------*/

class trackedSurface
:
    public IOdictionary
{
    // Private data

        //- Prefix for naming
        word prefix_;

        //- Upper case prefix for naming
        word Prefix_;

// TEST: Volume conservation
        scalar vol0_;

        //- Reference to fvMesh
        dynamicFvMesh& mesh_;

        //- Reference to density field
        const volScalarField& rho_;

        //- Reference to velocity field
        volVectorField& U_;

        //- Reference to pressure field
        volScalarField& p_;

        //- Reference to fluid flux field
        surfaceScalarField& phi_;

        //- Pointer to temperature field
        volScalarField* TPtr_;

        //- Pointer to concentration field
        volScalarField* cPtr_;

        //- Pointer to gravitational acceleration
        const uniformDimensionedVectorField* gPtr_;

        //- Pointer to the transportModel
        const twoPhaseMixture* transportPtr_;

        //- Pointer to turbulence model
        const incompressible::turbulenceModel* turbulencePtr_;

        //- Pointer to external pressure field
        const volScalarField* p0Ptr_;

        //- Current time index for surface
        label curTimeIndex_;

        //- Has the mesh moved?
        bool meshMoved_;

        //- Are there exist two fluids in the surface model
        Switch twoFluids_;

        //- Is it surface points displacement direction
        //  parallel with surface point normals
        Switch normalMotionDir_;

        //- Surface points displacement direction
        //  if not normal motion direction
        vector motionDir_;

        //- Is it surface clean of surfactant
        Switch cleanInterface_;

        //- ID of the fvMesh patch which represents
        //  A side of the surface
        label aPatchID_;

        //- ID of the fvMesh patch which represents
        //  B side of the surface
        label bPatchID_;

        //- Viscosity of fluid A
        dimensionedScalar muFluidA_;

        //- Viscosity of fluid B
        dimensionedScalar muFluidB_;

        //- Density of fluid A
        dimensionedScalar rhoFluidA_;

        //- Density of fluid B
        dimensionedScalar rhoFluidB_;

        //- Heat conductivity of fluid A
        dimensionedScalar kFluidA_;

        //- Heat conductivity of fluid B
        dimensionedScalar kFluidB_;

        //- Heat capacity of fluid A
        dimensionedScalar CpFluidA_;

        //- Heat capacity of fluid B
        dimensionedScalar CpFluidB_;

        //- Reference to gravitational acceleration
        dimensionedVector g_;

        //- Surface tension for the clean surface
        dimensionedScalar cleanInterfaceSurfTension_;

        //- Surface patches which do not move
        wordList fixedTrackedSurfacePatches_;

        //- Surface patches for witch point normals
        //  must be corrected
        wordList pointNormalsCorrectionPatches_;

         //- Number of surface correctors
        label nTrackedSurfCorr_;

// TEST: Fixed interface
        //- Fixed interface
        Switch fixedInterface_;

// TEST: Sub-mesh
        //- Use sub-mesh
        Switch useSubMesh_;

// TEST: Volume conservation
        //- Apply correction to flux to conserve total volume
        Switch correctVolume_;

// TEST: Move always from start
        Switch total0Update_;

        //- Reset control points at the begining of time step
        Switch resetControlPoints_;

        //- Interface smoothing at the begining of time step
        Switch smoothing_;

        //- Free contact angle
        Switch freeContactAngle_;

        //- Correct dispacement for points next to fixed points
        Switch correctPointNormals_;

        //- Correct dispacement for points next to fixed points
        Switch correctDisplacement_;

        //- Correct curvature for faces next to fixed points
        Switch correctCurvature_;

        //- Smooth curvature
        Switch smoothCurvature_;

        //- Extrapolation order for curvature correction
        label curvExtrapOrder_;

        //- Calc normal derivative of normal velocity using fvc
        Switch fvcNGradUn_;

        //- Implicit coupling
        Switch implicitCoupling_;

        //- Relative interface deformation limit
        scalar interfaceDeformationLimit_;

    // Demand-driven data

        //- Patch to patch interpolation object which deals with
        //  interpolation of flow variables between two surface
        //  (interface) patches
        mutable IOpatchToPatchInterpolation* interpolatorABPtr_;

        //- Patch to patch interpolation object which deals with
        //  interpolation of flow variables between two surface
        //  (interface) patches
        mutable IOpatchToPatchInterpolation* interpolatorBAPtr_;

        //- Points which are attached to the surface A side faces
        //  and which defines the surface shape
        mutable vectorIOField* controlPointsPtr_;

        //- Field which additionally determines
        //  the motion of surface points
        mutable labelList* motionPointsMaskPtr_;

        //- Displacement direction of surface points
        mutable vectorField* pointsDisplacementDirPtr_;

        //- Displacement direction of surface control points
        mutable vectorField* facesDisplacementDirPtr_;

        //- Total displacement of surface points
        //  in one time step
        mutable vectorIOField* totalDisplacementPtr_;

// TEST: Move always from start
        mutable pointField* points0Ptr_;

// TEST: Move always from start
        mutable vectorIOField* total0DisplacementPtr_;

        //- Finite area mesh which consists
        //  surface side A faces
        mutable faMesh* aMeshPtr_;

// TEST: Sub-mesh
        //- Finite area sub-mesh
        mutable faSubMesh* aSubMeshPtr_;

        //- Surface velocity field
        mutable areaVectorField* UsPtr_;

        //- Surface fluid flux
        mutable edgeScalarField* phisPtr_;

        //- Surface surfactant concetration
        mutable areaScalarField* surfactConcPtr_;

        //- Surface tension field
        mutable areaScalarField* surfaceTensionPtr_;

        //- Surfactant properties
        mutable surfactantProperties* surfactantPtr_;

        //- Fluid indicator
        mutable volScalarField* fluidIndicatorPtr_;

        //- Effiective viscosity value field of fluid A
        mutable scalarField* muEffFluidAvalPtr_;

        //- Effiective viscosity value field of fluid B
        mutable scalarField* muEffFluidBvalPtr_;

        //- Contact angle
        mutable edgeScalarField* contactAnglePtr_;

        //- Surface temperature
        mutable areaScalarField* temperaturePtr_;

        //- Surface concentration
        mutable areaScalarField* concentrationPtr_;

        //- Surface tension force
        mutable areaVectorField* surfaceTensionForcePtr_;

        //- Normal derivative of normal velocity at the interface
        mutable scalarField* nGradUnPtr_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        trackedSurface(const trackedSurface&);

        //- Disallow default bitwise assignment
        void operator=(const trackedSurface&);

        //- Demand-driven data

            void makeInterpolators();
            void makeControlPoints() const;
            void makeMotionPointsMask();
            void makeDirections();
            void makeTotalDisplacement();
// TEST: Move always from start
            void makeTotal0Displacement();
            void makeFaMesh() const;
// TEST: Sub-mesh
            void makeFaSubMesh() const;
            void makeUs() const;
            void makePhis();
            void makePhi();
            void makeDdtPhi();
            void makeSurfactConc() const;
            void makeSurfaceTension() const;
            void makeSurfactant() const;
            void makeFluidIndicator();
            void makeMuEffFluidAval() const;
            void makeMuEffFluidBval() const;
            void makeContactAngle() const;
            void makeTemperature() const;
            void makeConcentration() const;
            void makeSurfaceTensionForce() const;
            void makeNGradUn() const;

            //- Clear all demand-driven data
            void clearOut();

        //- Initialization

            //- Init properties
            void initProperties();

            //- Init/Check point normals correction
            void initCheckPointNormalsCorrection();

            //- Init/Check surface patches
            void initCheckSurfacePatches();

            //- Init motion point mask
            void initMotionPointMask();

            //- Init control points position
            void initControlPointsPosition();

        //- General operations

            //- Reset control points
            void resetControlPoints(bool force = false);

            //- Update displacement directions
            void adjustDisplacementDirections();

        //- Move calculations

// TEST: Volume conservation
            //- Correct interface phi for domain volume conservation
            void correctInterfacePhi(scalarField& interfacePhi);

            //- Calculate swept volume correction for given interface phi
            tmp<scalarField> calcSweptVolCorr(const scalarField& interfacePhi);

            //- Calculate shift of interface control points from swept volume field
            tmp<scalarField> calcDeltaH(const scalarField& sweptVolCorr);

            //- Calculate and add correction for interface control points shift
            void calcAddDeltaHcorrection(scalarField& deltaH);

            //- Calculate displacement
            tmp<vectorField> calcDisplacement(const scalarField& interfacePhi);

            //- Calculate new mesh points from interface displacement field
            tmp<pointField> calcNewMeshPoints(const vectorField& displacement);

            //- Calculate new mesh points from interface phi
            tmp<pointField> calcNewMeshPoints(const scalarField& interfacePhi);

            //- Calculate new mesh points from negative total displacement field
            tmp<pointField> calcNewMeshPoints();

        //- Move operations

            //- Move only surface points
            void movePoints(const scalarField& interfacePhi);

            //- Move and update all mesh points with motion solver for
            //  total displacement and delete total displacement pointer
            void moveMeshPoints();

            //- Move and update all mesh points with motion solver for
            //  given displacement field
            void moveMeshPoints(const pointField& displacement);

            //- Move finite area sub-mesh
            void moveFaSubMesh();

            //- Move sub-mesh of all correctedFvPatchField of one field
            template<class Type>
            void
            moveCorrectedPatchSubMesh
            (
                GeometricField<Type, fvPatchField, volMesh>& vf
            ) const;

            //- Move sub-meshes of all relevant fields
            void moveCorrectedPatchSubMeshes();

            //- Update all dependencies after point movement
            void moveUpdate();

        //- Point displacment calculation

            //- Move control ponts for deltaH and calculate
            //  surface points displacement
            //  for given new control points position
            tmp<vectorField> pointDisplacement(const scalarField& deltaH);

            //- Calc least sqare plane point and normal
            tmp<vectorField> lsPlanePointAndNormal
            (
                const vectorField& points,
                const vector& origin,
                const vector& axis
            ) const;

            //- Point displacement correction
            //  (currently only for fixed patches)
            void correctPointDisplacement
            (
                const scalarField& sweptVolCorr,
                vectorField& displacement
            );

//             //- Test inverse distance point displacement
//             tmp<vectorField> pointDisplacementID(const scalarField& deltaH);


// TODO
        void correctCurvature();
        void correctPointNormals();
        void smoothCurvature();
        void correctContactLinePointNormals();

public:

    // Runtime type information
    TypeName("trackedSurface");


    // Constructors

        trackedSurface
        (
            dynamicFvMesh& m,
            const volScalarField& rho,
            volVectorField& Ub,
            volScalarField& Pb,
            surfaceScalarField& phi,
            volScalarField* TbPtr = NULL,
            volScalarField* cbPtr = NULL,
            const uniformDimensionedVectorField* gPtr = NULL,
            const twoPhaseMixture* transportPtr = NULL,
            const incompressible::turbulenceModel* turbulencePtr = NULL,
            const volScalarField* p0Ptr = NULL,
            word prefix = typeName
        );


    // Destructor

        ~trackedSurface();


    // Member Functions

        //- Database

            //- Return reference to DB
            const Time& DB() const
            {
                return U_.mesh().time();
            }

            //- Return reference to fvMesh
            dynamicFvMesh& mesh()
            {
                return mesh_;
            }

            //- Return reference to fvMesh
            const dynamicFvMesh& mesh() const
            {
                return mesh_;
            }

        //- Related fields

            //- Return reference to density field
            const volScalarField& rho() const
            {
                return rho_;
            }

            //- Return reference to velocity field
            volVectorField& U()
            {
                return U_;
            }

            //- Return reference to velocity field
            const volVectorField& U() const
            {
                return U_;
            }

            //- Return reference to pressure field
            volScalarField& p()
            {
                return p_;
            }

            //- Return reference to pressure field
            const volScalarField& p() const
            {
                return p_;
            }

            //- Return reference to flux field
            surfaceScalarField& phi()
            {
                return phi_;
            }

            //- Return reference to flux field
            const surfaceScalarField& phi() const
            {
                return phi_;
            }

            //- Return reference to temperature field
            volScalarField& T()
            {
                return *TPtr_;
            }

            //- Return reference to temperature field
            const volScalarField& T() const
            {
                return *TPtr_;
            }

            //- Return reference to concentration field
            volScalarField& c()
            {
                return *cPtr_;
            }

            //- Return reference to concentration field
            const volScalarField& c() const
            {
                return *cPtr_;
            }

            //- Return reference to transport model
            const twoPhaseMixture& transport() const
            {
                return *transportPtr_;
            }

            //- Return reference to turbulence model auto pointer
            const incompressible::turbulenceModel& turbulence() const
            {
                return *turbulencePtr_;
            }

            //- Return reference to external pressure field
            const volScalarField& p0() const
            {
                return *p0Ptr_;
            }

        //- Has the mesh moved?
        bool meshMoved() const
        {
            return meshMoved_;
        }

        //- Are there exist two fluid
        //  in the surface model
        bool twoFluids() const
        {
            return twoFluids_;
        }

        //- Motion direction swithc
        bool normalMotionDir() const
        {
            return normalMotionDir_;
        }

        //- Is it interface clean
        bool cleanInterface() const
        {
            return cleanInterface_;
        }

        //- ID of patch which represent A side
        //  of surface
        const label& aPatchID() const
        {
            return aPatchID_;
        }

        //- ID of patch which represent B side
        //  of surface
        const label& bPatchID() const
        {
            return bPatchID_;
        }

        //- Viscosity of fluid A
        const dimensionedScalar& muFluidA() const
        {
            return muFluidA_;
        }

        //- Viscosity of fluid B
        const dimensionedScalar& muFluidB() const
        {
            return muFluidB_;
        }

        //- Density of fluid A
        const dimensionedScalar& rhoFluidA() const
        {
            return rhoFluidA_;
        }

        //- Densiti of fluid B
        const dimensionedScalar& rhoFluidB() const
        {
            return rhoFluidB_;
        }

        //- Heat conductivity of fluid A
        const dimensionedScalar& kFluidA() const
        {
            return kFluidA_;
        }

        //- Heat conductivity of fluid B
        const dimensionedScalar& kFluidB() const
        {
            return kFluidB_;
        }

        //- Heat capacity of fluid A
        const dimensionedScalar& CpFluidA() const
        {
            return CpFluidA_;
        }

        //- Heat capacity of fluid B
        const dimensionedScalar& CpFluidB() const
        {
            return CpFluidB_;
        }

        //- Gravity
        const dimensionedVector& g() const
        {
            return g_;
        }

        //- Surface tension for clean interface
        const dimensionedScalar& cleanInterfaceSurfTension() const
        {
            return cleanInterfaceSurfTension_;
        }

        // Marangoni stress
        bool MarangoniStress() const
        {
            return (!cleanInterface_ || TPtr_ || cPtr_);
        }

        //- Is the interface fixed?
        bool fixedInterface() const
        {
            return fixedInterface_;
        }

        //- Return reference to interpolator from surface side A to
        //  surface side B
        const IOpatchToPatchInterpolation& interpolatorAB();

        //- Return reference to interpolator from surface side B to
        //  surface side A
        const IOpatchToPatchInterpolation& interpolatorBA();

        //- Return control points
        vectorField& controlPoints();

        //- Return control points
        const vectorField& controlPoints() const;

        //- Return reference to motion points mask field
        labelList& motionPointsMask();

        //- Return reference to point displacement direction field
        vectorField& pointsDisplacementDir();

        //- Return reference to control points displacement direction field
        vectorField& facesDisplacementDir();

        //- Return total points displacement
        vectorField& totalDisplacement();

        //- Return points0
        vectorField& points0()
        {
            return *points0Ptr_;
        }

        //- Return total0 points displacement
        vectorField& total0Displacement();

        //- Return reference to finite area mesh
        faMesh& aMesh();

        //- Return reference to finite area mesh
        const faMesh& aMesh() const;

// TEST: Sub-mesh
        //- Return reference to sub-mesh
        faSubMesh& aSubMesh();

// TEST: Sub-mesh
        //- Return reference to sub-mesh
        const faSubMesh& aSubMesh() const;

// TEST: Sub-mesh
        //- Return reference to curvature
        const areaScalarField& curvature() const
        {
            if(aSubMeshPtr_)
            {
                return aSubMesh().faceCurvatures();
            }
            else
            {
                return aMesh().faceCurvatures();
            }
        }

        //- Return surface velocity field
        areaVectorField& Us();

        //- Return surface velocity field
        const areaVectorField& Us() const;

        //- Return surface fluid flux field
        edgeScalarField& Phis();

        //- Return surface surfactant concentration field
        areaScalarField& surfactantConcentration();

        //- Return surface surfactant concentration field
        const areaScalarField& surfactantConcentration() const;

        //- Return surface tension field
        areaScalarField& surfaceTension();

        //- Return surface tension field
        const areaScalarField& surfaceTension() const;

        //- Return surfactant properties
        const surfactantProperties& surfactant() const;

        //- Return reference to fluid indicator
        const volScalarField& fluidIndicator();

        //- Return effiective viscosity value field of fluid A
        scalarField& muEffFluidAval();

        //- Return effiective viscosity value field of fluid A
        const scalarField& muEffFluidAval() const;

        //- Return effiective viscosity value field of fluid B
        scalarField& muEffFluidBval();

        //- Return effiective viscosity value field of fluid B
        const scalarField& muEffFluidBval() const;

        //- Return reference to contact angle
        edgeScalarField& contactAngle();

        //- Return reference to contact angle
        const edgeScalarField& contactAngle() const;

        //- Return reference to interface temperature
        areaScalarField& temperature();

        //- Return reference to interface temperature
        const areaScalarField& temperature() const;

        //- Return reference to interface concentration
        areaScalarField& concentration();

        //- Return reference to interface concentration
        const areaScalarField& concentration() const;

        //- Return reference to surface tension force
        const areaVectorField& surfaceTensionForce() const;

        //- Return normal derivative of normal velicty at the interface
        scalarField& nGradUn();

        //- Return normal derivative of normal velicty at the interface
        const scalarField& nGradUn() const;

        //- Update properties
        void updateProperties();

        //- Move all mesh points from total displacement
        bool updateMesh();

        //- Update interface points position
        void updatePoints();

        //- Correct fixed surface patches
        void moveFixedPatches(const vectorField& displacement);

        //- Smoothing interface
        void smoothing();

        //- Smooth surface mesh using Laplacian smoothing
        void smoothMesh();

        //- Update surface boundary conditions
        void updateBoundaryConditions();

        //- Update contact angle
        void updateContactAngle();

        //- Correct effective kinematic viscosity on surface
        void updateMuEff();

        //- Update surface flux
        void updateSurfaceFlux();

        //- Update surface surfactant concentration
        void updateSurfactantConcentration();

        //- Update surface temperature
        void updateTemperature();

        //- Update surface concentration
        void updateConcentration();

        //- Correct surface velocity
        void updateVelocity();

        //- Correct surface pressure
        void updatePressure();

        //- Correct normal derivative of normal velocity
        void updateNGradUn();

        //- Calculate pressure force
        vectorField pressureForce() const;

        //- Calculate viscous force
        vectorField viscousForce() const;

        //- Calculate total pressure force
        vector totalPressureForce() const;

        //- Calculate total viscous force
        vector totalViscousForce() const;

        //- Calculate normal viscous force
        vector normalViscousForce() const;

        //- Calculate tangential viscous force
        vector tangentialViscousForce() const;

        //- Calculate total normal force
        vector totalNormalForce() const;

        //- Calculate total tangential force
        vector totalTangentialForce() const;

        //- Calculate total force
        vector totalForce() const;

        //- Calculate total normal surface tension force
        vector totalNormalSurfaceTensionForce() const;

        //- Calculate total tangential surface tension force
        vector totalTangentialSurfaceTensionForce() const;

        //- Calculate total surface tension force
        vector totalSurfaceTensionForce() const;

        //- Maximal surface tension based Courant number
        scalar maxCourantNumber();

        //- Writing

            //- Write VTK trackedSurface mesh
            void writeVTK() const;

            //- Write VTK point field
            void writeVTKpoints
            (
                const word fieldName,
                const vectorField& pf
            ) const;

            //- Write VTK trackedSurface control points
            void writeVTKControlPoints();


            //- Write area fields as boundary of corresponding volume fields
            template<class Type>
            void
            writeVol
            (
                const GeometricField<Type, faPatchField, areaMesh>& af
            ) const;

            //- Write selected fields of tracked surface side A
            void writeA();

            //- Write selected fields as volFields of tracked surface side A
            void writeVolA();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
