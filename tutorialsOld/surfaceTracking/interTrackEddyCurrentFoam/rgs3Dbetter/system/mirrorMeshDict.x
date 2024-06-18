/*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | foam-extend: Open Source CFD                    |
|  \\    /   O peration     | Version:     3.1                                |
|   \\  /    A nd           | Web:         http://www.extend-project.de       |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      mirrorMeshDict;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

/* A plane can be defined in different ways
 * 1. planeEquation
 *     planeType    planeEquation;
 *     planeEquationDict
 *     {
 *        a  1.0;
 *        b  2.0;
 *        c  3.0;
 *        d  0.0;
 *     }
 *
 * 2. embeddedPoints
 *     planeType    planeEquation;
 *     planeEquationDict
 *     {
 *        point1  (0 1 0);
 *        point2  (1 0 0);
 *        point3  (0 0 1);
 *     }
 *
 * 3. pointAndNormal
 *     planeType           pointAndNormal;
 *     pointAndNormalDict
 *     {
 *        basePoint       (0 0 0);
 *        normalVector    (0 1 0);
 *     }
 */

planeType           pointAndNormal;

pointAndNormalDict
{
    basePoint       (0 0 0);
    normalVector    (1 0 0);
}

planeTolerance      1e-6;

// ************************************************************************* //
