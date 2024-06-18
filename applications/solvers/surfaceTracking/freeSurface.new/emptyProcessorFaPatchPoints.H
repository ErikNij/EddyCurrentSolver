// Boundary processor patch points
{
    const labelList& curPointEdges = pointEdges[curPoint];

    label patchID = -1;
    label edgeID = -1;

    forAll (curPointEdges, edgeI)
    {
        label curEdge = curPointEdges[edgeI];

        if (edgeFaces[curEdge].size() == 1)
        {
            forAll(aMesh().boundary(), patchI)
            {
                const labelList& curEdges = aMesh().boundary()[patchI];

                forAll(curEdges, eI)
                {
                    if (curEdges[eI] == curEdge)
                    {
                        if 
                        (
                            aMesh().boundary()[patchI].type()
                         != processorFaPatch::typeName
                        )
                        {
                            patchID = patchI;
                            edgeID = curEdge;
                            break;
                        }
                    }
                }
            }
        }
    }

    if (patchID != -1)
    {
        vector nE = 
            pointNormals[edges[edgeID].start()]
          + pointNormals[edges[edgeID].end()];

        nE /= mag(nE);

        vector eP =
            controlPoints()[edgeFaces[edgeID][0]]
          - edges[edgeID].centre(points);

        vector mirrorCtrlPoint =
            edges[edgeID].centre(points)
          - ((I - 2.0*nE*nE)&eP);

        lsPoints[curPatchPoint].setSize(lsPoints[curPatchPoint].size() + 1);
        lsPoints[curPatchPoint][lsPoints[curPatchPoint].size() - 1] =
            mirrorCtrlPoint;
    }
}
