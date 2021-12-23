#include <QCoreApplication>

#include "hacdHACD.h"
#include "cd_wavefront.h"

bool fun(const char *, double, double, size_t) {
    return true;
}

void test_HACD(const char* filename) {
    ConvexDecomposition::WavefrontObj wo;

    unsigned int tcount = wo.loadObj(filename);

    printf("old numTriangles= %d\n",wo.mTriCount);
    printf("old numIndices = %d\n",wo.mTriCount*3);
    printf("old numVertices = %d\n",wo.mVertexCount);

    if (tcount) {
        std::vector< HACD::Vec3<HACD::Real> > points;
        std::vector< HACD::Vec3<long> > triangles;

        for(int i=0; i<wo.mVertexCount; i++ )
        {
            int index = i*3;
            HACD::Vec3<HACD::Real> vertex(wo.mVertices[index], wo.mVertices[index+1],wo.mVertices[index+2]);
            points.push_back(vertex);
        }

        for(int i=0;i<wo.mTriCount;i++)
        {
            int index = i*3;
            HACD::Vec3<long> triangle(wo.mIndices[index], wo.mIndices[index+1], wo.mIndices[index+2]);
            triangles.push_back(triangle);
        }

        HACD::HACD myHACD;
        myHACD.SetPoints(&points[0]);
        myHACD.SetNPoints(points.size());
        myHACD.SetTriangles(&triangles[0]);
        myHACD.SetNTriangles(triangles.size());
        myHACD.SetCompacityWeight(0.1);
        myHACD.SetVolumeWeight(0.0);

        // HACD parameters
        // Recommended parameters: 2 100 0 0 0 0
        size_t nClusters = 2;
        double concavity = 100;
        bool invert = false;
        bool addExtraDistPoints = false;
        bool addNeighboursDistPoints = false;
        bool addFacesPoints = false;

        myHACD.SetNClusters(nClusters);                     // minimum number of clusters
        myHACD.SetNVerticesPerCH(100);                      // max of 100 vertices per convex-hull
        myHACD.SetConcavity(concavity);                     // maximum concavity
        myHACD.SetAddExtraDistPoints(addExtraDistPoints);
        myHACD.SetAddNeighboursDistPoints(addNeighboursDistPoints);
        myHACD.SetAddFacesPoints(addFacesPoints);

        myHACD.SetCallBack(fun);
        myHACD.Compute();
        nClusters = myHACD.GetNClusters();

        myHACD.Save("output.wrl", false);
    }
}

int main(int argc, char *argv[])
{
    test_HACD("./toyplane.obj");
    return 0;
}
