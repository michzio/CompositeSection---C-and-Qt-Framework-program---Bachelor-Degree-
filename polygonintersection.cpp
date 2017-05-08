#include "polygonintersection.h"
//zalaczenie bibliotki gpc (General Polygon Clipping)
#include "gpc.h"
#include "gpc.c"
#include "flags.h"
#include<sstream>

/****************************************************
 * intesect(Surface *, Surface *) is public static  *
 * method that enables intersecting two Surfaces    *
 * using gpc (General Polygon Clipping) library     *
 * function returns QVarLengthArray<Surface *>      *
 * because there can be generated >1 new surface as *
 * the result of intersection (concave polygons)    *
 ****************************************************/
QVarLengthArray<Surface *> *PolygonIntersection::intersect(Surface *A, Surface *B)
{

    gpc_polygon *subject, *clip, *result;
    //Converts Surfaces object to gpc_polygon structures
    subject = convertToGPCPolygon(A);
    clip = convertToGPCPolygon(B);
    result = new gpc_polygon;

    //perform polygon intersections
    if(DEBUG)
        fprintf(stderr, "performing gpc_polygon_clip()\n");
    gpc_polygon_clip(GPC_INT, subject, clip, result);

    QVarLengthArray<Surface *> *resultSurfaces =
            new QVarLengthArray<Surface *>();

    for(int i = 0; i < result->num_contours; i++)
    {
        std::stringstream stream;
        stream << i;
        Surface *newSurface = extractFromGPCContour(result->contour[i],
        std::string("intersection_") + A->getName() + "_" + B->getName() + "_" + stream.str(),
                                                      B->getMaterial());
        resultSurfaces->append( newSurface );
    }

    return resultSurfaces;
}

gpc_polygon* PolygonIntersection::convertToGPCPolygon(Surface *surface)
{
        gpc_polygon *polygon;

        polygon = new gpc_polygon;
        //polygon initialization to empty polygon
        polygon->num_contours = 0;
        polygon->hole = NULL;
        polygon->contour = NULL;

        //fprintf(stderr, "gpc_polygon has been created\n");

        //we are getting array of points
        Point **point = surface->getPointsArray();
        //get number of points
        int numberOfPoints = surface->numberOfPoints();

        //creating gpc_vertex_list which will be added to gpc_polygon
        //as new gpc_contour
        gpc_vertex_list vertexList;
        vertexList.num_vertices = numberOfPoints;
        //vertex array allocation
        if(vertexList.num_vertices > 0) {
            vertexList.vertex = (gpc_vertex*)
                    malloc(vertexList.num_vertices*sizeof(gpc_vertex));
            if(!vertexList.vertex) {
                fprintf(stderr, "gpc malloc failure: vertex creation\n" );
                exit(0);
            }
        } else {
            vertexList.vertex = NULL;
        }

        for(int i=0; i<vertexList.num_vertices; i++) {
            vertexList.vertex[i].x = point[i]->getX();
            vertexList.vertex[i].y = point[i]->getY();
        }

        gpc_add_contour(polygon,
                        &vertexList,
                        0);

        return polygon;
}

Surface *PolygonIntersection::extractFromGPCContour(const gpc_vertex_list &contour,
                                                    std::string name, Material *material)
{
    Surface *surface = new Surface(name, material);
    //it finds only intersection of polygons so angle will be always 0.0
    //new polygon created as the result of intersection will have straight edges
    double omega = 0.0;

    if(DEBUG)
        fprintf(stderr, "Intersection vertices: \n");

    //GPC library generate vertices in clockwise (CW) order, we change it to counter-clockwise (CCW)
    for(int i=contour.num_vertices -1; i>=0 ;i--)
    {
        if(DEBUG)
            fprintf(stderr, "{%g, %g}\n", contour.vertex[i].x, contour.vertex[i].y);
        surface->addPointAndAngel(contour.vertex[i].x, contour.vertex[i].y, omega);
    }

    return surface;
}
