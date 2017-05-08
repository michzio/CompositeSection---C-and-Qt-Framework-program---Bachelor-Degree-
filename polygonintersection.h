#ifndef POLYGONINTERSECTION_H
#define POLYGONINTERSECTION_H

#include "gpc.h"
#include <QVarLengthArray>
#include "surface.h"

class PolygonIntersection
{
public:
    static QVarLengthArray<Surface *> *intersect(Surface *A, Surface *B);
private:
    static gpc_polygon* convertToGPCPolygon(Surface *);
    static Surface* extractFromGPCContour(const gpc_vertex_list& contour,
                                          std::string name, Material *material);
};

#endif // POLYGONINTERSECTION_H
