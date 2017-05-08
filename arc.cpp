#include "arc.h"

void Arc::setStartPoint(Point *start) {
    this->startPoint = start;
}

void Arc::setEndPoint(Point *end) {
    this->endPoint = end;
}

void Arc::setAngle(double angle) {
    this->angle = angle;
}

Point *Arc::getStartPoint() {
    return startPoint;
}

Point *Arc::getEndPoint() {
    return endPoint;
}

double Arc::getAngle() {
    return angle;
}
