#include "mygraphicsview.h"

MyGraphicsView::MyGraphicsView(QWidget *parent) :
    QGraphicsView(parent)
{
}

void MyGraphicsView::mouseMoveEvent(QMouseEvent *event)
{
    fprintf(stderr, "Mouse move event in QGraphicsView...\n");

}

void MyGraphicsView::mousePressEvent(QMouseEvent *event)
{
    fprintf(stderr, "Mouse press event in QGraphicsView...\n");
    //double rad = 1;
    QPointF pt = mapToScene(event->pos());
    startX = pt.x(); startY = pt.y();
    //this->scene()->addEllipse(pt.x()-rad, pt.y()-rad, rad*2.0, rad*2.0,
                //QPen(), QBrush(Qt::SolidPattern));
}

void MyGraphicsView::mouseReleaseEvent(QMouseEvent *event)
{
    fprintf(stderr, "Moude release event in QGraphicsView...\n");
    QPointF pt = mapToScene(event->pos());
    double endX = pt.x();
    double endY = pt.y();

    if((endX-startX) != 0 || (endY-startY) != 0) {

        //draw line
        this->scene()->addLine(QLineF(startX, startY, endX, endY), QPen());
    }

}
