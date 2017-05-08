#include "gridgraphicsscene.h"
#include <QPainter>

inline qreal round(qreal val, int step) {
   int tmp = int(val) + step /2;
   tmp -= tmp % step;
   return qreal(tmp);
}

//constructor
GridGraphicsScene::GridGraphicsScene(QObject *parent ) : QGraphicsScene(parent)
{}

GridGraphicsScene::GridGraphicsScene(QObject *parent, const QString& xLabel, const QString& yLabel) : QGraphicsScene(parent)
{
    this->xLabel = xLabel;
    this->yLabel = yLabel;
}


void GridGraphicsScene::drawBackground(QPainter *painter, const QRectF &rect)
{

   fprintf(stderr, "RECT: %g, %g", rect.width(), rect.height());
   //Scale calculeted as current drawing width divided by initial width
   double scale = (int)rect.width()/600;
   scale = (scale == 0) ? 1 : scale;

   //setting painter Font size adjusted to size of drawing based on calculeted scale
   QFont font; //size of text labels
   font.setPixelSize((int) (12*scale));
   painter->setFont(font);
   //setting thickness of lines drawing by painter using pen
   QPen pen1(QColor(200, 200, 255, 125)); //pen for drawng grid
   pen1.setWidth((int) (1*scale));
   pen1.setCosmetic(false);
   QPen pen2(QColor(100, 100, 255, 125)); //pen for drawng text
   pen2.setWidth((int) (1*scale));
   pen2.setCosmetic(false);
   QPen pen3(QColor(0, 0, 0, 100));
   pen3.setWidth(1*scale);
   pen3.setCosmetic(false);

   int stepY = gridStepY*scale;
   int stepX = gridStepX*scale;
   painter->setPen(pen1);

   // draw horizontal grid
   qreal start = round(rect.top(), stepY);
   if (start > rect.top()) {
      start -= stepY;
   }
   painter->setPen(pen2);
   painter->drawText(-50*scale, rect.top()+20*scale, getYLabel());
   painter->setPen(pen1);
   for (qreal y = start - stepY; y < rect.bottom(); ) {
      y += stepY;
      painter->drawLine(rect.left(), y, rect.right(), y);
      painter->setPen(pen2);
      painter->drawText(10*scale, y-5*scale, QString::number(-y/gridScaleY));
      painter->setPen(pen1);
   }

   // now draw vertical grid
   start = round(rect.left(), stepX);
   if (start > rect.left()) {
      start -= stepX;
   }

   qreal x;
   for (x = start - stepX; x < rect.right(); ) {
      x += stepX;
      painter->drawLine(x, rect.top(), x, rect.bottom());
      painter->setPen(pen2);
      painter->drawText(x+10*scale,-5*scale, QString::number(x/gridScaleX));
      painter->setPen(pen1);
   }
   painter->setPen(pen2);
   painter->drawText(rect.right()-50*scale, 30*scale, getXLabel());
   painter->setPen(pen1);

   //drawing coordinate system axies
   painter->setPen(pen3);
   qreal xLeft = rect.left(); qreal xRight = rect.right();
   painter->drawLine(xLeft, 0, xRight, 0);

   qreal yBottom = rect.bottom(); qreal yTop = rect.top();
   painter->drawLine(0, yBottom, 0, yTop);
}

void GridGraphicsScene::setXLabel(const QString& xLabel)
{
    this->xLabel = xLabel;
}

void GridGraphicsScene::setYLabel(const QString &yLabel)
{
    this->yLabel = yLabel;
}

const QString& GridGraphicsScene::getXLabel(void) const
{
    return this->xLabel;
}

const QString& GridGraphicsScene::getYLabel(void) const
{
    return this->yLabel;
}

void GridGraphicsScene::setGridStepX(int step)
{
    gridStepX = step;
}

void GridGraphicsScene::setGridStepY(int step)
{
    gridStepY = step;
}

void GridGraphicsScene::setGridScaleX(double scale)
{
    gridScaleX = scale;
}

void GridGraphicsScene::setGridScaleY(double scale)
{
    gridScaleY = scale;
}

