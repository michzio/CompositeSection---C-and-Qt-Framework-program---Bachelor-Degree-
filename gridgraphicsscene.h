#ifndef GRIDGRAPHICSSCENE_H
#define GRIDGRAPHICSSCENE_H

#include <QGraphicsScene>

class GridGraphicsScene : public QGraphicsScene
{
   static const int GRID_STEP = 100;
   int gridStepX = GRID_STEP;
   int gridStepY = GRID_STEP;
   static constexpr double GRID_SCALE = 1.0;
   double gridScaleX = GRID_SCALE;
   double gridScaleY = GRID_SCALE;
   QString xLabel;
   QString yLabel;


public:
    GridGraphicsScene(QObject *parent);
    GridGraphicsScene(QObject *parent, const QString& xLabel, const QString& yLabel);

    void setXLabel(const QString& xLabel);
    void setYLabel(const QString& yLabel);
    const QString& getXLabel(void) const;
    const QString& getYLabel(void) const;

    void setGridStepX(int step);
    void setGridStepY(int step);
    void setGridScaleX(double scale);
    void setGridScaleY(double scale);

protected:
   void drawBackground(QPainter * painter, const QRectF & rect );
};

#endif // GRIDGRAPHICSSCENE_H
