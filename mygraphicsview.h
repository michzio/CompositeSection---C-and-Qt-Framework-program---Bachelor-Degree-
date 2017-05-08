#ifndef MYGRAPHICSVIEW_H
#define MYGRAPHICSVIEW_H

#include <QGraphicsView>
#include <QMouseEvent>

class MyGraphicsView : public QGraphicsView
{
    Q_OBJECT

    double startX, startY;
public:
    explicit MyGraphicsView(QWidget *parent = 0);
    
signals:
    
public slots:

protected:
    void mouseMoveEvent(QMouseEvent *event);
    void mousePressEvent(QMouseEvent *event);
    void mouseReleaseEvent(QMouseEvent *event);
};

#endif // MYGRAPHICSVIEW_H
