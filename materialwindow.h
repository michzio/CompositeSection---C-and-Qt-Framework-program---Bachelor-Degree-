#ifndef MATERIALWINDOW_H
#define MATERIALWINDOW_H

#include <QMainWindow>
#include "material.h"
#include "frp.h"
#include "reinforcement.h"
#include "gridgraphicsscene.h"

namespace Ui {
class MaterialWindow;
}

class MaterialWindow : public QMainWindow
{
    Q_OBJECT

    Material *material;

public:
    explicit MaterialWindow(QWidget *parent = 0);
    MaterialWindow(Material *m, QWidget *parent = 0);
    ~MaterialWindow();

    void setMaterial(Material *m);
    void graphicsViewFitScene(Qt::AspectRatioMode aspectRatioMode);

private:
    Ui::MaterialWindow *ui;
    QWidget *parentWindow;
    GridGraphicsScene *scene;

    void plotSigmaEpsilon(void);
};

#endif // MATERIALWINDOW_H
