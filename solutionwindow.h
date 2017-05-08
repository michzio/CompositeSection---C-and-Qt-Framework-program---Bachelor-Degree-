#ifndef SOLUTIONWINDOW_H
#define SOLUTIONWINDOW_H

#include <QMainWindow>
#include "gridgraphicsscene.h"


namespace Ui {
class SolutionWindow;
}

class SolutionWindow : public QMainWindow
{
    Q_OBJECT

    const int N_RANGE = 10;

public:
    explicit SolutionWindow(QWidget *parent = 0);
    ~SolutionWindow();

private slots:

    void on_pushButtonDrawForMx_clicked();

    void on_pushButtonDrawForMy_clicked();

    void on_pushButtonDrawForN_clicked();

    void on_pushButtonDrawForAlfa_clicked();

    void on_pushButtonExport_clicked();

    void on_pushButtonImport_clicked();

private:
    Ui::SolutionWindow *ui;
    QWidget *parentWindow;
    GridGraphicsScene *sceneMxMy;
    GridGraphicsScene *sceneNMy;
    GridGraphicsScene *sceneNMx;
    GridGraphicsScene *sceneNM;

    void drawInteractionCurveMxMyForN(int N = 0);
    void drawInteractionCurveMyNForMx(int Mx = 0);
    void drawInteractionCurveMxNForMy(int My = 0);
    void drawInteractionCurveMNForAlfa(int alfa = 0);
};

#endif // SOLUTIONWINDOW_H
