#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QGraphicsScene>
#include "sectionsolver.h"
#include "mouseeventshandler.h"
#include <QString>
#include <string>
#include <QDomElement>
#include "solutionwindow.h"
#include "materialwindow.h"


namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT
    
public:
    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();
    SectionSolver *getSectionSolver(void);

public slots:
    void calculatePushButtonClicked();

private slots:
    void on_addRowPushButton_clicked();

    void on_removeRowPushButton_clicked();

    void on_surfaceTableWidget_cellChanged(int row, int column);

    void on_isOpeningCheckBox_stateChanged(int arg1);

    void on_surfaceMaterialComboBox_currentIndexChanged(int index);

    void on_saveSurfacePushButton_clicked();

    void on_addRowPushButton_2_clicked();

    void on_removeRowPushButton_2_clicked();

    void on_lineTableWidget_cellChanged(int row, int column);

    void on_lineMaterialComboBox_currentIndexChanged(int index);

    void on_addRowPushButton_3_clicked();

    void on_removeRowPushButton_3_clicked();

    void on_fibergroupTableWidget_cellChanged(int row, int column);

    void on_saveFibergroupPushButton_clicked();

    void on_saveLinePushButton_clicked();

    void on_tabWidget_currentChanged(int index);

    void on_sectionDefinitionTableWidget_cellDoubleClicked(int row, int column);

    void on_loadSurfaceFromFileButton_clicked();

    void on_loadLineFromFileButton_clicked();

    void on_loadFiberGroupFromFileButton_clicked();

    void on_loadSectionFromFileButton_clicked();

    void on_pushButton_clicked();

    void on_pushButton_2_clicked();

private:
    //interfejs
    Ui::MainWindow *ui;
    QGraphicsScene *scene;
    SolutionWindow *solutionWindow;

    SectionSolver *sectionSolver;
    Surface *surfaceToEdit = NULL;
    Line *lineToEdit = NULL;
    FiberGroup *fibergroupToEdit = NULL;
    void showEmptyCellMessage(void);
    void switchToSurfaceEditionMode(Surface *);
    void switchToLineEditionMode(Line *);
    void switchToFiberGroupEditionMode(FiberGroup *);
    int fckFromString(QString string);
    GradeOfStructuralSteel structuralSteelGradeFromString(QString string);
    GradeOfSteel steelGradeFromString(QString string);
    FiberGroup *loadFiberGroupDefinitionFromXML(const QDomElement& fibergroupElement);
    Line *loadLineDefinitionFromXML(const QDomElement& lineElement);
    Surface *loadSurfaceDefinitionFromXML(const QDomElement& surfaceElement);

    void drawSurface(Surface *, Qt::GlobalColor);
};

#endif // MAINWINDOW_H
