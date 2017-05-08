#include "mainwindow.h"
#include "ui_mainwindow.h"
#include "stdio.h"
#include <QMessageBox>
#include <QGraphicsView>
#include <QGraphicsPolygonItem>
#include <QFileDialog>
#include <QFile>
#include <QDomDocument>
#include <QDomNodeList>
#include <QList>
#include "strainprofile.h"
#define _USE_MATH_DEFINES

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow),
    sectionSolver(new SectionSolver)
{
    ui->setupUi(this);
    QObject::connect(ui->calculatePushButton, SIGNAL(clicked()), this, SLOT(calculatePushButtonClicked()) );

    //Surface **surfaces = sectionSolver->getSurfaces();
    fprintf(stderr, "Liczba zdefiniowanych powierzchni: %d.\n", sectionSolver->getNumberOfSurfaces());

    //setting GraphicsScene for drawing on GraphicsView
    scene = new QGraphicsScene(this);
    ui->graphicsView->setScene(scene);

    //initially after application launch we set visible only combo box which enable selecting concrete f_ck
    ui->fckComboBox->setVisible(true);
    ui->steelGradeComboBox->setVisible(false);

}

void MainWindow::drawSurface(Surface *surface, Qt::GlobalColor fillColor = Qt::red)
{
    //declaring QPolygonF object which will be
    //added to QGraphicsScene
    QPolygonF polygon;
    Point **point = surface->getPointsArray();

    if(surface->getGraphicsItem() != NULL)
        scene->removeItem(surface->getGraphicsItem());

    fprintf(stderr, "Liczba punktow: %d\n", surface->numberOfPoints());
    //filling QPolygonF with QPointF
    for(int i=0; i<surface->numberOfPoints(); i++)
    {
        fprintf(stderr, "{%g,%g}\n",  point[i]->getX(),
                                      point[i]->getY());
        polygon.append(QPointF( static_cast<int>(point[i]->getX() + 0.5),
                                -static_cast<int>(point[i]->getY() + 0.5)));
    }

    QGraphicsPolygonItem* polygonItem =
            scene->addPolygon(polygon, QPen(Qt::black), QBrush(fillColor));
    surface->setGraphicsItem(polygonItem);
}

MainWindow::~MainWindow()
{
    delete ui;
    delete sectionSolver;
}

SectionSolver *MainWindow::getSectionSolver(void)
{
    return sectionSolver;
}

void MainWindow::calculatePushButtonClicked() {
    fprintf(stderr, "Oblicz\n");
    Surface **surfaces = sectionSolver->getSurfaces();
    for(int i=0; i < sectionSolver->getNumberOfSurfaces(); i++) {
        surfaces[i]->performAdaptiveLinearization();
        Point *centroid = surfaces[i]->getGeometricalCenter();
        surfaces[i]->drawGeometricalCenter();
        fprintf(stderr, "Geometrical Center: (%g, %g).\n",
                centroid->getX(), centroid->getY());
    }

   /* //testing surface intersection!!!
    QVarLengthArray<Surface *> *intersected;
    intersected = surfaces[0]->intersection(surfaces[1]);
    for(int i=0; i< intersected->length(); i++) {
        drawSurface(intersected->at(i), Qt::green);
    }*/

    sectionSolver->drawCentoidOnScene(this->scene);


    double A_ = sectionSolver->weightedArea();
    double Sx_ = sectionSolver->weightedFirstMomentOfArea_Sx();
    double Sy_ = sectionSolver->weightedFirstMomentOfArea_Sy();
    fprintf(stderr, "Weighted A = %g, Sx = %g, Sy = %g\n", A_, Sx_, Sy_);
    fprintf(stderr, "Geometrical center Cx: %g, Cy: %g\n", Sy_/A_, Sx_/A_ );


    double y_MAX = sectionSolver->yMaxInCoordinateSystem(sectionSolver->getGeometricalCenter(),
                                          (315.0/180.0) *M_PI);
    double y_MIN = sectionSolver->yMinInCoordinateSystem(sectionSolver->getGeometricalCenter(),
                                          (315.0/180.0) *M_PI);

   // fprintf(stderr, "y_MIN: %g, y_MAX: %g\n", y_MIN, y_MAX);

    /*
    Point *p = new Point(0,100);
    p->setOrigin(sectionSolver->getGeometricalCenter());
    p->setRotation(M_PI*90/180);
    fprintf(stderr, "px: %g, py: %g\n", p->getLocalX(), p->getLocalY());
    */

    //FAFITIS EXAMPLE - STRAIN PROFILE TEST
    /*double eps_low = (0.0035/27.460) * (y_MAX - y_MIN - 27.460);
    double rotationAngle = (315.0/180.0)*M_PI;
    StrainProfile *profile = new StrainProfile(sectionSolver, rotationAngle, eps_low, -0.0035);
    double fi = profile->getCurvature_fi();
    double d = profile->getNeutralAxis_d();
    double e_o = profile->getStrainAtOrigin_eps_o();
    double e_p = profile->getStrainAtLocalPoint(new Point(0,-9.972));
    fprintf(stderr, "fi: %g, d: %g, e_o: %g, e_p: %g\n",
            fi, d, e_o, e_p);
    double N = sectionSolver->summateInternalActions_N(profile);
    double localeMx = sectionSolver->summateInternalActions_Mx(profile);
    double localeMy = sectionSolver->summateInternalActions_My(profile);
    double Mx = Section::transferBackFromLocalMx(localeMx, localeMy, rotationAngle); //[Nm]
    double My = Section::transferBackFromLocalMy(localeMx, localeMy, rotationAngle); //[Nm]

    fprintf(stderr, "N = %g\n", N);
    fprintf(stderr, "locale Mx = %g\n", localeMx);
    fprintf(stderr, "locale My = %g\n", localeMy);
    fprintf(stderr, "Mx = %g\n", Mx);
    fprintf(stderr, "My = %g\n", My);
    */
    //FAFITIS EXAMPLE END

   /* StrainProfile *profile = new StrainProfile(sectionSolver,
                                               M_PI*45/180, 0.01, -0.003);
    double fi = profile->getCurvature_fi();
    double d = profile->getNeutralAxis_d();
    double e_o = profile->getStrainAtOrigin_eps_o();
    double e_p = profile->getStrainAtPoint(new Point(0,200));
    fprintf(stderr, "fi: %g, d: %g, e_o: %g, e_p: %g\n",
            fi, d, e_o, e_p);
    */
      /*
      StrainProfile *profile = new StrainProfile(sectionSolver,
                                                M_PI*90/180, 0.0035, -0.0035);
      Surface *s = sectionSolver->getSurfaces()[0];
      double area = s->internalAction_N(profile);
      double S_x = s->internalAction_Mx(profile);
      double S_y = s->internalAction_My(profile);
      fprintf(stderr, "Area = %g, Sx = %g, Sy = %g\n", area, S_x, S_y);
     */
       //double fi = profile->getCurvature_fi();
       //double d = profile->getNeutralAxis_d();
       //double e_o = profile->getStrainAtOrigin_eps_o();
      // double e_p = profile->getStrainAtPoint(new Point(0,200));


    ///TEST ONLY CODE - REMOVE AFTER TEST
    //testing Fafitis formula for Sx, Sy, Ix, Iy, Ixy
    //Surface *s = surfaces[0]; //getting concrete surface
    //double area = s->internalAction_N(new StrainProfile(sectionSolver, M_PI*90.0/180.0, 0.0035, -0.0035));
    //double s_x = s->internalAction_Mx(new StrainProfile(sectionSolver, 0.0, 0.0035, -0.0035));
    //double s_y = s->internalAction_My(new StrainProfile(sectionSolver, 0.0, 0.0035, -0.0035));
    //fprintf(stderr, "TESTED AREA: %g\n", area);
    //fprintf(stderr, "TESTED S_X: %g\n", s_x);
    //fprintf(stderr, "TESTED S_Y: %g\n", s_y);
    ///TEST ONLY CODE - REMOVE AFTER TEST


    //method solving section and printing calculated N, Mx, My for permissible Strain Profiles
    sectionSolver->solveSection();
    solutionWindow = new SolutionWindow(this);
    solutionWindow->show();

    FiberGroup **fibergroup = sectionSolver->getFiberGroups();
    for(int i=0; i < sectionSolver->getNumberOfFiberGroups(); i++) {
        fibergroup[i]->drawFiberGroup(fibergroup[i]->getMaterial()->materialColor());
    }

    Line **line = sectionSolver->getLines();
    for(int i=0; i < sectionSolver->getNumberOfLines(); i++) {
        line[i]->drawLine(line[i]->getMaterial()->materialColor());
    }

}

//on click "+" button in Surface Tab view we add one additional row at the bottom of the table widget
void MainWindow::on_addRowPushButton_clicked()
{
    ui->surfaceTableWidget->insertRow(ui->surfaceTableWidget->rowCount());
}

//on click "-" button in Surface Tab view we remove last row at the table widget
void MainWindow::on_removeRowPushButton_clicked()
{
    ui->surfaceTableWidget->removeRow(ui->surfaceTableWidget->rowCount()-1);
}

//cell value changed in Surface Table
void MainWindow::on_surfaceTableWidget_cellChanged(int row, int column)
{
    fprintf(stderr, "New value in cell: %d, %d.\n", row, column);
}

void MainWindow::on_isOpeningCheckBox_stateChanged(int arg1)
{
    fprintf(stderr, "Opening Check Boc state: %d.\n", arg1);
}

void MainWindow::on_surfaceMaterialComboBox_currentIndexChanged(int index)
{
    fprintf(stderr, "Surface Material selected: %d.\n", index);
    if(index == 0) {
        //selected material is CONCRETE
        //we show the user new comboBox to select f_ck of this concrete material
        ui->steelGradeComboBox->setVisible(false);
        ui->fckComboBox->setVisible(true);
    } else if(index == 1) {
        //selected material is STRUCTURAL STEEL
        //we show the user ne comboBox to select steel grade of this structural steel material
        ui->fckComboBox->setVisible(false);
        ui->steelGradeComboBox->setVisible(true);
    } else {
        //UNDEFINDED MATERIAL!
        //we hide all additional comboBoxes
        ui->fckComboBox->setVisible(false);
        ui->steelGradeComboBox->setVisible(false);
    }
}

void MainWindow::on_saveSurfacePushButton_clicked()
{
    Surface *surface;
    Material *material;

    fprintf(stderr, "Add/Save Surface Button Clicked.\n");
    QString surfaceName = ui->surfaceNameLineEdit->text();
    fprintf(stderr, "New surface name: %s.\n", surfaceName.toStdString().c_str());

    int materialIdx = ui->surfaceMaterialComboBox->currentIndex();

    if(materialIdx == 0) {
        //material defined by user in comboBox is CONCRETE
        //we create Concrete material object with default value of f_ck
        int concrete_fck = fckFromString( ui->fckComboBox->currentText() );
        material = new Concrete(concrete_fck);

    } else if (materialIdx == 1) {
        //material defined by user in comboBox is STRUCTURAL STEEL
        //we create StructuralSteel material object with default value of steel grade
        GradeOfStructuralSteel grade = structuralSteelGradeFromString( ui->steelGradeComboBox->currentText());
        material = new StructuralSteel(grade);
    } else if(materialIdx == 2) {
        //material defined by user in comboBox is FAFITIS CONCRETE
        //we create FafitisConcrete material object
        material = new FafitisConcrete();
    } else if(materialIdx == 3) {
        //material defined by user in comboBox is HOLLOW
        //thas mean it is openning, empty surface
        material = new Hollow();
    } else {
        throw new MaterialNotFoundException("Given Material doesn't exists");
    }

    if(surfaceToEdit != NULL) {
        surface = surfaceToEdit;
        surface->setName(surfaceName.toStdString());
        surface->setMaterial(material);
        surface->clearPointsAndAngles();
    } else {
        surface = new Surface(surfaceName.toStdString(), material);
        surface->setGraphicsScene(scene);
    }
    surface->setIsOpening((bool)ui->isOpeningCheckBox->checkState());

    double x = 0,y = 0;
    double omega = 0;

    for(int i = 0; i< ui->surfaceTableWidget->rowCount(); i++) {
        if(ui->surfaceTableWidget->item(i,0)) {
            x = ui->surfaceTableWidget->item(i,0)->text().toDouble();
        } else {
            showEmptyCellMessage(); return;
        }
        if(ui->surfaceTableWidget->item(i,1)) {
            y = ui->surfaceTableWidget->item(i,1)->text().toDouble();
        } else {
            showEmptyCellMessage(); return;
        }
        if(ui->surfaceTableWidget->item(i,2)) {
            omega = ui->surfaceTableWidget->item(i,2)->text().toDouble();
            //converting degrees to radians
            omega = (omega/360) * 2*M_PI;
        } else {
            showEmptyCellMessage(); return;
        }
        surface->addPointAndAngel(x,y,omega);
    }
    if(!surfaceToEdit)
        sectionSolver->addSurface(surface);

    surface->drawSurface(surface->getMaterial()->materialColor());
    ui->graphicsView->fitInView(scene->sceneRect().x() + 50, scene->sceneRect().y() + 50,
                    scene->sceneRect().width()+100, scene->sceneRect().height() + 100, Qt::KeepAspectRatio);
}

void MainWindow::showEmptyCellMessage() {
    fprintf(stderr, "Empty Cell");
    QMessageBox Msgbox;
    Msgbox.setText("There are empty cells.");
    Msgbox.exec();
}

//on click "+" button in Line Tab view we add one additional row at the bottom of the table widget
void MainWindow::on_addRowPushButton_2_clicked()
{
    ui->lineTableWidget->insertRow(ui->lineTableWidget->rowCount());
}

//on click "-" button in Line Tab view we remove last row at the table widget
void MainWindow::on_removeRowPushButton_2_clicked()
{
    ui->lineTableWidget->removeRow(ui->lineTableWidget->rowCount()-1);
}

void MainWindow::on_lineTableWidget_cellChanged(int row, int column)
{
     fprintf(stderr, "New value in cell: %d, %d.\n", row, column);
}

void MainWindow::on_lineMaterialComboBox_currentIndexChanged(int index)
{
    fprintf(stderr, "Line Material selected: %d.\n", index);
    if(index == 0) {
        //selected material is REINFORCEMENT
        //we show the user new comboBox to select f_ck of this concrete material
         ui->reinforcementGradeComboBox->setVisible(true);
    } else if(index == 1) {
        //selected material is FRP
         ui->reinforcementGradeComboBox->setVisible(false);
    } else {
        //UNDEFINDED MATERIAL!
        //we hide all additional comboBoxes
        ui->reinforcementGradeComboBox->setVisible(false);
    }
}

//on click "+" button in Fibergroup Tab view we add one additional row at the bottom of the table widget
void MainWindow::on_addRowPushButton_3_clicked()
{
     ui->fibergroupTableWidget->insertRow(ui->fibergroupTableWidget->rowCount());
}

//on click "-" button in Fibergroup Tab view we remove last row at the table widget
void MainWindow::on_removeRowPushButton_3_clicked()
{
    ui->fibergroupTableWidget->removeRow(ui->fibergroupTableWidget->rowCount()-1);
}

void MainWindow::on_fibergroupTableWidget_cellChanged(int row, int column)
{
     fprintf(stderr, "New value in cell: %d, %d.\n", row, column);
}

void MainWindow::on_saveFibergroupPushButton_clicked()
{
    FiberGroup *fibergroup;
    Material *material;

    fprintf(stderr, "Add/Save Fibergroup Button Clicked.\n");
    QString fibergroupName = ui->fibergroupLineEdit->text();

    fprintf(stderr, "New fibergroup name: %s.\n", fibergroupName.toStdString().c_str() );

    int materialIdx = ui->fibergroupMaterialComboBox->currentIndex();

    if(materialIdx == 0) {
        //material defined by user in comboBox is REINFORCEMENT BARs
        //we get grade of reinforcement steel defined by user in comboBox
        QString grade_label = ui->reinforcementOfFiberGroupsComboBox->currentText();
        GradeOfSteel grade = steelGradeFromString(grade_label);
        material = new Reinforcement(grade);
    } else if(materialIdx == 1) {
        //material defined by user in comboBox is FAFITIES_REINFORCEMENT
        //we create this material
        material = new FafitisReinforcement();
    } else {
        throw new MaterialNotFoundException("Given Material doesn't exists");
    }

    if(fibergroupToEdit != NULL) {
        fibergroup = fibergroupToEdit;
        fibergroup->setName(fibergroupName.toStdString());
        fibergroup->setMaterial(material);
        fibergroup->clearPointsAndAreas();
    } else {
        fibergroup = new FiberGroup(fibergroupName.toStdString(), material);
        fibergroup->setGraphicsScene(scene);
    }

    double x = 0, y = 0, area = 0;

    for(int i = 0; i< ui->fibergroupTableWidget->rowCount(); i++) {
        if(ui->fibergroupTableWidget->item(i,0)) {
            x = ui->fibergroupTableWidget->item(i,0)->text().toDouble();
        } else {
            showEmptyCellMessage(); return;
        }
        if(ui->fibergroupTableWidget->item(i,1)) {
            y = ui->fibergroupTableWidget->item(i,1)->text().toDouble();
        } else {
            showEmptyCellMessage(); return;
        }
        if(ui->fibergroupTableWidget->item(i,2)) {
            area = ui->fibergroupTableWidget->item(i,2)->text().toDouble();
        } else {
            showEmptyCellMessage(); return;
        }
        fibergroup->addPointAndArea(x,y,area);
    }

    if(!fibergroupToEdit)
        sectionSolver->addFiberGroup(fibergroup);

    fibergroup->drawFiberGroup(fibergroup->getMaterial()->materialColor());
}

void MainWindow::on_saveLinePushButton_clicked()
{
    Line *line;
    Material *material;

    fprintf(stderr, "Add/Save Line Button Clicked.\n");
    QString lineName = ui->lineNameLineEdit->text();

    fprintf(stderr, "New line name: %s.\n", lineName.toStdString().c_str() );

    int materialIdx = ui->lineMaterialComboBox->currentIndex();

    if(materialIdx == 0) {
        //material defined by user in comboBox is REINFORCEMENT
        //we get grade of reinforcement steel defined by user in comboBox
        QString grade_label = ui->reinforcementGradeComboBox->currentText();
        GradeOfSteel grade = steelGradeFromString(grade_label);
        material = new Reinforcement(grade);
    } else if(materialIdx == 1) {
        //material defined by user in comboBox is FIBER_REINFORCED_PLASTIC
        material = new FRP();
    } else {
        throw new MaterialNotFoundException("Given Material doesn't exists");
    }

    if(lineToEdit != NULL) {
        line = lineToEdit;
        line->setName(lineName.toStdString());
        line->setMaterial(material);
        line->clearPointsAndAreas();
    } else {
        line = new Line(lineName.toStdString(), material);
        line->setGraphicsScene(scene);
    }

    double x = 0, y = 0, area = 0;

    for(int i = 0; i< ui->lineTableWidget->rowCount(); i++) {
        if(ui->lineTableWidget->item(i,0)) {
            x = ui->lineTableWidget->item(i,0)->text().toDouble();
        } else {
            showEmptyCellMessage(); return;
        }
        if(ui->lineTableWidget->item(i,1)) {
            y = ui->lineTableWidget->item(i,1)->text().toDouble();
        } else {
            showEmptyCellMessage(); return;
        }
        if(ui->lineTableWidget->item(i,2)) {
            area = ui->lineTableWidget->item(i,2)->text().toDouble();
        } else {
            showEmptyCellMessage(); return;
        }
        line->addPointAndArea(x,y,area);
    }
    if(!lineToEdit)
        sectionSolver->addLine(line);
    line->drawLine(line->getMaterial()->materialColor());
}

void MainWindow::on_tabWidget_currentChanged(int index)
{
    fprintf(stderr, "Aktualna zakładka: %d.\n", index);
    if(index == 0) {
        //Now user selected Section Definition Tab
        while (ui->sectionDefinitionTableWidget->rowCount() > 0)
        {
            //we remove all rows in table widget
            ui->sectionDefinitionTableWidget->removeRow(0);
        }

        int numOfSurfaces = sectionSolver->getNumberOfSurfaces();
        Surface **surfaces = sectionSolver->getSurfaces();

        for(int i=0; i<numOfSurfaces; i++) {
            int rowIdx = ui->sectionDefinitionTableWidget->rowCount();
            ui->sectionDefinitionTableWidget->insertRow(rowIdx);
            ui->sectionDefinitionTableWidget->setItem(rowIdx, 0, new QTableWidgetItem(surfaces[i]->getName().c_str()));
            ui->sectionDefinitionTableWidget->setItem(rowIdx, 1, new QTableWidgetItem("Surface"));
            ui->sectionDefinitionTableWidget->setItem(rowIdx, 2,
                                            new QTableWidgetItem(surfaces[i]->getMaterial()->getMaterialName().c_str()));
        }

        int numOfLines = sectionSolver->getNumberOfLines();
        Line **lines = sectionSolver->getLines();

        for(int i=0; i<numOfLines; i++) {
            int rowIdx = ui->sectionDefinitionTableWidget->rowCount();
            ui->sectionDefinitionTableWidget->insertRow(rowIdx);
            ui->sectionDefinitionTableWidget->setItem(rowIdx, 0, new QTableWidgetItem(lines[i]->getName().c_str()));
            ui->sectionDefinitionTableWidget->setItem(rowIdx, 1, new QTableWidgetItem("Line"));
            ui->sectionDefinitionTableWidget->setItem(rowIdx, 2,
                                            new QTableWidgetItem(lines[i]->getMaterial()->getMaterialName().c_str()));
        }

        int numOfFiberGroups = sectionSolver->getNumberOfFiberGroups();
        FiberGroup **fiberGroups = sectionSolver->getFiberGroups();

        for(int i=0; i<numOfFiberGroups; i++) {
            int rowIdx = ui->sectionDefinitionTableWidget->rowCount();
            ui->sectionDefinitionTableWidget->insertRow(rowIdx);
            ui->sectionDefinitionTableWidget->setItem(rowIdx, 0, new QTableWidgetItem(fiberGroups[i]->getName().c_str()));
            ui->sectionDefinitionTableWidget->setItem(rowIdx, 1, new QTableWidgetItem("Fibergroup"));
            ui->sectionDefinitionTableWidget->setItem(rowIdx, 2,
                                            new QTableWidgetItem(fiberGroups[i]->getMaterial()->getMaterialName().c_str()));
        }

        ui->sectionDefinitionTableWidget->setEditTriggers(QAbstractItemView::NoEditTriggers);
    } else if(index ==1) {
        ui->saveSurfacePushButton->setText("Add Surface");
        surfaceToEdit = NULL;
    } else if(index ==2) {
        ui->saveLinePushButton->setText("Add Line");
        lineToEdit = NULL;
    } else if(index == 3) {
        ui->saveFibergroupPushButton->setText("Add Fibergroup");
        fibergroupToEdit = NULL;
    }
}

void MainWindow::on_sectionDefinitionTableWidget_cellDoubleClicked(int row, int column)
{
    fprintf(stderr, "Cell: %d clicked.\n", row);
    QString type = ui->sectionDefinitionTableWidget->item(row,1)->text();
    if(type == "Surface") {
        Surface **surfaces = sectionSolver->getSurfaces();
        switchToSurfaceEditionMode(surfaces[row]);
    } else if(type == "Line") {
        Line **lines = sectionSolver->getLines();
        switchToLineEditionMode(lines[row - sectionSolver->getNumberOfSurfaces()]);
    } else if(type == "Fibergroup") {
        FiberGroup **groups = sectionSolver->getFiberGroups();
        switchToFiberGroupEditionMode(groups[row - sectionSolver->getNumberOfSurfaces() - sectionSolver->getNumberOfLines()]);
    }
}

void MainWindow::switchToSurfaceEditionMode(Surface *surface)
{
    QTableWidget *table = ui->surfaceTableWidget;
    double *angles = surface->getAngels();
    Point **points = surface->getPointsArray();

    //we remove all existing rows in table widget with surface definition
    while (table->rowCount() > 0)
    {
        table->removeRow(0);
    }
   //adding rows with point (x,y) and angle omega values in surface table widget
    for(int i=0; i< surface->numberOfPoints(); i++) {
        int rowIdx = table->rowCount();
        table->insertRow(rowIdx);
        fprintf(stderr, "x: %g\n", points[i]->getX());
        table->setItem(rowIdx, 0, new QTableWidgetItem(QString::number(points[i]->getX())));
        table->setItem(rowIdx, 1, new QTableWidgetItem(QString::number(points[i]->getY())));
        double omega = (angles[i]*360) / (2*M_PI);
        table->setItem(rowIdx, 2, new QTableWidgetItem(QString::number(omega)));
    }
    ui->surfaceNameLineEdit->setText(QString(surface->getName().c_str()));
    switch(surface->getMaterial()->getMaterialType()) {
    case CONCRETE: {
        ui->surfaceMaterialComboBox->setCurrentIndex(0);
        ui->steelGradeComboBox->setVisible(false);
        //setting current surface material characteristic strength under compression f_ck to correct value in comboBox
        Concrete *concrete = static_cast<Concrete *>(surface->getMaterial());
        int f_ck = std::abs(concrete->getCharacteristicStrength_fck());
        std::string f_ck_label = std::to_string(f_ck) + " MPa";
        ui->fckComboBox->setCurrentText(QString(f_ck_label.c_str()));
        ui->fckComboBox->setVisible(true);
        break;
    }
    case STRUCTURAL_STEEL: {
        ui->surfaceMaterialComboBox->setCurrentIndex(1);
        ui->fckComboBox->setVisible(false);
        //selecting current surface material grade of steel to correct value in comboBox
        StructuralSteel *steel = static_cast<StructuralSteel *>(surface->getMaterial());
        std::string grade_label = steel->getGradeAsString();
        ui->steelGradeComboBox->setCurrentText(QString(grade_label.c_str()));
        ui->steelGradeComboBox->setVisible(true);
        break;
    }
    case FAFITIS_CONCRETE: {
        ui->surfaceMaterialComboBox->setCurrentIndex(2);
        ui->fckComboBox->setVisible(false);
        ui->steelGradeComboBox->setVisible(false);
        break;
    }
    case HOLLOW: {
        ui->surfaceMaterialComboBox->setCurrentIndex(3);
        ui->fckComboBox->setVisible(false);
        ui->steelGradeComboBox->setVisible(false);
        break;
    }
    default:
        ui->surfaceMaterialComboBox->setCurrentIndex(0);
        ui->steelGradeComboBox->setVisible(false);
        ui->fckComboBox->setCurrentIndex(0);
        ui->fckComboBox->setVisible(true);
    }
    ui->isOpeningCheckBox->setChecked(surface->getIsOpening());
    ui->tabWidget->setCurrentIndex(1);
    ui->saveSurfacePushButton->setText("Save Surface");
    surfaceToEdit = surface;
}

void MainWindow::switchToLineEditionMode(Line *line)
{
    QTableWidget *table = ui->lineTableWidget;
    double *areas = line->getAreas();
    Point **points = line->getPointsArray();

    while (table->rowCount() > 0)
    {
        table->removeRow(0);
    }

    for(int i=0; i< line->numberOfPoints(); i++) {
        int rowIdx = table->rowCount();
        table->insertRow(rowIdx);
        //fprintf(stderr, "x: %g\n", points[i]->getX());
        table->setItem(rowIdx, 0, new QTableWidgetItem(QString::number(points[i]->getX())));
        table->setItem(rowIdx, 1, new QTableWidgetItem(QString::number(points[i]->getY())));
        table->setItem(rowIdx, 2, new QTableWidgetItem(QString::number(areas[i])));
    }

    ui->lineNameLineEdit->setText(QString(line->getName().c_str()));
    switch(line->getMaterial()->getMaterialType()) {
    case REINFORCEMENT: {
        ui->lineMaterialComboBox->setCurrentIndex(0);
        //hiding other materials comboBoxes
        Reinforcement *reinforcement = static_cast<Reinforcement *>(line->getMaterial());
        std::string grade_label = reinforcement->getGradeAsString();
        ui->reinforcementGradeComboBox->setCurrentText(QString(grade_label.c_str()));
        ui->reinforcementGradeComboBox->setVisible(true);
        break;
    }
    case FIBER_REINFORCED_PLASTIC: {
        ui->lineMaterialComboBox->setCurrentIndex(1);
        ui->reinforcementGradeComboBox->setVisible(false);
        break;
    }
    default:
        ui->lineMaterialComboBox->setCurrentIndex(0);
    }

    ui->tabWidget->setCurrentIndex(2);

    ui->saveLinePushButton->setText("Save Line");
    lineToEdit = line;
}

void MainWindow::switchToFiberGroupEditionMode(FiberGroup *fibergroup)
{
    QTableWidget *table = ui->fibergroupTableWidget;
    double *areas = fibergroup->getAreas();
    Point **points = fibergroup->getPointsArray();

    while (table->rowCount() > 0)
    {
        table->removeRow(0);
    }

    for(int i=0; i< fibergroup->numberOfPoints(); i++) {
        int rowIdx = table->rowCount();
        table->insertRow(rowIdx);
        //fprintf(stderr, "x: %g\n", points[i]->getX());
        table->setItem(rowIdx, 0, new QTableWidgetItem(QString::number(points[i]->getX())));
        table->setItem(rowIdx, 1, new QTableWidgetItem(QString::number(points[i]->getY())));
        table->setItem(rowIdx, 2, new QTableWidgetItem(QString::number(areas[i])));
    }

    ui->fibergroupLineEdit->setText(QString(fibergroup->getName().c_str()));
    switch(fibergroup->getMaterial()->getMaterialType()) {
    case REINFORCEMENT: {
        ui->fibergroupMaterialComboBox->setCurrentIndex(0);
        //hiding other materials comboBoxes
        Reinforcement *reinforcement = static_cast<Reinforcement *>(fibergroup->getMaterial());
        std::string grade_label = reinforcement->getGradeAsString();
        ui->reinforcementOfFiberGroupsComboBox->setCurrentText(QString(grade_label.c_str()));
        ui->reinforcementOfFiberGroupsComboBox->setVisible(true);
        break;
    }
    case FAFITIS_REINFORCEMENT: {
        ui->fibergroupMaterialComboBox->setCurrentIndex(1);
        ui->reinforcementOfFiberGroupsComboBox->setVisible(false);
        break;
    }
    default:
        ui->fibergroupMaterialComboBox->setCurrentIndex(0);
    }

    ui->tabWidget->setCurrentIndex(3);

    ui->saveFibergroupPushButton->setText("Save Fibergroup");
    fibergroupToEdit = fibergroup;
}

void MainWindow::on_loadSurfaceFromFileButton_clicked()
{

    QString fileName = QFileDialog::getOpenFileName(this, tr("Open File"),
                                                     "", tr("Files (*.xml)"));
    fprintf(stderr, "Loading surface from file: %s\n", fileName.toStdString().c_str());

    //opening file with xml with surface definition
    QFile xmlFile(fileName);
    if(!xmlFile.open(QIODevice::ReadOnly | QIODevice::Text)) {
        QMessageBox::critical(this, "Load XML file with surface definition problem",
                              "Couldn't open " + fileName + " to load surface definition!",
                              QMessageBox::Ok);
        return;
    }

    //creating DOM document object from xml file
    QDomDocument xmlDocument;
    if (!xmlDocument.setContent(xmlFile.readAll()))
    {
        QMessageBox::critical(this, "Creating DOM Document object Error",
                              "Couldn't create DOM Document object from xml file!",
                              QMessageBox::Ok);
        return;
    }

    //accessing root <surface> element of xml file
    QDomElement surfaceElement;
    surfaceElement = xmlDocument.documentElement();

    Surface * surface = loadSurfaceDefinitionFromXML(surfaceElement);
    switchToSurfaceEditionMode(surface);
    sectionSolver->addSurface(surface);
    surface->setGraphicsScene(scene);
    surface->drawSurface(surface->getMaterial()->materialColor());

    ui->graphicsView->fitInView(scene->sceneRect().x() + 50, scene->sceneRect().y() + 50,
                    scene->sceneRect().width()+100, scene->sceneRect().height() + 100, Qt::KeepAspectRatio);
}

/**
 * @brief MainWindows::fckFromString
 * @param string - takes in string e.g. "12 MPa" from comboBox which defines f_ck value in MPa
 * @return integer representation of this string f_ck characteristic strength of concrete under compression in MPa
 */
int MainWindow::fckFromString(QString fck_string)
{
    int f_ck = fck_string.split(" MPa").at(0).toInt();
    fprintf(stderr, "Selected concrete f_ck value is: %d.\n", f_ck);
    return f_ck; //[MPa]
}

/**
 * @brief MainWindow::structuralSteelGradeFromString
 * @param grade_string - text from selected field of comboBox with steel grades
 * @return steel grade represented as enumeration type GradeOfStructuralSteel
 */
GradeOfStructuralSteel MainWindow::structuralSteelGradeFromString(QString grade_string)
{
    if(grade_string == "S235")
    {
        fprintf(stderr, "Selected grade of steel is: S235.\n");
        return S235;
    } else if(grade_string == "S275")
    {
        fprintf(stderr, "Selected grade of steel is: S275.\n");
        return S275;
    } else if(grade_string == "S355")
    {
        fprintf(stderr, "Selected grade of steel is: S355.\n");
        return S355;
    } else if(grade_string == "S450")
    {
        fprintf(stderr, "Selected grade of steel is: S450.\n");
        return S450;
    }

    fprintf(stderr, "Not found selected grade of steel, used S235 by default.\n");
    return S235;
}

/**
 * @brief MainWindow::steelGradeFromString
 * @param grade_string - text from selected field of comboBox with steel grades
 * @return steel grade represented as enumeration type GradeOfSteel
 */
GradeOfSteel MainWindow::steelGradeFromString(QString grade_string)
{
    if(grade_string == "St3SYb500") {
        fprintf(stderr, "Selected grade of reinforcement is: St3SYb500.\n");
        return St3SYb500;
    } else if(grade_string == "RB500W") {
        fprintf(stderr, "Selected grade of reinforcement is: RB500W.\n");
        return RB500W;
    } else if(grade_string == "Bst500S") {
        fprintf(stderr, "Selected grade of reinforcement is: Bst500S.\n");
        return Bst500S;
    } else if(grade_string == "B500SP") {
        fprintf(stderr, "Selected grade of reinforcement is: B500SP.\n");
        return B500SP;
    }

    fprintf(stderr, "Not found selected grade of reinforcement, used RB500W by default");
    return RB500W;
}

/**
 * @brief MainWindow::on_loadLineFromFileButton_clicked
 * Method enables loading Line Element from xml file
 */
void MainWindow::on_loadLineFromFileButton_clicked()
{
    QString fileName = QFileDialog::getOpenFileName(this, tr("Open File"),
                                                    "", tr("Files (*.xml)"));
    fprintf(stderr, "Loading line from file: %s\n", fileName.toStdString().c_str());

    //opening file with xml with line definition
    QFile xmlFile(fileName);
    if(!xmlFile.open(QIODevice::ReadOnly | QIODevice::Text)) {
        QMessageBox::critical(this, "Load XML file with line definition problem",
                              "Couldn't open " + fileName + " to load line definition!",
                              QMessageBox::Ok);
        return;
    }

    //creating DOM document object from xml file
    QDomDocument xmlDocument;
    if(!xmlDocument.setContent(xmlFile.readAll()))
    {
        QMessageBox::critical(this, "Creating DOM Document object Error",
                              "Couldn't create DOM Document object from xml file!",
                              QMessageBox::Ok);
        return;
    }

    //accessing root <line> element of xml file
    QDomElement lineElement;
    lineElement = xmlDocument.documentElement();

    Line *line = loadLineDefinitionFromXML(lineElement);
    switchToLineEditionMode(line);
    sectionSolver->addLine(line);
    line->setGraphicsScene(scene);
    line->drawLine(line->getMaterial()->materialColor());
    //line->drawLine();
}

/**
 * @brief MainWindow::on_loadFiberGroupFromFileButton_clicked
 * Method enables loading Fiber Group from xml file
 */
void MainWindow::on_loadFiberGroupFromFileButton_clicked()
{
    QString fileName = QFileDialog::getOpenFileName(this, tr("Open File"),
                                                    "", tr("Files (*.xml)"));
    fprintf(stderr, "Loading fiber group from file: %s\n", fileName.toStdString().c_str());

    //opening file with xml with fiber group definition
    QFile xmlFile(fileName);
    if(!xmlFile.open(QIODevice::ReadOnly | QIODevice::Text)) {
        QMessageBox::critical(this, "Load XML file with fiber group definition problem",
                              "Couldn't open " + fileName + " to load fiber group definition!");
        return;
    }

    //creating DOM document object from xml file
    QDomDocument xmlDocument;
    if(!xmlDocument.setContent(xmlFile.readAll()))
    {
        QMessageBox::critical(this, "Creating DOM Document object Error",
                              "Couldn't create DOM Document object from xml file!",
                              QMessageBox::Ok);
        return;
    }

    //accessing root <fibergroup> element of xml file
    QDomElement fibergroupElement;
    fibergroupElement = xmlDocument.documentElement();

    //using helper function to loadFiberGroupDefinitionFromXML(QDomElement)
    FiberGroup *fibergroup = loadFiberGroupDefinitionFromXML(fibergroupElement);
    switchToFiberGroupEditionMode(fibergroup);
    sectionSolver->addFiberGroup(fibergroup);
    //fibergroup->drawFibergoup();
    fibergroup->setGraphicsScene(scene);
    fibergroup->drawFiberGroup(fibergroup->getMaterial()->materialColor());
}

void MainWindow::on_loadSectionFromFileButton_clicked()
{
    QString fileName = QFileDialog::getOpenFileName(this, tr("Open File"),
                                                    "", tr("Files (*.xml)"));
    fprintf(stderr, "Loading section from file: %s\n", fileName.toStdString().c_str());

    //opening file with xml with section definition (that consist of surfaces, lines and fibergroups)
    QFile xmlFile(fileName);
    if(!xmlFile.open(QIODevice::ReadOnly | QIODevice::Text))
    {
        QMessageBox::critical(this, "Load XML file with section definition problem",
                              "Couldn't open " + fileName + " to load section definition!",
                              QMessageBox::Ok);

        return;
    }

    //creating DOM document object form xml file
    QDomDocument xmlDocument;
    if(!xmlDocument.setContent(xmlFile.readAll()))
    {
        QMessageBox::critical(this, "Creating DOM Document object Error",
                              "Couldn't create DOM Document object from xml file!",
                              QMessageBox::Ok);
        return;
    }
    //accessing root <section> element of xml file
    QDomElement sectionElement;
    sectionElement = xmlDocument.documentElement();

    //this root <section> element contains multiple
    //<surface>, <line>, <fibergroup> elements
    QDomNodeList surfaceList = sectionElement.elementsByTagName("surface");
    QDomNodeList lineList = sectionElement.elementsByTagName("line");
    QDomNodeList fibergroupList = sectionElement.elementsByTagName("fibergroup");

    //loops through each list and gathering information about each component
    for(int i=0; i < surfaceList.count(); i++)
    {
        QDomElement surfaceElement = surfaceList.at(i).toElement();
        Surface *surface = loadSurfaceDefinitionFromXML(surfaceElement);
        sectionSolver->addSurface(surface);
        surface->setGraphicsScene(scene);
        surface->drawSurface(surface->getMaterial()->materialColor());
        ui->graphicsView->fitInView(scene->sceneRect().x() + 50, scene->sceneRect().y() + 50,
                        scene->sceneRect().width()+100, scene->sceneRect().height() + 100, Qt::KeepAspectRatio);
    }

    for(int i=0; i<lineList.count(); i++)
    {
        QDomElement lineElement = lineList.at(i).toElement();
        Line *line = loadLineDefinitionFromXML(lineElement);
        sectionSolver->addLine(line);
        //line->drawLine();
        line->setGraphicsScene(scene);
        line->drawLine(line->getMaterial()->materialColor());
    }

    for(int i=0; i<fibergroupList.count(); i++)
    {
        QDomElement fibergroupElement = fibergroupList.at(i).toElement();
        FiberGroup *fibergroup = loadFiberGroupDefinitionFromXML(fibergroupElement);
        sectionSolver->addFiberGroup(fibergroup);
        //fibergroup->drawFiberGroup();
        fibergroup->setGraphicsScene(scene);
        fibergroup->drawFiberGroup(fibergroup->getMaterial()->materialColor());
    }

    //reloading section TableWidget
    on_tabWidget_currentChanged(0);
}

/**
 * @brief MainWindow::loadSurfaceDefinitionFromXML
 * @param surfaceElement - QDomElement pointing to node (element) containing definition of surface component
 * Method having parent node <surface> passed as function argument
 * reads into defintion of surface contained in its child nodes (elements)
 */
Surface *MainWindow::loadSurfaceDefinitionFromXML(const QDomElement &surfaceElement)
{
    //constructing Surface object
    Surface *surface;
    Material *material;

    QDomElement nameElement = surfaceElement.firstChildElement("name");
    fprintf(stderr, "Loaded surface name is: %s\n", nameElement.text().toStdString().c_str());

    QDomNodeList materialList = surfaceElement.elementsByTagName("material");
    if(materialList.count() !=1) {
        QMessageBox::critical(this, "XML file doesn't contain material definition",
                              "XML file doesn't contain material definiton!",
                              QMessageBox::Ok);
        return NULL;
    }
    //there should be only one <material> element in xml file
    QDomElement materialElement = materialList.at(0).toElement();
    fprintf(stderr, "Material: %s\n", materialElement.text().toStdString().c_str());

    //depending on material type: 'Concrete' or 'Structural Steel' we need to read fck or grade attribute
    if(materialElement.text() == "Concrete") {

        //we need to read fck attribute correspondig to Concrete material
        int fck = materialElement.attribute("fck", "12").toInt(); //default value of fck is minimal concrete strength under compression
        fprintf(stderr, "Concrete characteristic strength under compression fck = %d.\n", fck);
        material = new Concrete(fck);

    } else if(materialElement.text() == "Structural Steel")
    {
        //we need to read structural steel grade
        QString grade_label = materialElement.attribute("grade", "S235");
        GradeOfStructuralSteel grade = structuralSteelGradeFromString(grade_label);
        fprintf(stderr, "Structural Steel grade is: %s.\n", grade_label.toStdString().c_str());
        material = new StructuralSteel(grade);
    } else if(materialElement.text() == "Fafitis Concrete")
    {
        material = new FafitisConcrete();
    } else if(materialElement.text() == "Hollow" || materialElement.text() == "Empty"
              || materialElement.text() == "Openning") {
        material = new Hollow();
    }

    QDomElement openingElement = surfaceElement.firstChildElement("opening");
    int isOpening = openingElement.text().toInt(); //0 - isn't opening 1 - is opening
    fprintf(stderr, "Is surface an opening? %s\n", (isOpening == 1) ? "YES" : "NO");

    surface = new Surface(nameElement.text().toStdString(), material);
    surface->setIsOpening((bool) isOpening);

    QDomNodeList pointsList = surfaceElement.elementsByTagName("points").at(0).toElement().elementsByTagName("point");

    for(int i=0; i<pointsList.count(); i++)
    {
        QDomElement pointElement = pointsList.at(i).toElement();
        double x = pointElement.attribute("x", "0").toDouble();
        double y = pointElement.attribute("y", "0").toDouble();
        double omega = pointElement.attribute("omega", "0").toDouble();
        fprintf(stderr, "Point: (%d,%d) Omega: %d\n", x, y, omega);
        omega = (omega/360) * 2*M_PI;
        surface->addPointAndAngel(x,y,omega);
    }

    return surface;
}

/**
 * @brief MainWindow::loadLineDefinitionFromXML
 * @param lineElement - QDomElement pointing to node (element) containing definition of line component
 * Method having parent node <line> passed as function argument
 * reads into definition of line contained in its child nodes (elements)
 */
Line *MainWindow::loadLineDefinitionFromXML(const QDomElement &lineElement)
{
    //constructing Line object
    Line *line;
    Material *material;

    QDomElement nameElement = lineElement.firstChildElement("name");
    fprintf(stderr, "Loaded line name is: %s\n", nameElement.text().toStdString().c_str());

    QDomNodeList materialList = lineElement.elementsByTagName("material");
    if(materialList.count() != 1)
    {
        QMessageBox::critical(this, "XML file doesn't contain material definition",
                              "XML file doesn't contain material definition!",
                              QMessageBox::Ok);
        return NULL;
    }

    //there should be only one <material> element in xml file
    QDomElement materialElement = materialList.at(0).toElement();
    fprintf(stderr, "Material: %s\n", materialElement.text().toStdString().c_str());

    //depending on material type: 'Reinforcement' or other? we need to read grade or sth else
    if(materialElement.text() == "Reinforcement")
    {
        //we need to read reinforcement steel grade
        QString grade_label = materialElement.attribute("grade", "RB500W"); //default grade of reinforcement steel we assume as RB500W
        GradeOfSteel grade = steelGradeFromString(grade_label);
        fprintf(stderr, "Reinforcement steel grade is: %s.\n", grade_label.toStdString().c_str());

        material = new Reinforcement(grade);
    } else if(materialElement.text() == "FRP")
    {
        material = new FRP();
    }

    line = new Line(nameElement.text().toStdString(), material);

    QDomNodeList pointsList = lineElement.elementsByTagName("points").at(0).toElement().elementsByTagName("point");

    for(int i=0; i<pointsList.count(); i++)
    {
        QDomElement pointElement = pointsList.at(i).toElement();
        double x = pointElement.attribute("x","0").toDouble();
        double y = pointElement.attribute("y","0").toDouble();
        double area = pointElement.attribute("area", "0").toDouble();
        fprintf(stderr, "Point: (%d, %d) Area: %d\n", x, y, area);

        line->addPointAndArea(x,y,area);
    }

    return line;
}

/**
 * @brief MainWindow::loadFiberGroupDefinitionFromXML
 * @param fibergroupElement - QDomElement pointing to node (element) containing definition of fibergroup component
 * Method having parent node <fibergroup> passed as function argument
 * reads into definition of fibergroup contained in its childe nodes (elements)
 */
FiberGroup *MainWindow::loadFiberGroupDefinitionFromXML(const QDomElement &fibergroupElement)
{
    //constructing FiberGroup
    FiberGroup *fibergroup;
    Material *material;

    QDomElement nameElement = fibergroupElement.firstChildElement("name");
    fprintf(stderr, "Loaded fiber group name is: %s\n", nameElement.text().toStdString().c_str());

    QDomNodeList materialList = fibergroupElement.elementsByTagName("material");
    if(materialList.count() != 1)
    {
        QMessageBox::critical(this, "XML file doesn't contain material definition",
                              "XML file doesn't contain definition!",
                              QMessageBox::Ok);
        return NULL;
    }

    //there should be only one <material> element in xml file
    QDomElement materialElement = materialList.at(0).toElement();
    fprintf(stderr, "Material: %s\n", materialElement.text().toStdString().c_str());

    //depending on material type: 'Reinforcement Bars' or other? we need to read grade or sth else
    if(materialElement.text() == "Reinforcement Bars")
    {
        //we need to read reinforcement bars steel grade
        QString grade_label = materialElement.attribute("grade", "RB500W"); //default grade of reinforcement bars steel we assume as RB500W
        GradeOfSteel grade = steelGradeFromString(grade_label);
        fprintf(stderr, "Reinforcement bars steel grade is: %s.\n", grade_label.toStdString().c_str());

        material = new Reinforcement(grade);
    } else if(materialElement.text() == "Fafitis Reinforcement")
    {
        material = new FafitisReinforcement();
    }

    fibergroup = new FiberGroup(nameElement.text().toStdString(), material);

    QDomNodeList pointsList = fibergroupElement.elementsByTagName("points").at(0).toElement().elementsByTagName("point");

    for(int i=0; i<pointsList.count(); i++)
    {
        QDomElement pointElement = pointsList.at(i).toElement();
        double x = pointElement.attribute("x", "0").toDouble();
        double y = pointElement.attribute("y", "0").toDouble();
        double area = pointElement.attribute("area", "0").toDouble();
        fprintf(stderr, "Point: (%d, %d) Area: %d\n", x, y, area);

        fibergroup->addPointAndArea(x,y,area);
    }

    return fibergroup;
}

void MainWindow::on_pushButton_clicked()
{
    QTableWidgetItem *tableItem = ui->sectionDefinitionTableWidget->selectedItems().first();
    int rowIdx = tableItem->row();

    fprintf(stderr, "Selected row index: %d\n", rowIdx);

    if(rowIdx < this->sectionSolver->getNumberOfSurfaces())
    {
        //we get material from surfaces material definitions
        Material *material = this->sectionSolver->getSurfaces()[rowIdx]->getMaterial();
        fprintf(stderr, "Material name: %s\n", material->getMaterialName().c_str());
        MaterialWindow *materialWindow = new MaterialWindow(this);
        materialWindow->setMaterial(material);
        materialWindow->show();
        materialWindow->graphicsViewFitScene(Qt::KeepAspectRatio);

        return;
    }

    rowIdx -= this->sectionSolver->getNumberOfSurfaces();

    if(rowIdx < this->sectionSolver->getNumberOfLines())
    {
        //we get material from lines material definitions
        Material *material = this->sectionSolver->getLines()[rowIdx]->getMaterial();
        fprintf(stderr, "Material name: %s\n", material->getMaterialName().c_str());
        MaterialWindow *materialWindow = new MaterialWindow(this);
        materialWindow->setMaterial(material);
        materialWindow->show();
        materialWindow->graphicsViewFitScene(Qt::KeepAspectRatio);

        return;
    }

    rowIdx -= this->sectionSolver->getNumberOfLines();

    if(rowIdx < this->sectionSolver->getNumberOfFiberGroups())
    {
        //we get material from fibergroups material definitions
        Material *material = this->sectionSolver->getFiberGroups()[rowIdx]->getMaterial();
        fprintf(stderr, "Material name: %s\n", material->getMaterialName().c_str());
        MaterialWindow *materialWindow = new MaterialWindow(this);
        materialWindow->setMaterial(material);
        materialWindow->show();
        materialWindow->graphicsViewFitScene(Qt::KeepAspectRatio);

        return;
    }

    fprintf(stderr, "Row Index in SectionTableWidget is Invalid.\n");
}

void MainWindow::on_pushButton_2_clicked()
{
    solutionWindow = new SolutionWindow(this);
    solutionWindow->show();
}
