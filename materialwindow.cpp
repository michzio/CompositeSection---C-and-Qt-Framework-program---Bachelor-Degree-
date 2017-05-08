#include "materialwindow.h"
#include "ui_materialwindow.h"

MaterialWindow::MaterialWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MaterialWindow)
{
    ui->setupUi(this);
    parentWindow = parent;

    //setting GraphicsScene for plotting material definition's sigma-epsilon function
    scene = new GridGraphicsScene(this, QString("epsilon"), QString("sigma"));
    ui->graphicsViewMaterial->setScene(scene);
}

MaterialWindow::MaterialWindow(Material *m, QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MaterialWindow)
{
    ui->setupUi(this);
    parentWindow = parent;

    //setting GraphicsScene for plotting material definition's sigma-epsilon function
    scene = new GridGraphicsScene(this, QString("epsilon"), QString("sigma"));
    ui->graphicsViewMaterial->setScene(scene);

    setMaterial(m);
}

MaterialWindow::~MaterialWindow()
{
    delete ui;
}

void MaterialWindow::setMaterial(Material *m)
{
    material = m;
    //plot sigma-epsilon function
    plotSigmaEpsilon();
}


void MaterialWindow::plotSigmaEpsilon(void)
{
   double eps_cu = material->compressiveStrainLimit_eps_cu();
   double eps_tu = material->tensileStrainLimit_eps_tu();
   double epsilonScale = 20000;
   double sigmaScale = 0.5;

   if(material->getMaterialType() == CONCRETE) {
       epsilonScale = 200000;
       sigmaScale = 20;
   }

   fprintf(stderr, "Material Definition Strain Limits eps_cu: %g, eps_tu: %g\n", eps_cu, eps_tu);

   //setting GraphicsScene GRID_STEP_X and GRID_STEP_Y values;
   /*
   int rangeY = std::abs(material->sigmaEpsilon(eps_tu) - material->sigmaEpsilon(eps_cu))*sigmaScale;
   int stepY = rangeY/10 - (rangeY/10)%10;
   scene->setGridStepY(stepY);
   int rangeX = std::abs(eps_tu - eps_cu)*epsilonScale;
   int stepX = rangeX/10 - (rangeX/10)%10;
   scene->setGridStepX(stepX);
   fprintf(stderr, "StepX: %d, StepY: %d\n", stepX, stepY);
   */
   scene->setGridStepY(50);
   scene->setGridStepX(100);
   scene->setGridScaleX(epsilonScale);
   scene->setGridScaleY(sigmaScale);

   for(double eps_0 = eps_cu, eps_1 = eps_cu + 0.0001;
       eps_1 < eps_tu; eps_0 = eps_1, eps_1 += 0.0001)
   {
       scene->addLine(eps_0*epsilonScale, -material->sigmaEpsilon(eps_0)*sigmaScale, eps_1*epsilonScale, -material->sigmaEpsilon(eps_1)*sigmaScale, QPen(Qt::black));
   }
}

void MaterialWindow::graphicsViewFitScene(Qt::AspectRatioMode aspectRatioMode)
{
    ui->graphicsViewMaterial->fitInView(scene->sceneRect().x() + 50, scene->sceneRect().y() + 50,
                                        scene->sceneRect().width() + 50, scene->sceneRect().height() + 50,
                                        aspectRatioMode);
}
