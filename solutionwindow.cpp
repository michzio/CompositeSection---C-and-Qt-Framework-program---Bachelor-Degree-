#include "solutionwindow.h"
#include "ui_solutionwindow.h"
#include <QTableWidget>
#include <QStackedWidget>
#include <QGraphicsView>
#include <QGraphicsPolygonItem>
#include "mainwindow.h"
#include <QtAlgorithms>
#include "flags.h"
#include <QFileDialog>
#include <QTextStream>


SolutionWindow::SolutionWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::SolutionWindow)
{
    ui->setupUi(this);
    parentWindow = parent;

    //setting GraphicsScene for drawing interaction curves
    //Mx-My, N-Mx, N-My
    sceneMxMy = new GridGraphicsScene(this, QString("Mx"), QString("My"));
    ui->graphicsViewMxMy->setScene(sceneMxMy);
    sceneNMx = new GridGraphicsScene(this, QString("Mx"), QString("N"));
    ui->graphicsViewNMx->setScene(sceneNMx);
    sceneNMy = new GridGraphicsScene(this, QString("My"), QString("N"));
    ui->graphicsViewNMy->setScene(sceneNMy);
    sceneNM = new GridGraphicsScene(this, QString("M"), QString("N"));
    ui->graphicsViewNM->setScene(sceneNM);
    /*
    MainWindow * mainwindow = static_cast<MainWindow *>(parent);
    QList<Vertex3D *> vertices = mainwindow->getSectionSolver()->getVertices();
    fprintf(stderr, "\n################ VERTICES UNSORTED ################\n");

    for(QList<Vertex3D *>::const_iterator itr = vertices.cbegin(); itr < vertices.cend(); itr++ )
    {
        fprintf(stderr, "[%d, %d, %d]\n", (*itr)->Mx(), (*itr)->My(), (*itr)->N());
    }

   mainwindow->getSectionSolver()->sortVerticesByN();
   vertices = mainwindow->getSectionSolver()->getVertices();
   fprintf(stderr, "\n################ VERTICES SORTED BY N ################\n");

   for(QList<Vertex3D *>::const_iterator itr = vertices.cbegin(); itr < vertices.cend(); itr++ )
   {
       fprintf(stderr, "[%d, %d, %d] fi: %g\n", (*itr)->Mx(), (*itr)->My(), (*itr)->N(),
                                                (*itr)->argumentMxMy());
   }

   mainwindow->getSectionSolver()->sortVerticesByMx();
   vertices = mainwindow->getSectionSolver()->getVertices();
   fprintf(stderr, "\n################ VERTICES SORTED BY Mx ################\n");

   for(QList<Vertex3D *>::const_iterator itr = vertices.cbegin(); itr < vertices.cend(); itr++ )
   {
       fprintf(stderr, "[%d, %d, %d]\n", (*itr)->Mx(), (*itr)->My(), (*itr)->N());
   }

   mainwindow->getSectionSolver()->sortVerticesByMy();
   vertices = mainwindow->getSectionSolver()->getVertices();
   fprintf(stderr, "\n################ VERTICES SORTED BY My ################\n");

   for(QList<Vertex3D *>::const_iterator itr = vertices.cbegin(); itr < vertices.cend(); itr++ )
   {
       fprintf(stderr, "[%d, %d, %d]\n", (*itr)->Mx(), (*itr)->My(), (*itr)->N());
   }
   */

}

void SolutionWindow::drawInteractionCurveMxMyForN(int N)
{
    QPolygonF polygon;

    MainWindow *mainWindow = static_cast<MainWindow *>(parentWindow);
    mainWindow->getSectionSolver()->sortVerticesByN();
    QList<Vertex3D *> vertices = mainWindow->getSectionSolver()->getVertices();

    QList<Vertex3D *> polygon_vertices;

    //searching vertices for axial force equal to N
    fprintf(stderr, "Constructing interaction curve for N = %d.\n", N);
    for(QList<Vertex3D *>::const_iterator itr = vertices.cbegin(); itr < vertices.cend(); itr++ )
    {
        if(DEBUG)
            fprintf(stderr, "[%d, %d, %d]\n", (*itr)->Mx(), (*itr)->My(), (*itr)->N());

        if((*itr)->N() < N-N_RANGE) continue;
        else if((*itr)->N() > N+N_RANGE) continue;
        else {
           fprintf(stderr, "[%d, %d, %d]\n", (*itr)->Mx(), (*itr)->My(), (*itr)->N());
           polygon_vertices.push_back(new Vertex3D((*itr)->N(), (*itr)->Mx(), (*itr)->My() ));
        }
    }

    qSort(polygon_vertices.begin(), polygon_vertices.end(), Vertex3DComparatorByAlfaInMxMy());

    int maxValue = 0;

    for(QList<Vertex3D *>::const_iterator itr = polygon_vertices.cbegin(); itr < polygon_vertices.cend(); itr++ )
    {
        if((*itr)->My() > maxValue) maxValue = (*itr)->My();
        if((*itr)->Mx() > maxValue) maxValue = (*itr)->Mx();
        //we mirror the plot around x-x axis because Graphics scene assume increasing values
        //in top-to-bottom direction and left-to-right direction
        polygon.append(QPointF((*itr)->Mx(), -(*itr)->My()));
    }

    QPen pen(Qt::black);
    int scale = 2*maxValue/600;
    pen.setWidth(1*scale);
    pen.setCosmetic(false);

    QGraphicsPolygonItem *polygonItem =
              sceneMxMy->addPolygon(polygon, pen, QBrush(Qt::transparent));

    ui->graphicsViewMxMy->fitInView(sceneMxMy->sceneRect(), Qt::KeepAspectRatio);
}

void SolutionWindow::drawInteractionCurveMyNForMx(int Mx)
{
    QPolygonF polygon;

    MainWindow *mainWindow = static_cast<MainWindow *>(parentWindow);
    mainWindow->getSectionSolver()->sortVerticesByMx();
    QList<Vertex3D *> vertices = mainWindow->getSectionSolver()->getVertices();

    QList<Vertex3D *> polygon_vertices;

    //searching vertices for moment Mx
    fprintf(stderr, "Constructing interaction curve for Mx = %d.\n", Mx);
    for(QList<Vertex3D *>::const_iterator itr = vertices.cbegin(); itr < vertices.cend(); itr++ )
    {
        if(DEBUG)
            fprintf(stderr, "[%d, %d, %d]\n", (*itr)->Mx(), (*itr)->My(), (*itr)->N());

        if((*itr)->Mx() < Mx) continue;
        else if((*itr)->Mx() > Mx) continue;
        else {
           polygon_vertices.push_back(new Vertex3D((*itr)->N(), (*itr)->Mx(), (*itr)->My() ));
        }
    }

    qSort(polygon_vertices.begin(), polygon_vertices.end(), Vertex3DComparatorByAlfaInMyN());

    int maxValue = 0;

    for(QList<Vertex3D *>::const_iterator itr = polygon_vertices.cbegin(); itr < polygon_vertices.cend(); itr++ )
    {
        if((*itr)->My()*10 > maxValue) maxValue = (*itr)->My()*10;
        if((*itr)->N() > maxValue) maxValue = (*itr)->N();
        polygon.append(QPointF((*itr)->My()*10, -(*itr)->N()));
    }

    QPen pen(Qt::black);
    int scale = 2*maxValue/600;
    pen.setWidth(1*scale);
    pen.setCosmetic(false);

    QGraphicsPolygonItem *polygonItem =
              sceneNMy->addPolygon(polygon, pen, QBrush(Qt::transparent));

    ui->graphicsViewNMy->fitInView(sceneNMy->sceneRect(), Qt::KeepAspectRatioByExpanding);

}

void SolutionWindow::drawInteractionCurveMxNForMy(int My)
{
    QPolygonF polygon;

    MainWindow *mainWindow = static_cast<MainWindow *>(parentWindow);
    mainWindow->getSectionSolver()->sortVerticesByMy();
    QList<Vertex3D *> vertices = mainWindow->getSectionSolver()->getVertices();

    QList<Vertex3D *> polygon_vertices;

    //searching vertices for moment My
    fprintf(stderr, "Constructing interaction curve for My = %d.\n", My);
    for(QList<Vertex3D *>::const_iterator itr = vertices.cbegin(); itr < vertices.cend(); itr++ )
    {
        fprintf(stderr, "[%d, %d, %d]\n", (*itr)->Mx(), (*itr)->My(), (*itr)->N());

        if((*itr)->My() < My) continue;
        else if((*itr)->My() > My) continue;
        else {
           polygon_vertices.push_back(new Vertex3D((*itr)->N(), (*itr)->Mx(), (*itr)->My() ));
        }
    }

    qSort(polygon_vertices.begin(), polygon_vertices.end(), Vertex3DComparatorByAlfaInMxN());

    int maxValue = 0;

    for(QList<Vertex3D *>::const_iterator itr = polygon_vertices.cbegin(); itr < polygon_vertices.cend(); itr++ )
    {
        if((*itr)->Mx()*10 > maxValue) maxValue = (*itr)->Mx()*10;
        if((*itr)->N() > maxValue) maxValue = (*itr)->N();
        //x-y ratio 10:1
        polygon.append(QPointF((*itr)->Mx()*10, -(*itr)->N()));
    }

    QPen pen(Qt::black);
    int scale = 2*maxValue/600;
    pen.setWidth(1*scale);
    pen.setCosmetic(false);

    QGraphicsPolygonItem *polygonItem =
              sceneNMx->addPolygon(polygon, pen, QBrush(Qt::transparent));

     ui->graphicsViewNMx->fitInView(sceneNMx->sceneRect(), Qt::KeepAspectRatioByExpanding);
     //ui->graphicsViewNMx->fitInView(0, 0, sceneNMx->width(), sceneNMx->height());
}

void SolutionWindow::drawInteractionCurveMNForAlfa(int alfa)
{
    QPolygonF polygon;

    MainWindow *mainWindow = static_cast<MainWindow *>(parentWindow);
    mainWindow->getSectionSolver()->sortVerticesByAlfa();
    QList<Vertex3D *>vertices = mainWindow->getSectionSolver()->getVertices();

    QList<Vertex3D *> polygon_vertices;

    //searching vertices for moment My
    fprintf(stderr, "Constructing interaction curve for Alfa = %d.\n", alfa);
    for(QList<Vertex3D *>::const_iterator itr = vertices.cbegin(); itr < vertices.cend(); itr++)
    {
        int vertexAlfa = (int)((*itr)->argumentMxMy() *180/M_PI)%180;
        fprintf(stderr, "[%d, %d, %d]\n", vertexAlfa, (*itr)->M(), (*itr)->N());
        if(vertexAlfa < alfa) continue;
        else if(vertexAlfa > alfa) continue;
        else {
            polygon_vertices.push_back(new Vertex3D((*itr)->N(), (*itr)->Mx(), (*itr)->My() ));
        }
    }

    int maxValue = 0;

    //aditional sorting if values will be selected for range of alfa values [alfa-x, alfa+x]
    qSort(polygon_vertices.begin(), polygon_vertices.end(), Vertex3DComparatorByBetaInMN() );
    for(QList<Vertex3D *>::const_iterator itr = polygon_vertices.cbegin(); itr < polygon_vertices.cend(); itr++)
    {
        if((*itr)->M()*10 > maxValue) maxValue = (*itr)->M()*10;
        if((*itr)->N() > maxValue) maxValue = (*itr)->N();

        polygon.append( QPointF( (*itr)->M()*10, -(*itr)->N() ));
    }

    QPen pen(Qt::black);
    int scale = 2*maxValue/600;
    pen.setWidth(1*scale);
    pen.setCosmetic(false);

    QGraphicsPolygonItem *polygonItem =
            sceneNM->addPolygon(polygon, pen, QBrush(Qt::transparent));
    ui->graphicsViewNM->fitInView(sceneNM->sceneRect(), Qt::KeepAspectRatio);
}

SolutionWindow::~SolutionWindow()
{
    delete ui;
}

void SolutionWindow::on_pushButtonDrawForMx_clicked()
{
    int Mx =  ui->lineEditMxValue->text().toInt();
    drawInteractionCurveMyNForMx(Mx);
}

void SolutionWindow::on_pushButtonDrawForMy_clicked()
{
    int My =  ui->lineEditMyValue->text().toInt();
    drawInteractionCurveMxNForMy(My);
}

void SolutionWindow::on_pushButtonDrawForN_clicked()
{
    int N =  ui->lineEditNValue->text().toInt();
    drawInteractionCurveMxMyForN(N);
}

void SolutionWindow::on_pushButtonDrawForAlfa_clicked()
{
    int alfa = ui->lineEditAlfaValue->text().toInt();
    drawInteractionCurveMNForAlfa(alfa%180); //[180,360] -> [0,180]
}

void SolutionWindow::on_pushButtonExport_clicked()
{
    //exporting vertices to .txt file

    QString fileName;
    int mode = 0;

    if(ui->radioButtonN->isChecked()) {
         fileName = "N.txt";
         mode = 0;
    } else if(ui->radioButtonMx->isChecked()) {
         fileName = "Mx.txt";
         mode = 1;
    } else if(ui->radioButtonMy->isChecked()) {
         fileName = "My.txt";
         mode = 2;
    } else if(ui->radioButtonAll->isChecked()) {
         fileName = "MxMyN.txt";
         mode = 3;
    } else {
         fprintf(stderr, "Any radiobutton hasn't been selected while exporting file.\n");
    }

   //creating Save As... QFileDialog
   fileName = QFileDialog::getSaveFileName(this, tr("Save File"),
                              fileName,
                              tr("Text files (*.txt)"));
   fprintf(stderr, "File name selected by user is: %s\n", fileName.toStdString().c_str());

   //creating QFile object
   QFile file(fileName);

   //opening file to write text only
   if(file.open(QIODevice::WriteOnly | QIODevice::Text)) {
       //creating output stream to file object
       QTextStream out_stream( &file);

       //getting access to list of vertices
       MainWindow *mainWindow = static_cast<MainWindow *>(parentWindow);
       mainWindow->getSectionSolver()->sortVerticesByAlfa();
       QList<Vertex3D *>vertices = mainWindow->getSectionSolver()->getVertices();

       //outputing data to file depending on selected mode
       switch (mode) {
       case 0: {
           //output only N values
            for(QList<Vertex3D *>::const_iterator itr = vertices.cbegin(); itr < vertices.cend(); itr++)
            {
                    //loop through each element in QList<Vertex3d *> using const_iterator
                    out_stream << (*itr)->N() << " ";
            }
           break;
       }
       case 1: {
           //output only Mx values
            for(QList<Vertex3D *>::const_iterator itr = vertices.cbegin(); itr < vertices.cend(); itr++)
            {
                    //loop through each element in QList<Vertex3d *> using const_iterator
                    out_stream << (*itr)->Mx() << " ";
            }
           break;
       }
       case 2: {
           //output only My values
           for(QList<Vertex3D *>::const_iterator itr = vertices.cbegin(); itr < vertices.cend(); itr++)
           {
                   //loop through each element in QList<Vertex3d *> using const_iterator
                   out_stream << (*itr)->My() << " ";
           }
           break;
       }
       case 3: {
           //output oall MxMyN values
           for(QList<Vertex3D *>::const_iterator itr = vertices.cbegin(); itr < vertices.cend(); itr++)
           {

                   //loop through each element in QList<Vertex3d *> using const_iterator
                   out_stream << (*itr)->Mx() << "," <<  (*itr)->My() << "," << (*itr)->N() << endl;
           }
           break;
       }
       default:
           fprintf(stderr, "Inappropriate writing mode in Export vertices values method.\n");
           break;
       }
       //closing file object
       file.close();
   }

}

void SolutionWindow::on_pushButtonImport_clicked()
{
    //importing verticies from .txt file
    QString fileName = QFileDialog::getOpenFileName(this, tr("Open File"),
                                                    "",  tr("Text files (*.txt)"));

    fprintf(stderr, "File name selected by user is: %s\n", fileName.toStdString().c_str());

    //creating QFile object with handle to file fileName
    QFile file(fileName);

    //opening file to read text only
   if( file.open(QIODevice::ReadOnly | QIODevice::Text) ) {
       //creating input stream from file object
        QTextStream in_stream(&file);

        //getting access to list of vertices
        MainWindow *mainWindow = static_cast<MainWindow *>(parentWindow);
        mainWindow->getSectionSolver()->sortVerticesByAlfa();
        QList<Vertex3D *> &vertices = mainWindow->getSectionSolver()->getVertices();

        //loop through each line in input text
        while(!in_stream.atEnd()) {
                //reading line to QString object
                 QString line = in_stream.readLine();
                 //spliting each line Mx,My,N into list
                 QStringList forces = line.split(",");
                 Vertex3D *vertex = new Vertex3D(forces[2].toDouble(),
                                                 forces[0].toDouble(),
                                                 forces[1].toDouble());
                 fprintf(stderr, "Adding vertex: N = %g, Mx = %g, My = %g\n", forces[2].toDouble(), forces[0].toDouble(), forces[1].toDouble());
                 vertices.push_back(vertex);

        }
   }

}
