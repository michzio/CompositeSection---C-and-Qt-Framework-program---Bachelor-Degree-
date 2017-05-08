#ifndef GLWIDGET_H
#define GLWIDGET_H

#include <QGLWidget>

class GLWidget : public QGLWidget
{
    Q_OBJECT
    float m_lineWidth = 1.0f;
public:
    explicit GLWidget(QWidget *parent = 0);

protected:
void initializeGL();
void resizeGL(int width, int height);
void paintGL();
void keyPressEvent(QKeyEvent * e);

private:
void drawStar(float fX, float fY);
void drawStars();
void drawStripes();

};

#endif // GLWIDGET_H
