#include "glwidget.h"
#include "math.h"

GLWidget::GLWidget(QWidget *parent) :
    QGLWidget(parent)
{
    setMouseTracking(true);
}

void GLWidget::initializeGL()
{
    //sets the clear color to black
    glClearColor(0.0, 0.0, 102.0/255.0, 0.0);
    //Set the drawing color to green
    glColor3f( 0.0f, 1.0f, 0.0f );

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    //glOrtho(0.0, 1.0, 0.0, 1.0, -1.0, 1.0);


    /*
    QSize viewport_size = size();
    glViewport(0, 0, viewport_size.width(), viewport_size.height());

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glFrustum(-1, 1, -1, 1, 5, 7); // near and far match your triangle Z distance

    glMatrixMode(GL_MODELVIEW);*/
}

void GLWidget::resizeGL(int width, int height)
{
    glViewport(0, 0, (GLint)width, (GLint)height);

    // reset the coordinate system
       glMatrixMode( GL_PROJECTION );
       glLoadIdentity();

       // Establish the clipping volume by setting up an orthographic projection
       double range = 100.0;
       float m_aspectRatio = double( width ) / double( height );
       if ( width <=height )
           glOrtho( -range, range, -range / m_aspectRatio, range / m_aspectRatio, range, -range );
       else
           glOrtho( -range * m_aspectRatio, range * m_aspectRatio, -range, range, range, -range );

       glMatrixMode( GL_MODELVIEW );
       glLoadIdentity();
}

void GLWidget::paintGL()
{
    glClear(GL_COLOR_BUFFER_BIT);
   // drawStripes();
   // drawStars();
    glFlush();

    // Clear the buffer with the current clearing color
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
      //glLoadIdentity();

    // Save matrix state and do the custom rotation
       float m_theta = 0.0f;
       float m_phi = 0.0f;
       glPushMatrix();
       glRotatef( m_theta, 1.0f, 0.0f, 0.0f );
       glRotatef( m_phi,   0.0f, 1.0f, 0.0f );

       // Draw some Lines in a helix
       const float pi = 3.141592653f;
       const float twoPi = 2.0f * pi;
       const float piBy2 = 0.5f * pi;
       const float degToRad = pi / 180.0f;
       const float radToDeg = 180.0f / pi;
       glLineWidth( m_lineWidth );
       glBegin( GL_LINE_STRIP );
       float z = -50.0f;
       float angle = 0.0f;
       for ( angle = 0.0f; angle <= twoPi * 3.0f; angle += 0.1f, z += 0.5f )
       {
           float x = 50.0 * sin( angle );
           float y = 50.0 * cos( angle );
           glVertex3f( x, y, z );
       }
       glEnd();

       // Restore the matrix state
       glPopMatrix();

      //glTranslatef(-1.5f,0.0f,-6.0f);

     /* glBegin(GL_TRIANGLES);
        glVertex3f( 0.0f, 1.0f, 0.0f);
        glVertex3f(-1.0f,-1.0f, 0.0f);
        glVertex3f( 1.0f,-1.0f, 0.0f);
      glEnd();
      */

      //glTranslatef(3.0f,0.0f,0.0f);

      glBegin(GL_TRIANGLES);
        glVertex3f(-0.5f, 0.5f, 0.0f);
        glVertex3f( 0.5f,-0.5f, 0.0f);
        glVertex3f(-0.5f,-0.5f, 0.0f);

        glVertex3f( 0.0f, 1.0f, 0.0f);
        glVertex3f(-1.0f,-1.0f, 0.0f);
        glVertex3f( 1.0f,-1.0f, 0.0f);
      glEnd();

   glFlush();
}

void GLWidget::keyPressEvent(QKeyEvent *e)
{

}


void GLWidget::drawStar(float fX, float fY)
{
const float kfPi = 3.1415926535897932384626433832795;
// draw ten triangles
const float kfRadius = 0.0616/2.0;
const float kfInnerRadius = kfRadius*(1.0/(sin((2.0*kfPi)/5.0)*2.0*cos(kfPi/10.0) + sin((3.0*kfPi)/10.0)));
glColor3f(1.0, 1.0, 1.0);

glBegin(GL_TRIANGLE_FAN);
glVertex3f(fX, fY, 0.0);
for (int iVertIndex = 0; iVertIndex < 10; ++iVertIndex)
{
float fAngleStart    = kfPi/2.0 + (iVertIndex*2.0*kfPi)/10.0;
float fAngleEnd        = fAngleStart + kfPi/5.0;
if (iVertIndex % 2 == 0)
{
glVertex3f(fX + kfRadius*cos(fAngleStart)/1.9, fY + kfRadius*sin(fAngleStart), 0.0);
glVertex3f(fX + kfInnerRadius*cos(fAngleEnd)/1.9, fY + kfInnerRadius*sin(fAngleEnd), 0.0);
} else
{
glVertex3f(fX + kfInnerRadius*cos(fAngleStart)/1.9, fY + kfInnerRadius*sin(fAngleStart), 0.0);
glVertex3f(fX + kfRadius*cos(fAngleEnd)/1.9, fY + kfRadius*sin(fAngleEnd), 0.0);
}
}
glEnd();
}

void GLWidget::drawStars()
{
for (int iStarRow = 0; iStarRow < 9; ++iStarRow)
{
float fY = 6.0/13.0 + (iStarRow + 1)*((7.0/13.0)/10);
// alternate between rows of five or six stars
if (iStarRow % 2 == 0)
{
for (int iStarCol = 0; iStarCol < 6; ++iStarCol)
{
drawStar(iStarCol*((0.76/1.9)/6.0) + (0.76/1.9)/12.0, fY);
}
} else
{
for (int iStarCol = 0; iStarCol < 5; ++iStarCol)
{
drawStar((iStarCol + 1)*((0.76/1.9)/6.0), fY);
}
}
}
}

void GLWidget::drawStripes()
{
for (int iStripeIndex = 0; iStripeIndex < 13; ++iStripeIndex)
{
// Alternate stripe colors
if (iStripeIndex % 2 == 0)
{
glColor3f(204.0/255.0, 0.0, 0.0);
} else
{
glColor3f(1.0, 1.0, 1.0);
}

float fStartX    = 0.0;
float fEndX    = 1.0;
float fStartY    = iStripeIndex*(1.0/13.0);
float fEndY    = (iStripeIndex + 1)*(1.0/13.0);

// the last seven stripes are shorter
if (iStripeIndex > 5)
{
fStartX = 0.76/1.9;
}

glBegin(GL_QUADS);
glVertex3f(fStartX, fStartY, 0.0);
glVertex3f(fEndX, fStartY, 0.0);
glVertex3f(fEndX, fEndY, 0.0);
glVertex3f(fStartX, fEndY, 0.0);
glEnd();
}
}
