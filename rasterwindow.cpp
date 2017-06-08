#include "rasterwindow.h"

#include <vector>
#include <cmath>
#include <string>
#include <cassert>

RasterWindow::RasterWindow(QWindow *parent)
    : QWindow(parent)
    , m_update_pending(false)
{
    m_backingStore = new QBackingStore(this);
    create();

    setGeometry(600, 200, 900, 700);
}

bool RasterWindow::event(QEvent *event)
{
    if (event->type() == QEvent::UpdateRequest)
    {
        m_update_pending = false;
        renderNow();
        return true;
    }
    return QWindow::event(event);
}

void RasterWindow::renderLater()
{
    if (!m_update_pending)
    {
        m_update_pending = true;
        QCoreApplication::postEvent(this, new QEvent(QEvent::UpdateRequest));
    }
}

void RasterWindow::resizeEvent(QResizeEvent *resizeEvent)
{
    m_backingStore->resize(resizeEvent->size());
    if (isExposed())
        renderNow();
}

void RasterWindow::exposeEvent(QExposeEvent *)
{
    if (isExposed())
    {
        renderNow();
    }
}


QImage imgScreen(800, 600, QImage::Format_RGBX8888);


struct Point2D
{
    float x;
    float y;
    Point2D()
    {
        x = 0;
        y = 0;
    }

    Point2D(float _x, float _y)
    {
        x = _x;
        y = _y;
    }
};

struct Point3D
{
    float x;
    float y;
    float z;
    float w;
    Point3D()
    {
        x = y = z = 0.0f;
        w = 1.0f;
    }
    Point3D(float _x, float _y, float _z)
    {
        x = _x;
        y = _y;
        z = _z;
        w = 1.0f;
    }

    Point2D getPoint2D() const
    {
        return Point2D(x, y);
    }
};


void putPixel(float _xc, float _yc, int pixel_size = 5, QColor color = QColor(200, 210, 0))
{
    int xc = (int)_xc;
    int yc = (int)_yc;

    for (int dy = -pixel_size/2; dy <= pixel_size/2; dy++)
    {
        for (int dx = -pixel_size/2; dx <= pixel_size/2; dx++)
        {
            int x = xc + dx;
            int y = yc + dy;

            if ((x >= 800) || (x < 0) || (y < 0) || (y >= 600))
                continue;

            imgScreen.setPixelColor( QPoint(x, y), color );
        }
    }
}

void drawPoint(const Point2D &p, int point_size = 5)
{
    putPixel(p.x, p.y, point_size);
}




float degToRad(float deg)
{
    return deg * 3.14159265 / 180.0f;
}

float cotangent(float angle_rad)
{
    return cosf(angle_rad)/sinf(angle_rad);
}

class Matrix
{
    float m[16];  // y, x
public:
    Point3D multiply(const Point3D &p)
    {
        Point3D outPoint;
        outPoint.x = m[0] * p.x + m[1] * p.y + m[2] * p.z + m[3] * p.w;
        outPoint.y = m[4] * p.x + m[5] * p.y + m[6] * p.z + m[7] * p.w;
        outPoint.z = m[8] * p.x + m[9] * p.y + m[10] * p.z + m[11] * p.w;
        outPoint.w = m[12] * p.x + m[13] * p.y + m[14] * p.z + m[15] * p.w;
        return outPoint;
    }

    void zero()
    {
        for (int i = 0; i < 16; i++)
        {
            m[i] = 0;
        }
    }

    void loadIdentity()
    {
        zero();
        m[0]  = 1.0f;
        m[5]  = 1.0f;
        m[10] = 1.0f;
        m[15] = 1.0f;
    }

    void loadTranslation(float tx, float ty, float tz)
    {
        loadIdentity();
        m[3] = tx;
        m[7] = ty;
        m[11] = tz;
    }

    void loadProjection(float fovy_deg = 90.0f, float aspect_wh = 800.0f/600.0f, float zNear = 1.0f, float zFar = 6.0f)
    {
        zero();

        // ...
        float fovy_rad = degToRad(fovy_deg);
        float f = cotangent(fovy_rad / 2.0f);

        m[0] = f / aspect_wh;
        m[5] = f;
        m[10] = (zFar + zNear) / (zNear - zFar);
        m[11] = 2.0f * zFar * zNear / (zNear - zFar);
        m[14] = -1.0f;
    }

    void loadRotationY(float angleY_deg)
    {
        float teta = degToRad(angleY_deg);
        loadIdentity();

        m[0] = cosf(teta);
        m[2] = sinf(teta);
        m[8] = -sinf(teta);
        m[10] = cosf(teta);
    }

    void loadRotationX(float angleX_deg)
    {
        float teta = degToRad(angleX_deg);
        loadIdentity();

        m[5] = cosf(teta);
        m[6] = -sinf(teta);
        m[9] = sinf(teta);
        m[10] = cosf(teta);
    }
};


struct Triangle3D
{
    Point3D v[3];
    // normal
    // texcoord
};


std::vector<Triangle3D> loadModelFromFile(const std::string &filename)
{
    std::vector<Triangle3D> vTriangles;

    float scale = 0.1;
    //for (int i = 0; i < 7; i++)
    for (int i = 0; i < 1; i++)
    {
        //float scalefactor = (float)(i) * scale + 0.3f;
        float scalefactor = 1.0f;
        Triangle3D tri;
        tri.v[0] = Point3D(-2.0f*scalefactor, 1.5f*scalefactor, 0);
        tri.v[1] = Point3D(0.7f*scalefactor,  2.0f*scalefactor, 0);
        tri.v[2] = Point3D(1.7f*scalefactor, -1.8f*scalefactor, 0);

        vTriangles.push_back(tri);
    }

    return vTriangles;
}

std::vector<Point3D> loadPlane(int W, int H, int numPointsX, int numPointsY)
{
    std::vector<Point3D> vPoints;
    vPoints.reserve((numPointsX+1)*(numPointsY+1));

    if (numPointsX < 1)
        numPointsX = 2;
    if (numPointsY < 1)
        numPointsY = 2;

    float delta_x = (float)W / (float)(numPointsX - 1);
    float delta_y = (float)H / (float)(numPointsY - 1);

    for (int iy = 0; iy <= numPointsY; iy++)
    {
        for (int ix = -numPointsX/2; ix <= numPointsX/2; ix++)
        {
            Point3D point;

            point.x = (float)ix * delta_x;
            point.y = (float)iy * delta_y;
            point.z = 0;
            point.w = 1.0f;

            vPoints.push_back(point);
        }
    }

    return vPoints;
}


Point2D transform(const Point3D &point3d, float angleY_deg)
{
    Point2D p;


    Matrix matModel;
    //matModel.loadIdentity();
    matModel.loadRotationY(angleY_deg);



    Point3D point3D_model = matModel.multiply(point3d);

    Matrix matView;
    matView.loadTranslation(0, 0, 5);
    Point3D point3D_modelview = matView.multiply(point3D_model);

    // NDC - clip coordinates
    Matrix matProjection;
    matProjection.loadProjection();
    Point3D point3D_modelviewprojection = matProjection.multiply(point3D_modelview);
    // p2d = P * V * M * p3d

    // Divide by W
    Point3D point3D_ndc = point3D_modelviewprojection;
    point3D_ndc.x = point3D_ndc.x/point3D_ndc.w;
    point3D_ndc.y = point3D_ndc.y/point3D_ndc.w;
    point3D_ndc.z = point3D_ndc.z/point3D_ndc.w;

    float SCREEN_WIDTH = 800.0f;
    float SCREEN_HEIGHT = 600.0f;

    // Window coords
    float Wx = (point3D_ndc.x + 1.0f)*SCREEN_WIDTH/2.0f;
    float Wy = (point3D_ndc.y + 1.0f)*SCREEN_HEIGHT/2.0f;
    float Wz = (point3D_ndc.z + 1.0f)/2.0f;

    p.x = Wx;
    p.y = Wy;

    return p;
}

const float K_EPSILON = 0.0005f;

class Line2D
{
    double k;
    double b;

    Line2D()
    {
        k = 0;
        b = 0;
    }

public:
    Line2D(const Point2D &v1, const Point2D &v2)
    {
        // y = k      * x + b
        k = (v2.y - v1.y) / (v2.x - v1.x);
        if (fabs(v2.x - v1.x) < K_EPSILON)
        {
            k = 0;
        }

        b = -k*v1.x + v1.y;
    }

    static Point2D getIntersectionPoint(const Line2D &l1, const Line2D &l2)
    {
        float x = (l2.b - l1.b) / (l1.k - l2.k);
        if (fabs(l1.k - l2.k) < K_EPSILON)
        {
            x = 0;
        }
        float y = l1.k * x + l1.b;

        return Point2D(x, y);
    }

    Line2D getPerpendicularLine(const Point2D &v)
    {
        Line2D linePerpendicular;

        if (this->k == 0)
        {
            linePerpendicular.k = 1.0f;
            linePerpendicular.b = v.y;
        }
        else
        {
            linePerpendicular.k = -1.0f/this->k;
            linePerpendicular.b = v.x/this->k + v.y;
        }
        return linePerpendicular;
    }
};


class Vector2D
{
    float x, y;

public:
    Vector2D(const Point2D &pStart, const Point2D &pEnd)
    {
        x = pEnd.x - pStart.x;
        y = pEnd.y - pStart.y;
    }

    void normalize()
    {
        float len = length();
        if (len < K_EPSILON)
        {
            //assert(false);
            x = 0;
            y = 0;
        }

        x = x / len;
        y = y / len;
    }

    float length() const
    {
        float len = sqrt(x * x + y * y);
        return len;
    }

    static float DotProduct(const Vector2D &v1, const Vector2D &v2)
    {
        float lenDotProduct = v1.x*v2.x + v1.y*v2.y;
        return lenDotProduct;
    }

    bool isCoDirected(const Vector2D &v2) const
    {
        bool isPointInside = false;

        if (this->length() <= K_EPSILON)
        {
            isPointInside = true;
        }
        else
        {
            Vector2D copyVecNorm = *this;
            copyVecNorm.normalize();

            float projectionLen = Vector2D::DotProduct(v2, copyVecNorm);
            if (fabs(projectionLen) < K_EPSILON)
            {
                isPointInside = true;
            }
            else
            {
                if (projectionLen > 0.0f)
                {
                    isPointInside = true;
                }
                else
                {
                    isPointInside = false;
                }
            }
        }
        return isPointInside;
    }

};


void renderLine(int _x1, int _y1, int _x2, int _y2)
{
    putPixel(_x1, _y1, 1);
    putPixel(_x2, _y2, 1);

    float x1 = (float)_x1;
    float y1 = (float)_y1;

    float x2 = (float)_x2;
    float y2 = (float)_y2;

    // y = k * x + b
    float k = (y2 - y1) / (x2 - x1);
    if (fabs(x2 - x1) < K_EPSILON)
    {
        k = 0;
    }
    float b = -k*x1 + y1;

    int left_x = std::min(_x1, _x2);
    int right_x = std::max(_x1, _x2);

    for (int my_x = left_x; my_x <= right_x; my_x++)
    {
        float x = (float)my_x;
        float y = k * x + b;
        putPixel(x, y, 1);
    }


    int left_y = std::min(_y1, _y2);
    int right_y = std::max(_y1, _y2);

    for (int my_y = left_y; my_y <= right_y; my_y++)
    {
        float y = (float)my_y;
        // float y = k * x + b;
        float x = (y - b) / k;
        if (k == 0)
        {
            x = 0;
        }
        putPixel(x, y, 1);
    }
}

void renderLine(const Point2D &v1, const Point2D &v2)
{
    int x1 = v1.x;
    int y1 = v1.y;
    int x2 = v2.x;
    int y2 = v2.y;
    renderLine(x1, y1, x2, y2);
}


Point3D calculateMedian(const Point3D &v1, const Point3D &v2)
{
    Point3D median = Point3D((v1.x+v2.x)/2.0f, (v1.y+v2.y)/2.0f, (v1.z+v2.z)/2.0f);
    return median;
}


void rasterize(const Triangle3D &tri)
{
#if 0
    // Draw triangle wireframe
    renderLine(tri.v[0].x, tri.v[0].y, tri.v[1].x, tri.v[1].y);
    renderLine(tri.v[1].x, tri.v[1].y, tri.v[2].x, tri.v[2].y);
    renderLine(tri.v[2].x, tri.v[2].y, tri.v[0].x, tri.v[0].y);
#endif

    Point3D m0 = calculateMedian(tri.v[1], tri.v[2]);
    Point3D m1 = calculateMedian(tri.v[0], tri.v[2]);
    Point3D m2 = calculateMedian(tri.v[0], tri.v[1]);

#if 0
    drawPoint(m0.getPoint2D(), 5);
    drawPoint(m1.getPoint2D(), 5);
    drawPoint(m2.getPoint2D(), 5);

    renderLine(tri.v[0].getPoint2D(), m0.getPoint2D());
    renderLine(tri.v[1].getPoint2D(), m1.getPoint2D());
    renderLine(tri.v[2].getPoint2D(), m2.getPoint2D());
#endif

    Line2D l0(tri.v[0].getPoint2D(), m0.getPoint2D());
    Line2D l1(tri.v[1].getPoint2D(), m1.getPoint2D());
    Point2D barycenentric = Line2D::getIntersectionPoint(l0, l1);

#if 0
    drawPoint(barycenentric, 7);
#endif

    // Calculate normal touch points
    Line2D L01(tri.v[0].getPoint2D(), tri.v[1].getPoint2D());
    Line2D perp_L01 = L01.getPerpendicularLine(barycenentric);
    Point2D n2 = Line2D::getIntersectionPoint(L01, perp_L01);

    Line2D L12(tri.v[1].getPoint2D(), tri.v[2].getPoint2D());
    Line2D perp_L12 = L12.getPerpendicularLine(barycenentric);
    Point2D n0 = Line2D::getIntersectionPoint(L12, perp_L12);

    Line2D L02(tri.v[0].getPoint2D(), tri.v[2].getPoint2D());
    Line2D perp_L02 = L02.getPerpendicularLine(barycenentric);
    Point2D n1 = Line2D::getIntersectionPoint(L02, perp_L02);

#if 0
    drawPoint(n0, 3);
    drawPoint(n1, 3);
    drawPoint(n2, 3);

    renderLine(barycenentric, n0);
    renderLine(barycenentric, n1);
    renderLine(barycenentric, n2);
#endif

    // Side 1-2

    bool isPointInside = false;

    Vector2D vecN0(n0, barycenentric);
    vecN0.normalize();

    Vector2D vecN1(n1, barycenentric);
    vecN1.normalize();

    Vector2D vecN2(n2, barycenentric);
    vecN2.normalize();

    float bbox_left = tri.v[0].x;
    float bbox_right = tri.v[0].x;
    float bbox_up = tri.v[0].y;
    float bbox_bottom = tri.v[0].y;

    for (int i = 1; i < 3; i++)
    {
        bbox_left = std::min(bbox_left, tri.v[i].x);
        bbox_right = std::max(bbox_right, tri.v[i].x);

        bbox_up = std::min(bbox_up, tri.v[i].y);
        bbox_bottom = std::max(bbox_bottom, tri.v[i].y);
    }

    int x_start = bbox_left - 1;
    int x_end = bbox_right + 1;

    x_start = std::max(x_start, 0);
    x_start = std::min(x_start, 800 - 1);
    x_end = std::max(x_end, 0);
    x_end = std::min(x_end, 800 - 1);

    int y_start = bbox_up - 1;
    int y_end = bbox_bottom + 1;

    y_start = std::max(y_start, 0);
    y_start = std::min(y_start, 600 - 1);
    y_end = std::max(y_end, 0);
    y_end = std::min(y_end, 600 - 1);

    for (int x = x_start; x <= x_end; x++)
    {
        //int y = 200;
        for (int y = y_start; y <= y_end; y++)
        {
            Point2D interestPoint(x, y);

            Vector2D vecInter0(n0, interestPoint);
            Vector2D vecInter1(n1, interestPoint);
            Vector2D vecInter2(n2, interestPoint);

            bool isPointInside0 = vecInter0.isCoDirected(vecN0);
            bool isPointInside1 = vecInter1.isCoDirected(vecN1);
            bool isPointInside2 = vecInter2.isCoDirected(vecN2);

            bool isRasterizeable = isPointInside0 && isPointInside1 && isPointInside2;
            if (isRasterizeable)
            {
                if (y % 2)
                {
                    putPixel(x, y, 1);
                }
                else
                {
                    putPixel(x, y, 1, QColor(30, 150, 70));
                }
                putPixel(x, y, 1);
            }
#if 0
            if (!isRasterizeable)
            {
                putPixel(x, y, 1, QColor(30, 150, 70));
            }
#endif
        }

    }


}




void render3D()
{
    static int frameNumber = 0;
    static float angleY = -10.0f;

    imgScreen.fill(QColor(0, 0, 96));

    //std::vector<Point3D> vPoints3D = loadPlane(1.0f, 2.0f, 7, 10);

    std::vector<Triangle3D> vTriangles = loadModelFromFile("nofile.3ds.lalalalal");

    std::vector<Point3D> vPoints3D;

    for (const Triangle3D &tri: vTriangles)
    {
        for (int idxVert = 0; idxVert < 3; idxVert++)
        {
            vPoints3D.push_back(tri.v[idxVert]);
        }
    }

    angleY += 10.0f;

    std::vector<Point2D> points2D;
    points2D.resize(vPoints3D.size());

    for (int i = 0; i < vPoints3D.size(); i++)
    {
        points2D[i] = transform(vPoints3D[i], angleY);
    }

    // Triangle Assembly
    std::vector<Triangle3D> vTrianglesTransformed;
    int num_triangles = vTriangles.size();
    for (int i = 0; i < num_triangles; i++)
    {
        Triangle3D tri;

        float x = 0;
        float y = 0;
        //....
        x = points2D[i*3+0].x;
        y = points2D[i*3+0].y;
        tri.v[0] = Point3D(x, y, 0);

        x = points2D[i*3+1].x;
        y = points2D[i*3+1].y;
        tri.v[1] = Point3D(x, y, 0);

        x = points2D[i*3+2].x;
        y = points2D[i*3+2].y;
        tri.v[2] = Point3D(x, y, 0);

        vTrianglesTransformed.push_back(tri);
    }

    // RasterizeTriangles
    for (const Triangle3D &tri: vTrianglesTransformed)
    {
        rasterize(tri);
    }



    /*
    for (Point2D & p : points)
    {
        p.x = rand() % 700;
        p.y = rand() % 500;
    }
    */

#if 0
    for (Point2D & p : points2D)
    {
        drawPoint(p);
    }
#endif

    /*
    for (int i = 0; i < 800; i++)
    {
        for (int y = 0; y < 500; y++)
        {
            imgScreen.setPixelColor(QPoint(i, y), QColor((10+i + j)%255, abs((256+i-y))%255, (128+i*2+y)%255));
        }
    }
    */
    frameNumber++;
}



void RasterWindow::renderNow()
{
    if (!isExposed())
        return;

    QRect rect(0, 0, width(), height());
    m_backingStore->beginPaint(rect);

    QPaintDevice *device = m_backingStore->paintDevice();
    QPainter painter(device);

    painter.fillRect(0, 0, width(), height(), Qt::red);

    render3D();

    //renderLine(200, 100, 700, 400);

    painter.drawImage(0, 0, imgScreen);
    render(&painter);

    m_backingStore->endPaint();
    m_backingStore->flush(rect);
}

void RasterWindow::render(QPainter *painter)
{
    //painter->drawText(QRectF(0, 0, width(), height()), Qt::AlignCenter, QStringLiteral("QWindow"));
}

