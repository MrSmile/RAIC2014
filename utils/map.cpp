#include <cmath>
#include <vector>
#include <cstdint>
#include <pnglite.h>

using namespace std;

constexpr double pi = 3.14159265358979323846264338327950288;


inline double sqr(double x)
{
    return x * x;
}

inline double rem(double x, double y)
{
    x /= y;  return y * (x - floor(x));
}


struct Vec2D
{
    double x, y;

    Vec2D()
    {
    }

    constexpr Vec2D(const Vec2D &v) : x(v.x), y(v.y)
    {
    }

    constexpr Vec2D(double x_, double y_) : x(x_), y(y_)
    {
    }

    Vec2D &operator = (const Vec2D &v)
    {
        x = v.x;  y = v.y;  return *this;
    }

    constexpr Vec2D operator + (const Vec2D &v) const
    {
        return Vec2D(x + v.x, y + v.y);
    }

    Vec2D &operator += (const Vec2D &v)
    {
        x += v.x;  y += v.y;  return *this;
    }

    constexpr Vec2D operator - (const Vec2D &v) const
    {
        return Vec2D(x - v.x, y - v.y);
    }

    Vec2D &operator -= (const Vec2D &v)
    {
        x -= v.x;  y -= v.y;  return *this;
    }

    constexpr Vec2D operator - () const
    {
        return Vec2D(-x, -y);
    }

    constexpr Vec2D operator * (double a) const
    {
        return Vec2D(a * x, a * y);
    }

    Vec2D &operator *= (double a)
    {
        x *= a;  y *= a;  return *this;
    }

    constexpr double operator * (const Vec2D &v) const
    {
        return x * v.x + y * v.y;
    }

    constexpr Vec2D operator / (double a) const
    {
        return (*this) * (1 / a);
    }

    Vec2D operator /= (double a)
    {
        return (*this) *= (1 / a);
    }

    constexpr Vec2D operator ~ () const
    {
        return Vec2D(y, -x);
    }

    constexpr double operator % (const Vec2D &v) const
    {
        return *this * ~v;
    }

    constexpr double sqr() const
    {
        return x * x + y * y;
    }

    constexpr double len() const
    {
        return std::sqrt(x * x + y * y);
    }
};

inline constexpr Vec2D operator * (double a, const Vec2D &v)
{
    return v * a;
}

inline constexpr Vec2D normalize(const Vec2D &v)
{
    return v / v.len();
}

inline constexpr Vec2D sincos(double angle)
{
    return Vec2D(cos(angle), sin(angle));
}

inline constexpr Vec2D conj(const Vec2D &v)
{
    return Vec2D(v.x, -v.y);
}



constexpr double rinkLeft = 65, rinkRight = 959, rinkTop = 170, rinkBottom = 723;
constexpr double goalCenter = 446.5, goalHalf = 100, hockeyistRad = 30, puckRad = 20, goalieSpd = 6, goalieRange = 70;
constexpr double hockeyistFrict = 0.02, puckFrict = 0.001, beta = -log(1 - puckFrict);
constexpr double wallCoeff = 0.25, goalieCoeff = 0.325, goalieFrict = 0.1;
constexpr double minDepth = 0.01, depthFactor = 0.8;


struct Sector
{
    double centerAngle, halfSpan;

    Sector(double angle, double span) : centerAngle(angle), halfSpan(span / 2)
    {
    }
};

bool iterateGoalEstimation(double dist, double offs, double invSpd, Vec2D &dir, bool first)
{
    const double rad = hockeyistRad + puckRad;
    double end = (dist - rad * dir.y) * invSpd, cmp;
    if(first)
    {
        cmp = dist * dir.y - rad;  if(!(goalieSpd * end < cmp))return false;
    }
    double mul = 1 / dir.x, beg = offs * invSpd / dir.y;  end *= mul;
    if(first)
    {
        cmp *= mul;  if(!(goalieSpd * (end - beg) + offs < cmp))return false;
    }
    beg = -log(1 - beta * beg) / beta;  end = -log(1 - beta * end) / beta;
    double len = goalieSpd * (end - beg) + offs;  if(first && !(len < cmp))return false;
    double w = dist * dist + len * len, z = sqrt(w - rad * rad);
    dir = Vec2D(dist * z - len * rad, len * z + dist * rad) / w;
    //cout << atan2(dir.y, dir.x) << ' ';
    return true;
}

Sector estimateGoalAngle(const Vec2D &pos, double spd, bool right)
{
    Vec2D dir = pos;
    if(right)dir.x = rinkLeft + rinkRight - dir.x;
    if(pos.y < goalCenter)dir.y = 2 * goalCenter - dir.y;
    dir -= Vec2D(rinkLeft, goalCenter - goalHalf);

    const double rad = hockeyistRad + puckRad;
    double dist = dir.x - hockeyistRad;  if(!(dist > rad))return Sector(0, 0);
    double offs = max(0.0, dir.y - 2 * goalHalf);  dir = normalize(dir);

    Vec2D start = dir;  double invSpd = 1 / spd;
    if(!iterateGoalEstimation(dist, offs, invSpd, dir, true))return Sector(0, 0);
    for(int i = 0; i < 3; i++)iterateGoalEstimation(dist, offs, invSpd, dir, false);
    //cout << atan2(start.y, start.x) << endl;

    double span = atan2(dir % start, dir * start);  dir += start;
    if(!right)dir.x = -dir.x;  if(!(pos.y < goalCenter))dir.y = -dir.y;
    return Sector(atan2(dir.y, dir.x), span);
}


bool createMap(const char *file, double spd)
{
    constexpr Vec2D corner(rinkLeft + puckRad, rinkTop + puckRad);
    constexpr Vec2D size = Vec2D(rinkRight - puckRad, rinkBottom - puckRad) - corner;
    constexpr int width = lround(size.x), height = lround(size.y);

    vector<uint8_t> img(width * height);
    for(int i = 0, k = 0; i < height; i++)for(int j = 0; j < width; j++, k++)
    {
        Sector score = estimateGoalAngle(corner + Vec2D(j + 0.5, i + 0.5), spd, true);
        img[k] = min(255l, lround(255 * (30 / pi) * score.halfSpan));
    }

    png_t png;
    if(png_open_file_write(&png, file))return false;
    bool res = !png_set_data(&png, width, height, 8, PNG_GREYSCALE, img.data());
    return !png_close_file(&png) && res;
}

int main()
{
    if(png_init(0, 0))return -1;
    createMap("map15.png", 15);  createMap("map20.png", 20);
    createMap("map25.png", 25);  createMap("map30.png", 30);
    return 0;
}
