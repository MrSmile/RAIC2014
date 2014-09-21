#include <cmath>
#include <vector>
#include <cstdint>
#include <pnglite.h>

using namespace std;

constexpr double pi = 3.14159265358979323846264338327950288;


inline constexpr double sqr(double x)
{
    return x * x;
}

inline constexpr double rem(double x, double y)
{
    return y * (x / y - floor(x / y));
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



constexpr double rinkLeft = 65, rinkRight = 1135, rinkTop = 150, rinkBottom = 770;
constexpr Vec2D rinkCenter = Vec2D(rinkLeft + rinkRight, rinkTop + rinkBottom) / 2;
constexpr double goalCenter = 460, goalHalf = 100, hockeyistRad = 30, puckRad = 20, goalieSpd = 6, goalieRange = 70;
constexpr double hockeyistFrict = 0.02, puckFrict = 0.001, puckBeta = -log(1 - puckFrict);
constexpr double wallCoeff = 0.25, goalieCoeff = 0.325, goalieFrict = 0.1;
constexpr double minDepth = 0.01, depthFactor = 0.8;


struct Sector
{
    double centerAngle, halfSpan;

    Sector(double angle, double span) : centerAngle(angle), halfSpan(span / 2)
    {
    }
};

bool iterateGoalEstimation(double dist, double offs, double spd, Vec2D &dir, bool first)
{
    const double rad = hockeyistRad + puckRad;
    double invSpd = 1 / spd, end = (dist - rad * dir.y) * invSpd, cmp = 0;
    if(first)
    {
        cmp = dist * dir.y - rad;  if(!(goalieSpd * end < cmp))return false;
    }
    double mul = 1 / dir.x, beg = offs * invSpd / dir.y;  end *= mul;
    if(first)
    {
        cmp *= mul;  if(!(goalieSpd * (end - beg) + offs < cmp))return false;
    }
    beg = -log(1 - puckBeta * beg) / puckBeta;
    end = -log(1 - puckBeta * end) / puckBeta;
    double len = goalieSpd * (end - beg) + offs;  if(first && !(len < cmp))return false;
    double w = dist * dist + len * len, z = sqrt(w - rad * rad);
    dir = Vec2D(dist * z - len * rad, len * z + dist * rad) / w;
    return true;
}

Sector estimateGoalAngle(const Vec2D &pos, const Vec2D &spd, double power, bool right)
{
    Vec2D dir = pos, bonus = -spd;
    if(right)
    {
        dir.x = 2 * rinkCenter.x - dir.x;  bonus.x = -bonus.x;
    }
    if(pos.y < goalCenter)
    {
        dir.y = 2 * goalCenter - dir.y;  bonus.y = -bonus.y;
    }
    dir -= Vec2D(rinkLeft, goalCenter - goalHalf);

    const double rad = hockeyistRad + puckRad;
    double dist = dir.x - hockeyistRad;  if(!(dist > rad))return Sector(0, 0);
    double offs = max(0.0, dir.y - 2 * goalHalf);  dir = normalize(dir);

    Vec2D start = dir;
    if(!iterateGoalEstimation(dist, offs, power + bonus * dir, dir, true))return Sector(0, 0);
    for(int i = 0; i < 3; i++)iterateGoalEstimation(dist, offs, power + bonus * dir, dir, false);

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
        Sector score = estimateGoalAngle(corner + Vec2D(j + 0.5, i + 0.5), {0, 0}, spd, true);
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
    createMap("map15.png", 15);  createMap("map20.png", 20);  createMap("map25.png", 25);
    return 0;
}
