#include <cmath>
#include <vector>
#include <iostream>
#include <iomanip>

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


struct Info
{
    static constexpr int n = 1;
    static constexpr double accelF = 25.0 / 216, accelB = -0.75 * accelF;
    static constexpr double frict = 0.02, turn = pi / 60;

    Vec2D pos, spd;
    double angle;

    void nextStep(int acc, int rot)
    {
        for(int i = 0; i < n; i++)
        {
            spd += (acc ? accelB : accelF) * sincos(angle);  spd -= frict * spd;
            pos += spd;  angle += turn * (rot - 1);
        }
    }

    void init()
    {
        pos = spd = {0, 0};  angle = 0;
    }

    void next(double acc, double turn)
    {
        spd += acc * sincos(angle);  spd -= frict * spd;  pos += spd;  angle += turn;
    }
};

void sweep(const Info &info, int acc, int rot, int steps, int nAcc, int nRot)
{
    Info cur = info;  cur.nextStep(acc, rot);
    for(int i = 1; i < steps; i++)
    {
        if(nAcc)
        {
            sweep(cur, acc ^ 1, rot, steps - i, nAcc - 1, nRot);
            if(nRot)
            {
                sweep(cur, acc ^ 1, (rot + 1) % 3, steps - i, nAcc - 1, nRot - 1);
                sweep(cur, acc ^ 1, (rot + 2) % 3, steps - i, nAcc - 1, nRot - 1);
            }
        }
        if(nRot)
        {
            sweep(cur, acc, (rot + 1) % 3, steps - i, nAcc, nRot - 1);
            sweep(cur, acc, (rot + 2) % 3, steps - i, nAcc, nRot - 1);
        }
        cur.nextStep(acc, rot);
    }
    cout << cur.pos.x << ' ' << cur.pos.y << endl;
}

void sweep(int steps, int nAcc, int nRot)
{
    steps /= Info::n;  Info info;  info.init();
    sweep(info, 0, 0, steps, nAcc, nRot);
    sweep(info, 0, 1, steps, nAcc, nRot);
    sweep(info, 0, 2, steps, nAcc, nRot);
    sweep(info, 1, 0, steps, nAcc, nRot);
    sweep(info, 1, 1, steps, nAcc, nRot);
    sweep(info, 1, 2, steps, nAcc, nRot);
    cout << "==========================================================" << endl;
}


void generatePoint(const Info &info, int steps, double acc, double turn)
{
    Info cur = info;
    for(int i = 0; i < steps; i++)
    {
        double stepTurn = abs(turn) > Info::turn ? copysign(Info::turn, turn) : turn;
        cur.next(acc, stepTurn);  turn -= stepTurn;
    }
    cout << cur.pos.x << ' ' << cur.pos.y << endl;
}

void generateEdge(int steps)
{
    int n = min(steps, (int)floor(pi / 2 / Info::turn + 0.001));

    Info info;  info.init();
    generatePoint(info, steps, Info::accelF, 0);
    generatePoint(info, steps, Info::accelB, 0);

    info.init();  info.next(Info::accelF, Info::turn);
    for(int i = 1; i < n; i++)
    {
        generatePoint(info, steps - i, Info::accelF, 0);
        generatePoint(info, steps - i, Info::accelB, pi / 2);
        info.next(Info::accelF, Info::turn);
    }
    info.init();  info.next(Info::accelF, -Info::turn);
    for(int i = 1; i < n; i++)
    {
        generatePoint(info, steps - i, Info::accelF, 0);
        generatePoint(info, steps - i, Info::accelB, -pi / 2);
        info.next(Info::accelF, -Info::turn);
    }
    info.init();  info.next(Info::accelB, Info::turn);
    for(int i = 1; i < n; i++)
    {
        generatePoint(info, steps - i, Info::accelB, 0);
        generatePoint(info, steps - i, Info::accelF, pi / 2);
        info.next(Info::accelB, Info::turn);
    }
    info.init();  info.next(Info::accelB, -Info::turn);
    for(int i = 1; i < n; i++)
    {
        generatePoint(info, steps - i, Info::accelB, 0);
        generatePoint(info, steps - i, Info::accelF, -pi / 2);
        info.next(Info::accelB, -Info::turn);
    }
}


int main()
{
    //sweep(128, 0, 1);
    generateEdge(64);
    return 0;
}
