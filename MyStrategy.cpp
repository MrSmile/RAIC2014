#include "MyStrategy.h"

#include <cmath>
#include <cstdlib>

#include <iostream>  // DEBUG
#include <iomanip>  // DEBUG
#include <ctime>  // DEBUG

using namespace model;
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



struct Predictor
{
    Vec2D pos, spd;
    double angle;

    void set(const Unit& unit)
    {
        pos = Vec2D(unit.getX(), unit.getY());
        spd = Vec2D(unit.getSpeedX(), unit.getSpeedY());
        angle = unit.getAngle();
    }

    void predict(double accel, double turn)
    {
        spd += accel * sincos(angle);  spd -= 0.02 * spd;  pos += spd;  angle += turn;
    }
};


Predictor hockeyist[6];

void MyStrategy::move(const Hockeyist& self, const World& world, const Game& game, Move& move)
{
    /*
    if(!self.getTeammateIndex())
    {
        auto list = world.getHockeyists();  Vec2D goalie[2];
        for(auto &hockeyist : list)if(hockeyist.getType() == GOALIE)
            goalie[hockeyist.isTeammate() ? 0 : 1] = Vec2D(hockeyist.getX(), hockeyist.getY());

        auto puck = world.getPuck();
        Vec2D puckPos(puck.getX(), puck.getY()), puckSpd(puck.getSpeedX(), puck.getSpeedY());

        cout << puckPos.x << ' ' << puckPos.y << "   ";
        cout << puckSpd.x << ' ' << puckSpd.y << "   ";
        cout << goalie[0].x << ' ' << goalie[0].y << "   ";
        cout << goalie[1].x << ' ' << goalie[1].y << endl;
    }
    */

    auto &hock = hockeyist[self.getTeammateIndex()];
    if(world.getTick())
    {
        cout << (hock.pos.x - self.getX()) << ' ' << (hock.pos.y - self.getY()) << ' ';
        cout << (hock.spd.x - self.getSpeedX()) << ' ' << (hock.spd.y - self.getSpeedY()) << ' ';
        cout << (rem(hock.angle - self.getAngle() + pi, 2 * pi) - pi) << endl;
    }

    double accel = rand() * (2.0 / RAND_MAX) - 1;
    double turn = (rand() * (2.0 / RAND_MAX) - 1) * game.getHockeyistTurnAngleFactor();
    move.setSpeedUp(accel);  move.setTurn(turn);  move.setAction(NONE);

    accel *= accel > 0 ? game.getHockeyistSpeedUpFactor() : game.getHockeyistSpeedDownFactor();
    hock.set(self);  hock.predict(accel, turn);
}

MyStrategy::MyStrategy()
{
    cout << setprecision(16);
}
