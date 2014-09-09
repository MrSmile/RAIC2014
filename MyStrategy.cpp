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
    const auto &player = world.getMyPlayer();
    if(player.isJustMissedGoal() || player.isJustScoredGoal())
    {
        move.setTurn(self.getTeammateIndex() & 1 ? -pi : pi);
        move.setSpeedUp(0);  move.setAction(NONE);  return;
    }

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
    /*
    if(world.getTick())
    {
        cout << (hock.pos.x - self.getX()) << ' ' << (hock.pos.y - self.getY()) << ' ';
        cout << (hock.spd.x - self.getSpeedX()) << ' ' << (hock.spd.y - self.getSpeedY()) << ' ';
        cout << (rem(hock.angle - self.getAngle() + pi, 2 * pi) - pi) << endl;
    }
    */

    double accel = rand() * (2.0 / RAND_MAX) - 1;
    double turn = (rand() * (2.0 / RAND_MAX) - 1) * game.getHockeyistTurnAngleFactor();
    move.setSpeedUp(accel);  move.setTurn(turn);  move.setAction(NONE);

    accel *= accel > 0 ? game.getHockeyistSpeedUpFactor() : game.getHockeyistSpeedDownFactor();
    hock.set(self);  hock.predict(accel, turn);


    if(!self.getTeammateIndex())
    {
        const auto &puck = world.getPuck();
        static double x_min = game.getRinkLeft() + puck.getRadius();
        static double x_max = game.getRinkRight() - puck.getRadius();
        static double y_min = game.getRinkTop() + puck.getRadius();
        static double y_max = game.getRinkBottom() - puck.getRadius();

        static bool fly = false;
        static Vec2D oldPos, oldSpd;
        if(puck.getOwnerHockeyistId() < 0)
        {
            Vec2D pos(puck.getX(), puck.getY()), spd(puck.getSpeedX(), puck.getSpeedY());
            if(fly)
            {
                if(oldPos.x < x_min && oldSpd.x < 0)
                {
                    double delta = oldPos.x - x_min + 0.01;
                    cout << "HIT X: " << delta << ' ' << (pos.x + 0.25 * 0.999 * oldSpd.x - x_min) << ' ' << oldSpd.x << endl;
                    oldSpd.x *= -0.25;  if(oldSpd.x < -delta)oldPos.x -= 0.8 * delta;

                }
                if(oldPos.x > x_max && oldSpd.x > 0)
                {
                    double delta = oldPos.x - x_max - 0.01;
                    cout << "HIT X: " << delta << ' ' << (pos.x + 0.25 * 0.999 * oldSpd.x - x_max) << ' ' << oldSpd.x << endl;
                    oldSpd.x *= -0.25;  if(oldSpd.x > -delta)oldPos.x -= 0.8 * delta;
                }
                if(oldPos.y < y_min && oldSpd.y < 0)
                {
                    double delta = oldPos.y - y_min + 0.01;
                    cout << "HIT Y: " << delta << ' ' << (pos.y + 0.25 * 0.999 * oldSpd.y - y_min) << ' ' << oldSpd.y << endl;
                    oldSpd.y *= -0.25;  if(oldSpd.y < -delta)oldPos.y -= 0.8 * delta;
                }
                if(oldPos.y > y_max && oldSpd.y > 0)
                {
                    double delta = oldPos.y - y_max - 0.01;
                    cout << "HIT Y: " << delta << ' ' << (pos.y + 0.25 * 0.999 * oldSpd.y - y_max) << ' ' << oldSpd.y << endl;
                    oldSpd.y *= -0.25;  if(oldSpd.y > -delta)oldPos.y -= 0.8 * delta;
                }

                oldSpd -= 0.001 * oldSpd;  oldPos += oldSpd;
                Vec2D errPos = pos - oldPos, errSpd = spd - oldSpd;
                if(abs(errPos.x) > 1e-3 || abs(errPos.y) > 1e-3 || abs(errSpd.x) > 1e-5 || abs(errSpd.y) > 1e-5)
                {
                    cout << "Error: ";
                    cout << (pos.x - oldPos.x) << ' ' << (pos.y - oldPos.y) << ' ';
                    cout << (spd.x - oldSpd.x) << ' ' << (spd.y - oldSpd.y) << endl;
                }
            }
            oldPos = pos;  oldSpd = spd;  fly = true;
        }
        else fly = false;
    }
}

MyStrategy::MyStrategy()
{
    cout << setprecision(16);
}
