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
    const auto &player = world.getMyPlayer();
    if(player.isJustMissedGoal() || player.isJustScoredGoal())
    {
        move.setTurn(self.getTeammateIndex() & 1 ? -pi : pi);
        move.setSpeedUp(0);  move.setAction(NONE);  return;
    }
    */

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

    double accel = 0;//rand() * (2.0 / RAND_MAX) - 1;
    double turn = 0;//(rand() * (2.0 / RAND_MAX) - 1) * game.getHockeyistTurnAngleFactor();
    move.setSpeedUp(accel);  move.setTurn(turn);  move.setAction(NONE);

    accel *= accel > 0 ? game.getHockeyistSpeedUpFactor() : game.getHockeyistSpeedDownFactor();
    hock.set(self);  hock.predict(accel, turn);


    if(!self.getTeammateIndex())
    {
        const auto &puck = world.getPuck();
        static const double rinkLeft = game.getRinkLeft() + puck.getRadius();
        static const double rinkRight = game.getRinkRight() - puck.getRadius();
        static const double rinkTop = game.getRinkTop() + puck.getRadius();
        static const double rinkBottom = game.getRinkBottom() - puck.getRadius();
        static const double goalTop = game.getGoalNetTop() + self.getRadius();
        static const double goalBottom = game.getGoalNetTop() + game.getGoalNetHeight() - self.getRadius();
        static const double goalieSpd = game.getGoalieMaxSpeed();
        static const double rad = self.getRadius() + puck.getRadius();
        constexpr double /*hockeyistFrict = 0.02,*/ puckFrict = 0.001;
        constexpr double wallCoeff = 0.25, goalieCoeff = 0.325, goalieFrict = 0.1;
        constexpr double minDepth = 0.01, depthFactor = 0.8;

        static bool fly = false;
        static Vec2D oldPos, oldSpd, goalie[2];
        if(puck.getOwnerHockeyistId() < 0)
        {
            Vec2D pos(puck.getX(), puck.getY()), spd(puck.getSpeedX(), puck.getSpeedY());
            if(fly)
            {
                bool flag = false;
                if(oldPos.x < rinkLeft && oldSpd.x < 0)
                {
                    double delta = oldPos.x - rinkLeft + minDepth;
                    oldSpd.x *= -wallCoeff;  if(oldSpd.x < -delta)oldPos.x -= depthFactor * delta;
                    //if(!(oldSpd.x < -delta))flag = true;

                }
                if(oldPos.x > rinkRight && oldSpd.x > 0)
                {
                    double delta = oldPos.x - rinkRight - minDepth;
                    oldSpd.x *= -wallCoeff;  if(oldSpd.x > -delta)oldPos.x -= depthFactor * delta;
                    //if(!(oldSpd.x > -delta))flag = true;
                }
                if(oldPos.y < rinkTop && oldSpd.y < 0)
                {
                    double delta = oldPos.y - rinkTop + minDepth;
                    oldSpd.y *= -wallCoeff;  if(oldSpd.y < -delta)oldPos.y -= depthFactor * delta;
                    //if(!(oldSpd.y < -delta))flag = true;
                }
                if(oldPos.y > rinkBottom && oldSpd.y > 0)
                {
                    double delta = oldPos.y - rinkBottom - minDepth;
                    oldSpd.y *= -wallCoeff;  if(oldSpd.y > -delta)oldPos.y -= depthFactor * delta;
                    //if(!(oldSpd.y > -delta))flag = true;
                }
                for(int i = 0; i < 2; i++)
                {
                    Vec2D delta(goalie[i] - oldPos);  if(!(delta.sqr() < rad * rad))continue;
                    double dot = oldSpd * delta;  if(!(dot > 0))continue;
                    double impDot = (1 + goalieCoeff) * dot, impCross = goalieFrict * impDot;
                    impCross = max(-impCross, min(impCross, oldSpd % delta));

                    double invLen2 = 1 / delta.sqr();
                    double depth = ((rad - minDepth) * sqrt(invLen2) - 1);
                    if((impDot - dot) * invLen2 < depth)oldPos -= depthFactor * depth * delta;
                    oldSpd -= (impDot * delta + impCross * ~delta) * invLen2;
                    //if(!((impDot - dot) * invLen2 < depth))flag = true;
                }
                if(flag)cout << "!!! ";

                oldSpd -= puckFrict * oldSpd;  oldPos += oldSpd;
                Vec2D errPos = pos - oldPos, errSpd = spd - oldSpd;
                if(abs(errPos.x) > 1e-3 || abs(errPos.y) > 1e-3 || abs(errSpd.x) > 1e-4 || abs(errSpd.y) > 1e-4)
                {
                    cout << "Error: ";
                    cout << (pos.x - oldPos.x) << ' ' << (pos.y - oldPos.y) << ' ';
                    cout << (spd.x - oldSpd.x) << ' ' << (spd.y - oldSpd.y) << endl;
                }
                else if(flag)cout << "OK" << endl;
            }
            oldPos = pos;  oldSpd = spd;  fly = true;

            const auto &hockeyists = world.getHockeyists();
            for(auto &hockeyist : hockeyists)if(hockeyist.getType() == GOALIE)
            {
                Vec2D cur = Vec2D(hockeyist.getX(), hockeyist.getY());
                cur.y = max(cur.y - goalieSpd, min(cur.y + goalieSpd, pos.y));
                cur.y = max(goalTop, min(goalBottom, cur.y));
                goalie[hockeyist.isTeammate() ? 0 : 1] = cur;
            }
        }
        else fly = false;
    }
}

MyStrategy::MyStrategy()
{
    cout << setprecision(16);
}
