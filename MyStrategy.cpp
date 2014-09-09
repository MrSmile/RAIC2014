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



double rinkLeft, rinkRight, rinkTop, rinkBottom;
double goalCenter, goalHalf, hockeyistRad, puckRad, goalieSpd;
constexpr double hockeyistFrict = 0.02, puckFrict = 0.001;
constexpr double wallCoeff = 0.25, goalieCoeff = 0.325, goalieFrict = 0.1;
constexpr double minDepth = 0.01, depthFactor = 0.8;

void initConsts(const Game& game, const World& world)
{
    rinkLeft = game.getRinkLeft();  rinkRight = game.getRinkRight();
    rinkTop = game.getRinkTop();  rinkBottom = game.getRinkBottom();
    goalHalf = game.getGoalNetHeight() / 2;  goalCenter = game.getGoalNetTop() + goalHalf;
    hockeyistRad = world.getHockeyists()[0].getRadius();  puckRad = world.getPuck().getRadius();
    goalieSpd = game.getGoalieMaxSpeed();
}

struct HockeyistPredictor
{
    Vec2D pos, spd;
    double angle;

    void set(const Hockeyist& hockeyist)
    {
        pos = Vec2D(hockeyist.getX(), hockeyist.getY());
        spd = Vec2D(hockeyist.getSpeedX(), hockeyist.getSpeedY());
        angle = hockeyist.getAngle();
    }

    void predict(double accel, double turn)
    {
        spd += accel * sincos(angle);  spd -= hockeyistFrict * spd;  pos += spd;  angle += turn;
    }
};

struct PuckPredictor
{
    string hitType;  // DEBUG
    Vec2D pos, spd, goalie[2];
    int goalieTime;

    void set(const Puck& puck, const vector<Hockeyist> &hockeyists)
    {
        pos = Vec2D(puck.getX(), puck.getY());
        spd = Vec2D(puck.getSpeedX(), puck.getSpeedY());

        for(auto &hockeyist : hockeyists)if(hockeyist.getType() == GOALIE)
        {
            Vec2D cur = Vec2D(hockeyist.getX(), hockeyist.getY());
            cur.y = max(cur.y - goalieSpd, min(cur.y + goalieSpd, pos.y));
            cur.y = max(goalCenter - goalHalf + hockeyistRad, min(goalCenter + goalHalf - hockeyistRad, cur.y));
            goalie[hockeyist.isTeammate() ? 0 : 1] = cur;
        }
        goalieTime = 1000000;  // TODO
    }

    void bounceCircle(const Vec2D &center, double rad)
    {
        Vec2D delta = center - pos;  if(!(delta.sqr() < rad * rad))return;
        double normSpd = max(0.0, spd * delta), normImp = (1 + goalieCoeff) * normSpd;
        double tanImp = goalieFrict * normImp;  tanImp = max(-tanImp, min(tanImp, spd % delta));

        double invLen2 = 1 / delta.sqr();
        double depth = ((rad - minDepth) * sqrt(invLen2) - 1);
        if((normImp - normSpd) * invLen2 < depth)pos -= depthFactor * depth * delta;
        spd -= (normImp * delta + tanImp * ~delta) * invLen2;

        hitType += "CIRCLE ";
    }

    int predict()
    {
        hitType.clear();

        double delta;
        if((delta = rinkLeft + puckRad - pos.x) > 0)
        {
            if(abs(pos.y - goalCenter) < goalHalf)
            {
                if(delta > puckRad)return -1;
                cout << "CIRCLES!!!" << endl;
                bounceCircle(Vec2D(rinkLeft, goalCenter - goalHalf), puckRad);
                bounceCircle(Vec2D(rinkLeft, goalCenter + goalHalf), puckRad);
            }
            else
            {
                delta -= minDepth;  if(spd.x < 0)spd.x *= -wallCoeff;
                if(spd.x < delta)pos.x += depthFactor * delta;
                hitType += "LEFT ";
            }
        }
        if((delta = rinkRight - puckRad - pos.x) < 0)
        {
            if(abs(pos.y - goalCenter) < goalHalf)
            {
                if(delta < -puckRad)return 1;
                cout << "CIRCLES!!!" << endl;
                bounceCircle(Vec2D(rinkRight, goalCenter - goalHalf), puckRad);
                bounceCircle(Vec2D(rinkRight, goalCenter + goalHalf), puckRad);
            }
            else
            {
                delta += minDepth;  if(spd.x > 0)spd.x *= -wallCoeff;
                if(spd.x > delta)pos.x += depthFactor * delta;
                hitType += "RIGHT ";
            }
        }
        if((delta = rinkTop + puckRad - pos.y) > 0)
        {
            delta -= minDepth;  if(spd.y < 0)spd.y *= -wallCoeff;
            if(spd.y < delta)pos.y += depthFactor * delta;
            hitType += "TOP ";
        }
        if((delta = rinkBottom - puckRad - pos.y) < 0)
        {
            delta += minDepth;  if(spd.y > 0)spd.y *= -wallCoeff;
            if(spd.y > delta)pos.y += depthFactor * delta;
            hitType += "BOTTOM ";
        }
        if(goalieTime > 0)for(int i = 0; i < 2; i++)bounceCircle(goalie[i], hockeyistRad + puckRad);
        spd -= puckFrict * spd;  pos += spd;  return 0;
    }
};


HockeyistPredictor hockeyistPred[6];
PuckPredictor puckPred;

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

    auto &hock = hockeyistPred[self.getTeammateIndex()];
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
        if(!world.getTick())initConsts(game, world);

        static bool fly = false;
        const auto &puck = world.getPuck();
        if(puck.getOwnerHockeyistId() < 0)
        {
            Vec2D pos(puck.getX(), puck.getY()), spd(puck.getSpeedX(), puck.getSpeedY());
            if(fly)
            {
                if(puckPred.predict())cout << "GOAL!!!" << endl;
                Vec2D errPos = pos - puckPred.pos, errSpd = spd - puckPred.spd;
                if(abs(errPos.x) > 1e-3 || abs(errPos.y) > 1e-3 || abs(errSpd.x) > 1e-4 || abs(errSpd.y) > 1e-4)
                {
                    cout << puckPred.hitType << "Error: ";
                    cout << errPos.x << ' ' << errPos.y << ' ';
                    cout << errSpd.x << ' ' << errSpd.y << endl;
                }
            }
            const auto &hockeyists = world.getHockeyists();
            puckPred.set(puck, hockeyists);  fly = true;
        }
        else fly = false;
    }
}

MyStrategy::MyStrategy()
{
    cout << setprecision(16);
}
