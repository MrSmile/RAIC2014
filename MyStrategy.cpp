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
double goalCenter, goalHalf, hockeyistRad, puckRad, goalieSpd, goalieRange;
constexpr double hockeyistFrict = 0.02, puckFrict = 0.001, beta = -log(1 - puckFrict);
constexpr double wallCoeff = 0.25, goalieCoeff = 0.325, goalieFrict = 0.1;
constexpr double minDepth = 0.01, depthFactor = 0.8;

void initConsts(const Game& game, const World& world)
{
    rinkLeft = game.getRinkLeft();  rinkRight = game.getRinkRight();
    rinkTop = game.getRinkTop();  rinkBottom = game.getRinkBottom();
    goalHalf = game.getGoalNetHeight() / 2;  goalCenter = game.getGoalNetTop() + goalHalf;
    hockeyistRad = world.getHockeyists()[0].getRadius();  puckRad = world.getPuck().getRadius();
    goalieSpd = game.getGoalieMaxSpeed();  goalieRange = goalHalf - hockeyistRad;
}


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
            cur.y = max(goalCenter - goalieRange, min(goalCenter + goalieRange, cur.y));
            goalie[hockeyist.isTeammate() ? 0 : 1] = cur;
        }
        goalieTime = 1000000;  // TODO
    }

    void bounceCircle(const Vec2D &center, double rad, double coeff, double frict)
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
                bounceCircle(Vec2D(rinkLeft, goalCenter - goalHalf), puckRad, 0.025, 0.0);
                bounceCircle(Vec2D(rinkLeft, goalCenter + goalHalf), puckRad, 0.025, 0.0);
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
                bounceCircle(Vec2D(rinkRight, goalCenter - goalHalf), puckRad, 0.025, 0.0);
                bounceCircle(Vec2D(rinkRight, goalCenter + goalHalf), puckRad, 0.025, 0.0);
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
        if(goalieTime > 0)for(int i = 0; i < 2; i++)
            bounceCircle(goalie[i], hockeyistRad + puckRad, goalieCoeff, goalieFrict);
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
    move.setSpeedUp(rand() * (2.0 / RAND_MAX) - 1);
    move.setTurn((rand() * (2.0 / RAND_MAX) - 1) * game.getHockeyistTurnAngleFactor());
    move.setAction(NONE);
    */

    static Vec2D targetPos;
    static double targetAngle;
    const auto &puck = world.getPuck();
    if(self.getTeammateIndex())
    {
        move.setSpeedUp(0);  move.setTurn(0);  move.setAction(NONE);
    }
    else if(self.getSwingTicks())
    {
        move.setSpeedUp(0);  move.setTurn(0);
        //move.setAction(rand() % (game.getMaxEffectiveSwingTicks() - self.getSwingTicks()) ? SWING : STRIKE);
        move.setAction(self.getSwingTicks() < game.getMaxEffectiveSwingTicks() ? SWING : STRIKE);
    }
    else if(puck.getOwnerHockeyistId() != self.getId())
    {
        if(!world.getTick())
        {
            initConsts(game, world);  srand(time(0));
            for(;;)
            {
                double x = rinkLeft + (rinkRight - rinkLeft) * rand() / RAND_MAX;
                double y = rinkTop  + (rinkBottom - rinkTop) * rand() / RAND_MAX;
                Sector target = estimateGoalAngle(Vec2D(x, y), 20, false);
                if(target.halfSpan < pi / 45)continue;

                targetAngle = target.centerAngle;
                targetPos = Vec2D(x, y) - 55 * sincos(targetAngle);
                cout << "Target: " << targetPos.x << ' ' << targetPos.y << ' ';
                cout << targetAngle * (180 / pi) << ' ' << target.halfSpan * (180 / pi) << endl;
                break;
            }
        }
        Vec2D delta(puck.getX() - self.getX(), puck.getY() - self.getY());
        double angle = rem(atan2(delta.y, delta.x) - self.getAngle() + pi, 2 * pi) - pi;
        move.setSpeedUp(abs(angle) < pi / 4 ? 1 : 0);  move.setTurn(angle);  move.setAction(TAKE_PUCK);
    }
    else
    {
        Vec2D delta = targetPos - Vec2D(self.getX(), self.getY());
        if(delta.sqr() < 1)
        {
            double angle = rem(targetAngle - self.getAngle() + pi, 2 * pi) - pi;
            if(abs(angle) > 1e-4)
            {
                move.setSpeedUp(0);  move.setTurn(angle);  move.setAction(NONE);
            }
            else
            {
                move.setSpeedUp(0);  move.setTurn(0);  move.setAction(SWING);
            }
        }
        else
        {
            double angle = rem(atan2(delta.y, delta.x) - self.getAngle() + pi / 2, 2 * pi) - pi;
            if(abs(abs(angle) - pi / 2) > 1e-4)
            {
                move.setSpeedUp(0);  move.setTurn(pi / 2 - abs(angle));  move.setAction(NONE);
            }
            else
            {
                double acc = min(1.0, delta.len() / 1e3);
                move.setSpeedUp(angle > 0 ? -acc : acc);  move.setTurn(0);  move.setAction(NONE);
            }
        }
    }

    auto &pred = hockeyistPred[self.getTeammateIndex()];
    if(world.getTick())
    {
        Vec2D pos(self.getX(), self.getY()), spd(self.getSpeedX(), self.getSpeedY());

        Vec2D errPos = pos - pred.pos, errSpd = spd - pred.spd;
        double errAngle = rem(pred.angle - self.getAngle() + pi, 2 * pi) - pi;
        if(abs(errPos.x) > 1e-3 || abs(errPos.y) > 1e-3 ||
            abs(errSpd.x) > 1e-4 || abs(errSpd.y) > 1e-4 || abs(errAngle) > 1e-5)
        {
            cout << "Error: " << errPos.x << ' ' << errPos.y << ' ';
            cout << errSpd.x << ' ' << errSpd.y << ' ' << errAngle << endl;
        }
    }

    double accel = max(-1.0, min(1.0, move.getSpeedUp()));
    accel *= accel > 0 ? game.getHockeyistSpeedUpFactor() : game.getHockeyistSpeedDownFactor();
    double turn = max(-game.getHockeyistTurnAngleFactor(), min(game.getHockeyistTurnAngleFactor(), move.getTurn()));
    pred.set(self);  pred.predict(accel, turn);

    /*
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
    */
}

MyStrategy::MyStrategy()
{
    cout << setprecision(16);
}
