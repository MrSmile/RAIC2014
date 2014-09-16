#include "MyStrategy.h"

#include <cmath>
#include <cstdlib>
#include <algorithm>
#include <limits>

#include <iostream>  // DEBUG
#include <iomanip>  // DEBUG
#include <ctime>  // DEBUG

using namespace model;
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



constexpr int maxLookahead = 256;

constexpr double hockeyistFrict = 0.02, puckFrict = 0.001, beta = -log(1 - puckFrict);
constexpr double wallBounce = 0.25, hockeyistBounce = 0.325, goalieFrict = 0.1;
constexpr double minDepth = 0.01, depthFactor = 0.8;

double rinkLeft, rinkRight, rinkTop, rinkBottom;
double goalCenter, goalHalf, hockeyistRad, puckRad, goalieSpd, goalieRange;
double accelMin, accelMax, turnAngle, maxTurnSteps;
double holdDist, strikeBase, strikeGrowth, strikeBeta;
double passBeta, passSector;  int maxSwing;

void initConsts(const Game& game, const World& world)
{
    rinkLeft = game.getRinkLeft();  rinkRight = game.getRinkRight();
    rinkTop = game.getRinkTop();  rinkBottom = game.getRinkBottom();
    goalHalf = game.getGoalNetHeight() / 2;  goalCenter = game.getGoalNetTop() + goalHalf;
    hockeyistRad = world.getHockeyists()[0].getRadius();  puckRad = world.getPuck().getRadius();
    goalieSpd = game.getGoalieMaxSpeed();  goalieRange = goalHalf - hockeyistRad;
    accelMin = -game.getHockeyistSpeedDownFactor();  accelMax = game.getHockeyistSpeedUpFactor();
    turnAngle = game.getHockeyistTurnAngleFactor();  maxTurnSteps = (pi / 2) / turnAngle;

    holdDist = game.getPuckBindingRange();
    strikeBase = 20 * game.getStrikePowerBaseFactor();
    strikeGrowth = 20 * game.getStrikePowerGrowthFactor();
    strikeBeta = 1 / (game.getStrikeAngleDeviation() * sqrt(2.0));
    passBeta = 1 / (game.getPassAngleDeviation() * sqrt(2.0));
    passSector = game.getPassSector() / 2;
    maxSwing = game.getMaxEffectiveSwingTicks();
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
    double end = (dist - rad * dir.y) * invSpd, cmp = 0;
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


struct UnitInfo
{
    Vec2D pos, spd;

    void set(const Unit& unit)
    {
        pos = Vec2D(unit.getX(), unit.getY());
        spd = Vec2D(unit.getSpeedX(), unit.getSpeedY());
    }

    bool bounceCircle(const Vec2D &center, double rad, double bounce, double frict = 0)
    {
        Vec2D delta = center - pos;  if(!(delta.sqr() < rad * rad))return false;
        double normSpd = max(0.0, spd * delta), normImp = (1 + bounce) * normSpd;
        double tanImp = frict * normImp;  tanImp = max(-tanImp, min(tanImp, spd % delta));  // TODO: incorrect tanImp

        double invLen2 = 1 / delta.sqr(), invLen = sqrt(invLen2), depth = rad * invLen - 1;
        if((normImp - normSpd) * invLen2 < depth)pos -= depthFactor * (depth - minDepth * invLen) * delta;
        spd -= (normImp * delta + tanImp * ~delta) * invLen2;  return true;
    }
};

struct HockeyistInfo : public UnitInfo
{
    double angle, angSpd;
    int cooldown;

    void set(const Hockeyist& hockeyist)
    {
        UnitInfo::set(hockeyist);
        angle = hockeyist.getAngle();  angSpd = hockeyist.getAngularSpeed();
        cooldown = hockeyist.getRemainingCooldownTicks();
    }

    void nextStep(double accel, double turn)
    {
        spd += accel * sincos(angle);  spd -= hockeyistFrict * spd;

        double delta;
        if((delta = rinkLeft + hockeyistRad - pos.x) > 0)
        {
            if(abs(pos.y - goalCenter) < goalHalf)
            {
                if(delta > hockeyistRad)return;
            }
            else
            {
                if(spd.x < 0)spd.x *= -hockeyistBounce;
                if(spd.x < delta)pos.x += depthFactor * (delta - minDepth);
            }
        }
        if((delta = rinkRight - hockeyistRad - pos.x) < 0)
        {
            double sign = pos.y - goalCenter, offs = abs(sign) - goalHalf;
            if(offs < 0)
            {
                if(delta < -1)pos.x = rinkRight - hockeyistRad - 1;
            }
            else if(!bounceCircle(Vec2D(rinkRight, goalCenter + copysign(goalHalf, sign)), hockeyistRad, hockeyistBounce) &&
                !bounceCircle(Vec2D(rinkRight, goalCenter + copysign(goalHalf + 25, sign)), hockeyistRad, hockeyistBounce / 10))
            {
                if(offs > 25)
                {
                    if(spd.x > 0)spd.x *= -hockeyistBounce;
                    if(spd.x > delta)pos.x += depthFactor * (delta + minDepth);
                }
                else
                {
                    if(spd.x > 0)spd.x *= -hockeyistBounce / 10;
                    if(spd.x > delta)pos.x += depthFactor * (delta + minDepth);
                }
            }
        }
        if((delta = rinkTop + hockeyistRad - pos.y) > 0)
        {
            if(spd.y < 0)spd.y *= -hockeyistBounce;
            if(spd.y < delta)pos.y += depthFactor * (delta - minDepth);
        }
        if((delta = rinkBottom - hockeyistRad - pos.y) < 0)
        {
            if(spd.y > 0)spd.y *= -hockeyistBounce;
            if(spd.y > delta)pos.y += depthFactor * (delta + minDepth);
        }

        pos += spd;  angSpd -= 0.0270190131 * angSpd;  angle += turn + angSpd;
    }
};

struct PuckInfo : public UnitInfo
{
    Vec2D goalie[2];
    int goalieTime;

    void set(const Puck& puck, const vector<Hockeyist> &hockeyists)
    {
        UnitInfo::set(puck);
        for(auto &hockeyist : hockeyists)if(hockeyist.getType() == GOALIE)
        {
            Vec2D cur = Vec2D(hockeyist.getX(), hockeyist.getY());
            cur.y = max(cur.y - goalieSpd, min(cur.y + goalieSpd, pos.y));
            cur.y = max(goalCenter - goalieRange, min(goalCenter + goalieRange, cur.y));
            goalie[hockeyist.isTeammate() ? 0 : 1] = cur;
        }
        goalieTime = 1000000;  // TODO
    }

    int nextStep()
    {
        spd -= puckFrict * spd;

        double delta;
        if((delta = rinkLeft + puckRad - pos.x) > 0)
        {
            if(abs(pos.y - goalCenter) < goalHalf)
            {
                if(delta > puckRad)return -1;
                //bounceCircle(Vec2D(rinkLeft, goalCenter - goalHalf), puckRad, 0.025, 0.0);
                //bounceCircle(Vec2D(rinkLeft, goalCenter + goalHalf), puckRad, 0.025, 0.0);
            }
            else
            {
                if(spd.x < 0)spd.x *= -wallBounce;
                if(spd.x < delta)pos.x += depthFactor * (delta - minDepth);
            }
        }
        if((delta = rinkRight - puckRad - pos.x) < 0)
        {
            if(abs(pos.y - goalCenter) < goalHalf)
            {
                if(delta < -puckRad)return 1;
                //bounceCircle(Vec2D(rinkRight, goalCenter - goalHalf), puckRad, 0.025, 0.0);
                //bounceCircle(Vec2D(rinkRight, goalCenter + goalHalf), puckRad, 0.025, 0.0);
            }
            else
            {
                if(spd.x > 0)spd.x *= -wallBounce;
                if(spd.x > delta)pos.x += depthFactor * (delta + minDepth);
            }
        }
        if((delta = rinkTop + puckRad - pos.y) > 0)
        {
            if(spd.y < 0)spd.y *= -wallBounce;
            if(spd.y < delta)pos.y += depthFactor * (delta - minDepth);
        }
        if((delta = rinkBottom - puckRad - pos.y) < 0)
        {
            if(spd.y > 0)spd.y *= -wallBounce;
            if(spd.y > delta)pos.y += depthFactor * (delta + minDepth);
        }
        if(goalieTime > 0)for(int i = 0; i < 2; i++)
            bounceCircle(goalie[i], hockeyistRad + puckRad, hockeyistBounce, goalieFrict);

        pos += spd;  return 0;
    }
};


double evaluateStrike(const HockeyistInfo &info, double &angle, double spd, double beta, double sector)
{
    Vec2D pos = info.pos + hockeyistFrict * info.spd + holdDist * sincos(info.angle - info.angSpd);
    //spd += info.spd * normalize(Vec2D(rinkLeft, goalCenter) - pos);
    Sector res = estimateGoalAngle(pos, spd, false);  // TODO: side
    if(!(res.halfSpan > 0))return 0;

    angle = rem(res.centerAngle - info.angle + pi, 2 * pi) - pi;
    double offs = max(0.0, angle - sector) + min(0.0, angle + sector);
    //return atan(beta * (offs + res.halfSpan)) - atan(beta * (offs - res.halfSpan));
    return erf(beta * (offs + res.halfSpan)) - erf(beta * (offs - res.halfSpan));
}

double evaluateStrike(const HockeyistInfo &info)
{
    HockeyistInfo cur = info;  double angle;
    for(int i = 0; i < maxSwing; i++)cur.nextStep(0, 0);
    return evaluateStrike(cur, angle, strikeBase + strikeGrowth * maxSwing, strikeBeta, 0);
}

double evaluatePass(const HockeyistInfo &info, double &angle)
{
    return evaluateStrike(info, angle, strikeBase, passBeta, passSector);
}

constexpr double stepValue(double val)
{
    return val < 0 ? 0 : (val > 1 ? 1 : val);
}

struct MovePlan
{
    enum Flags
    {
        FORW = 0, FIRST_LEFT  = 0, SECOND_LEFT  = 0,
        BACK = 1, FIRST_RIGHT = 2, SECOND_RIGHT = 4,
    };

    int flags, strikeTime;
    double firstTurnEnd, secondTurnStart, score, passAngle;
    bool pass;

    struct Helper
    {
        double accBase, accDelta, turn1, turn2;
        double accTime1, accTime2, turnTime1, turnTime2;

        Helper(const MovePlan &plan)
        {
            accBase = accelMin;  accDelta = accelMax;
            if(plan.flags & MovePlan::BACK)swap(accBase, accDelta);  accDelta -= accBase;
            turn1 = (plan.flags & MovePlan::FIRST_RIGHT  ? -turnAngle : turnAngle);
            turn2 = (plan.flags & MovePlan::SECOND_RIGHT ? -turnAngle : turnAngle);

            accTime2 = plan.secondTurnStart + maxTurnSteps;
            accTime1 = min(plan.firstTurnEnd - maxTurnSteps, accTime2);
            turnTime2 = plan.secondTurnStart;  turnTime1 = min(plan.firstTurnEnd, turnTime2);
        }

        double accel(int step) const
        {
            return accBase + accDelta * (stepValue(accTime2 - step) - stepValue(accTime1 - step));
        }

        double turn(int step) const
        {
            return turn1 * stepValue(turnTime1 - step) + turn2 * (1 - stepValue(turnTime2 - step));
        }
    };

    MovePlan() = default;

    MovePlan(int fl, double turn1, double turn2) : flags(fl), firstTurnEnd(turn1), secondTurnStart(turn2)
    {
    }

    MovePlan(const HockeyistInfo &info, int fl, double turn1) :
        flags(fl), firstTurnEnd(turn1), secondTurnStart(maxLookahead)
    {
        evaluate(info);  secondTurnStart = strikeTime;
    }

    MovePlan(const HockeyistInfo &info, const MovePlan &old, double step1, double step2) :
        flags(old.flags), firstTurnEnd(old.firstTurnEnd + step1), secondTurnStart(old.secondTurnStart + step2)
    {
        if(firstTurnEnd < 0)
        {
            flags ^= FIRST_RIGHT;  firstTurnEnd = -firstTurnEnd;
        }
        if(firstTurnEnd > 2 * maxTurnSteps)firstTurnEnd = 4 * maxTurnSteps - firstTurnEnd;

        if(secondTurnStart < 0)secondTurnStart = -secondTurnStart;
        if(secondTurnStart > old.strikeTime)
        {
            flags ^= SECOND_RIGHT;  secondTurnStart = 2 * old.strikeTime - secondTurnStart;
        }

        evaluate(info);
    }

    bool operator < (const MovePlan &plan) const
    {
        return score < plan.score;
    }

    void evaluate(const HockeyistInfo &info)
    {
        HockeyistInfo cur = info;  Helper helper(*this);
        score = -numeric_limits<double>::infinity();  strikeTime = -1;
        int n = min(maxLookahead, (int)ceil(helper.turnTime2 + maxTurnSteps));
        for(int i = 0; i < n; i++)
        {
            cur.nextStep(helper.accel(i), helper.turn(i));
            if(cur.pos.x < rinkLeft + 100)return;
            if(info.cooldown > i)continue;

            double val = evaluateStrike(cur), angle = 0;
            if(val > score)
            {
                score = val;  strikeTime = i + 1;  pass = false;
            }
            val = evaluatePass(cur, angle);
            if(val > score)
            {
                score = val;  strikeTime = i + 1;  passAngle = angle;  pass = true;
            }
        }
    }

    void execute(Move &move)
    {
        Helper helper(*this);  double acc = helper.accel(0);
        move.setSpeedUp(acc > 0 ? acc / accelMax : acc / accelMin);
        move.setTurn(helper.turn(0));  move.setAction(NONE);

        cout << "Best score: " << score << ", strike: " << strikeTime;
        if(pass)cout << ", pass: " << passAngle;  cout << endl;

        firstTurnEnd = max(0.0, firstTurnEnd - 1);
        secondTurnStart = max(0.0, secondTurnStart - 1);
    }
};

void optimizeMove(const HockeyistInfo &info, vector<MovePlan> &moves, double step)
{
    constexpr int survive = 4, offspring = 4;

    double step2 = 4 * step;  int n = 0;
    auto rend = moves.rend();  make_heap(moves.rbegin(), rend);
    for(; n < survive && rend != moves.rbegin(); n++, --rend)
        pop_heap(moves.rbegin(), rend);  moves.resize(n);
    for(int i = 0; i < n; i++)for(int j = 0; j < offspring; j++)
        moves.emplace_back(info, moves[i],
            step  * (rand() * (2.0 / RAND_MAX) - 1),
            step2 * (rand() * (2.0 / RAND_MAX) - 1));
}

MovePlan findBestMove(const HockeyistInfo &info, const MovePlan &old)
{
    constexpr int step = 4;

    vector<MovePlan> moves;  moves.push_back(old);
    moves.emplace_back(info, MovePlan::FORW, 0);
    moves.emplace_back(info, MovePlan::BACK, 0);
    int n = min(maxLookahead, (int)ceil(2 * maxTurnSteps));
    for(int i = step; i <= n; i += step)
    {
        moves.emplace_back(info, MovePlan::FORW | MovePlan::FIRST_LEFT,  i);
        moves.emplace_back(info, MovePlan::FORW | MovePlan::FIRST_RIGHT, i);
        moves.emplace_back(info, MovePlan::BACK | MovePlan::FIRST_LEFT,  i);
        moves.emplace_back(info, MovePlan::BACK | MovePlan::FIRST_RIGHT, i);
    }
    optimizeMove(info, moves, 4);
    optimizeMove(info, moves, 2);
    optimizeMove(info, moves, 1);
    optimizeMove(info, moves, 0.5);
    optimizeMove(info, moves, 0.25);
    optimizeMove(info, moves, 0.125);

    const MovePlan *res = &old;
    double best = -numeric_limits<double>::infinity();
    for(auto &move : moves)if(move.score > best)
    {
        res = &move;  best = move.score;
    }
    return *res;
}


int state[6] = {0, 0, 0, 0, 0, 0};
HockeyistInfo hockeyistInfo[6];
PuckInfo puckInfo;

void MyStrategy::move(const Hockeyist& self, const World& world, const Game& game, Move& move)
{
    const auto &player = world.getMyPlayer();
    if(player.isJustMissedGoal() || player.isJustScoredGoal())
    {
        move.setTurn(self.getTeammateIndex() & 1 ? -pi : pi);
        move.setSpeedUp(0);  move.setAction(NONE);  return;
    }

    /*
    if(!world.getTick() && self.getTeammateIndex())
    {
        initConsts(game, world);  srand(time(0));
    }
    if(!(rand() % 57))state[self.getTeammateIndex()] = rand() % 6;
    switch(state[self.getTeammateIndex()])
    {
    case 0:  move.setSpeedUp(+1);  move.setTurn(0);    break;
    case 1:  move.setSpeedUp(+1);  move.setTurn(+pi);  break;
    case 2:  move.setSpeedUp(-1);  move.setTurn(-pi);  break;
    case 3:  move.setSpeedUp(-1);  move.setTurn(0);    break;
    case 4:  move.setSpeedUp(-1);  move.setTurn(+pi);  break;
    case 5:  move.setSpeedUp(+1);  move.setTurn(-pi);  break;
    }
    move.setAction(NONE);
    */

    /*
    if(!self.getTeammateIndex())
    {
        static bool hold = false;
        static Vec2D oldPos, oldPuck, oldSpd;
        static double oldAngle, oldAngSpd;
        static int oldSwing;

        const auto &puck = world.getPuck();
        if(puck.getOwnerHockeyistId() >= 0)
        {
            const auto &hockeyists = world.getHockeyists();
            for(auto &hockeyist : hockeyists)if(hockeyist.getId() == puck.getOwnerHockeyistId())
            {
                oldPos = {hockeyist.getX(), hockeyist.getY()};
                oldSpd = {hockeyist.getSpeedX(), hockeyist.getSpeedY()};
                oldAngle = hockeyist.getAngle();  oldAngSpd = hockeyist.getAngularSpeed();
                oldSwing = hockeyist.getSwingTicks();  break;
            }
            oldPuck = {puck.getX(), puck.getY()};  hold = true;
        }
        else if(hold)
        {
            cout << "Strike: ";
            cout << oldPos.x << ' ' << oldPos.y << ' ';
            cout << oldPuck.x << ' ' << oldPuck.y << ' ';
            cout << oldSpd.x << ' ' << oldSpd.y << ' ';
            cout << oldAngle << ' ' << oldAngSpd << ' ' << oldSwing << ' ';
            cout << puck.getX() << ' ' << puck.getY() << ' ';
            cout << puck.getSpeedX() / 0.999 << ' ' << puck.getSpeedY() / 0.999 << endl;
            hold = false;
        }
    }
    move.setSpeedUp(0);  move.setTurn(0);  move.setAction(NONE);
    */

    //static Vec2D targetPos;
    //static double targetAngle;
    const auto &puck = world.getPuck();
    if(self.getTeammateIndex())
    {
        move.setSpeedUp(0);  move.setTurn(0);  move.setAction(NONE);
    }
    else if(self.getSwingTicks())
    {
        move.setSpeedUp(0);  move.setTurn(0);
        //move.setAction(rand() % (maxSwing - self.getSwingTicks()) ? SWING : STRIKE);
        move.setAction(self.getSwingTicks() < maxSwing ? SWING : STRIKE);
    }
    else if(puck.getOwnerHockeyistId() != self.getId())
    {
        if(!world.getTick())
        {
            initConsts(game, world);  srand(time(0));
            /*
            for(;;)
            {
                double x = rinkLeft + (rinkRight - rinkLeft) * rand() / RAND_MAX;
                double y = rinkTop  + (rinkBottom - rinkTop) * rand() / RAND_MAX;
                Sector target = estimateGoalAngle(Vec2D(x, y), 15, false);
                if(target.halfSpan < pi / 45)continue;

                targetAngle = target.centerAngle;
                targetPos = Vec2D(x, y) - 55 * sincos(targetAngle);
                cout << "Target: " << targetPos.x << ' ' << targetPos.y << ' ';
                cout << targetAngle * (180 / pi) << ' ' << target.halfSpan * (180 / pi) << endl;
                break;
            }
            */
            //targetPos = Vec2D(150, 300);  targetAngle = 0;
        }
        Vec2D delta(puck.getX() - self.getX(), puck.getY() - self.getY());
        double angle = rem(atan2(delta.y, delta.x) - self.getAngle() + pi, 2 * pi) - pi;
        move.setSpeedUp(abs(angle) < pi / 4 ? 1 : 0);  move.setTurn(angle);  move.setAction(TAKE_PUCK);
    }
    else
    {
        /*
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
        */

        HockeyistInfo info;  info.set(self);
        static MovePlan plan(MovePlan::FORW, 0, maxLookahead);
        plan.evaluate(info);  plan = findBestMove(info, plan);
        double angle = 0, strike = evaluateStrike(info), pass = evaluatePass(info, angle);
        if(plan.score < max(strike, pass))
        {
            move.setSpeedUp(0);  move.setTurn(0);
            if(pass > strike)
            {
                move.setPassPower(1);  move.setPassAngle(angle);  move.setAction(PASS);
            }
            else move.setAction(SWING);  plan = MovePlan(MovePlan::FORW, 0, maxLookahead);
        }
        else plan.execute(move);
    }

#if 0

    static Vec2D target[2];
    static int mode = 0, wait = 0;
    if(self.getTeammateIndex() != 1)
    {
        //move.setSpeedUp(0);  move.setTurn(0);  move.setAction(NONE);
    }
    else if(world.getTick() < wait)
    {
        move.setSpeedUp(0);  move.setTurn(0);  move.setAction(NONE);
    }
    else
    {
        if(!world.getTick())
        {
            initConsts(game, world);  srand(time(0));
            target[0] = Vec2D(rinkRight, goalCenter + goalHalf);
        }
        Vec2D delta = target[mode] - Vec2D(self.getX(), self.getY());
        if(delta.sqr() < sqr(60))
        {
            move.setSpeedUp(0);  move.setTurn(0);  move.setAction(NONE);
            if(mode)target[0] = Vec2D(rinkRight, 530 + rand() * (60.0 / RAND_MAX));
            //if(mode)target[0] = Vec2D(rinkRight, goalCenter + 0.8 * goalHalf + rand() * (0.6 * goalHalf / RAND_MAX));
            else target[1] = Vec2D(155 + double(869 - 155) * rand() / RAND_MAX, 350 + double(690 - 350) * rand() / RAND_MAX);
            mode = 1 - mode;  wait = world.getTick() + 30;
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
                double acc = 1;
                move.setSpeedUp(angle > 0 ? -acc : acc);  move.setTurn(0);  move.setAction(NONE);
            }
        }
    }

    auto &info = hockeyistInfo[self.getTeammateIndex()];
    if(world.getTick())
    {
        Vec2D pos(self.getX(), self.getY()), spd(self.getSpeedX(), self.getSpeedY());

        Vec2D errPos = pos - info.pos, errSpd = spd - info.spd;
        double errAngle = rem(info.angle - self.getAngle() + pi, 2 * pi) - pi;
        if(abs(errPos.x) > 1e-3 || abs(errPos.y) > 1e-3 ||
            abs(errSpd.x) > 1e-4 || abs(errSpd.y) > 1e-4 || abs(errAngle) > 1e-5)
        {
            cout << info.hitType << "Error: ";
            cout << errPos.x << ' ' << errPos.y << ' ';
            cout << errSpd.x << ' ' << errSpd.y << ' ';
            cout << errAngle << endl;
        }
    }

    double accel = max(-1.0, min(1.0, move.getSpeedUp()));
    accel *= accel > 0 ? game.getHockeyistSpeedUpFactor() : game.getHockeyistSpeedDownFactor();
    double turn = max(-game.getHockeyistTurnAngleFactor(), min(game.getHockeyistTurnAngleFactor(), move.getTurn()));
    info.set(self);  info.nextStep(accel, turn);

#endif

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
                if(puckInfo.nextStep())cout << "GOAL!!!" << endl;
                Vec2D errPos = pos - puckInfo.pos, errSpd = spd - puckInfo.spd;
                if(abs(errPos.x) > 1e-3 || abs(errPos.y) > 1e-3 || abs(errSpd.x) > 1e-4 || abs(errSpd.y) > 1e-4)
                {
                    cout << puckInfo.hitType << "Error: ";
                    cout << errPos.x << ' ' << errPos.y << ' ';
                    cout << errSpd.x << ' ' << errSpd.y << endl;
                }
            }
            const auto &hockeyists = world.getHockeyists();
            puckInfo.set(puck, hockeyists);  fly = true;
        }
        else fly = false;
    }
    */
}

MyStrategy::MyStrategy()
{
    cout << setprecision(16);
}
