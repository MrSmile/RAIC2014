#include "MyStrategy.h"

#include <cmath>
#include <cstdlib>
#include <cassert>
#include <algorithm>
#include <limits>
#include <map>

#include <iostream>  // DEBUG
#include <iomanip>  // DEBUG

using namespace model;
using namespace std;

constexpr double pi = 3.14159265358979323846264338327950288;


inline double frand()
{
    return rand() * (2.0 / RAND_MAX) - 1;
}


inline constexpr double sqr(double x)
{
    return x * x;
}

inline constexpr double rem(double x, double y)
{
    return y * (x / y - floor(x / y));
}

inline constexpr double relAngle(double angle)
{
    return rem(angle + pi, 2 * pi) - pi;
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

inline constexpr Vec2D sincos_fast(float angle)
{
    return Vec2D(cos(angle), sin(angle));
}

inline constexpr Vec2D conj(const Vec2D &v)
{
    return Vec2D(v.x, -v.y);
}




constexpr int maxLookahead = 256;
constexpr double timeGamma = 1 - 1.0 / maxLookahead;
constexpr double cellSize = 48;
constexpr int gridBorder = 2;

constexpr double hockeyistFrict = 0.02, puckFrict = 0.001, puckBeta = -log(1 - puckFrict);
constexpr double wallBounce = 0.25, hockeyistBounce = 0.325, goalieFrict = 0.1;
constexpr double minDepth = 0.01, depthFactor = 0.8;

double rinkLeft, rinkRight, rinkTop, rinkBottom;
double rinkWidth, rinkHeight;  Vec2D rinkCenter;
int gridHalfWidth, gridHalfHeight, gridLine, gridHeight, gridCenter, gridSize;
double goalCenter, goalHalf, hockeyistRad, puckRad, goalieSpd, goalieRange;
double accelMin, accelMax, turnAngle, maxTurnSteps;
double stickLength, stickSectorTan, holdDist, timeGammaSwing;  int maxSwing;
double strikeBase, strikeGrowth, passBase, strikeBeta, passBeta, passSector;
double minChance, maxChance, pickChance, strikeChance, chanceDrop;
double knockdownChance, takeAwayChance;

int globalTick = -1;
int defendCell, attackCell[2];
bool leftPlayer;


void initConsts(const Game& game, const World& world)
{
    rinkLeft = game.getRinkLeft();  rinkRight = game.getRinkRight();
    rinkTop = game.getRinkTop();  rinkBottom = game.getRinkBottom();
    rinkWidth = rinkRight - rinkLeft;  rinkHeight = rinkBottom - rinkTop;
    rinkCenter.x = (rinkLeft + rinkRight) / 2;
    rinkCenter.y = (rinkTop + rinkBottom) / 2;

    gridHalfWidth  = lround(rinkWidth  / (2 * cellSize)) + gridBorder + 1;
    gridHalfHeight = lround(rinkHeight / (2 * cellSize)) + gridBorder + 1;
    gridLine = 2 * gridHalfWidth + 1;  gridHeight = 2 * gridHalfHeight + 1;
    gridCenter = gridHalfHeight * gridLine + gridHalfWidth;
    gridSize = gridHeight * gridLine;

    goalHalf = game.getGoalNetHeight() / 2;  goalCenter = game.getGoalNetTop() + goalHalf;
    hockeyistRad = world.getHockeyists()[0].getRadius();  puckRad = world.getPuck().getRadius();
    goalieSpd = game.getGoalieMaxSpeed();  goalieRange = goalHalf - hockeyistRad;
    accelMin = -game.getHockeyistSpeedDownFactor();  accelMax = game.getHockeyistSpeedUpFactor();
    turnAngle = game.getHockeyistTurnAngleFactor();  maxTurnSteps = (pi / 2) / turnAngle;

    stickLength = game.getStickLength();
    stickSectorTan = tan(game.getStickSector() / 2);
    holdDist = game.getPuckBindingRange();
    maxSwing = game.getMaxEffectiveSwingTicks();
    timeGammaSwing = pow(timeGamma, maxSwing);
    strikeBase   = game.getStruckPuckInitialSpeedFactor() * game.getStrikePowerBaseFactor();
    strikeGrowth = game.getStruckPuckInitialSpeedFactor() * game.getStrikePowerGrowthFactor();
    passBase     = game.getStruckPuckInitialSpeedFactor() * game.getPassPowerFactor();
    strikeBeta = 1 / (game.getStrikeAngleDeviation() * sqrt(2.0));
    passBeta   = 1 / (game.getPassAngleDeviation()   * sqrt(2.0));
    passSector = game.getPassSector() / 2;

    minChance = game.getMinActionChance();
    maxChance = game.getMaxActionChance();
    pickChance   = 1 + game.getPickUpPuckBaseChance();
    strikeChance = 1 + game.getStrikePuckBaseChance();
    chanceDrop = 1 / game.getStruckPuckInitialSpeedFactor();
    knockdownChance = game.getKnockdownChanceFactor();


    leftPlayer = (2 * world.getMyPlayer().getNetBack() < rinkLeft + rinkRight);
    int defHOffs = lround(0.40 * rinkWidth  / cellSize);
    int attHOffs = lround(0.25 * rinkWidth  / cellSize);
    int attVOffs = lround(0.25 * rinkHeight / cellSize);
    if(leftPlayer)
    {
        defHOffs = -defHOffs;  attHOffs = -attHOffs;
    }
    defendCell = gridCenter + defHOffs;
    attackCell[0] = gridCenter - attHOffs - attVOffs * gridLine;
    attackCell[1] = gridCenter - attHOffs + attVOffs * gridLine;
}

inline int gridPos(const Vec2D &pos)
{
    Vec2D offs = (pos - rinkCenter) / cellSize;
    int x = lround(offs.x), y = lround(offs.y);
    return gridCenter + y * gridLine + x;
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
    Vec2D dir;
    double angle, angSpd;
    int knockdown, cooldown;

    void set(const Hockeyist& hockeyist)
    {
        UnitInfo::set(hockeyist);
        dir = sincos_fast(angle = hockeyist.getAngle());
        angSpd = hockeyist.getAngularSpeed();

        knockdown = hockeyist.getRemainingKnockdownTicks();
        cooldown = max(knockdown, hockeyist.getRemainingCooldownTicks());
    }

    void nextStep(double accel, double turn, int time = -1)
    {
        spd += accel * dir;  spd -= hockeyistFrict * spd;

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

        pos += spd;  angSpd -= 0.0270190131 * angSpd;
        dir = sincos_fast(angle += turn + angSpd);
    }
};

struct PuckInfo : public UnitInfo
{
    Vec2D goalie[2];
    int goalieTime;

    PuckInfo() = default;

    PuckInfo(const Vec2D &pos_, const Vec2D &spd_)
    {
        pos = pos_;  spd = spd_;
        goalie[0] = Vec2D(rinkLeft  + hockeyistRad, pos_.y);
        goalie[1] = Vec2D(rinkRight - hockeyistRad, pos_.y);
    }

    int nextStep()
    {
        double y = max(goalCenter - goalieRange, min(goalCenter + goalieRange, pos.y));
        for(int i = 0; i < 2; i++)goalie[i].y = max(goalie[i].y - goalieSpd, min(goalie[i].y + goalieSpd, y));

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

double interceptProbability(const UnitInfo &info, int time);

struct PuckState : public UnitInfo
{
    double intercept;

    void set(const PuckInfo &info, double &intercept_, int time)
    {
        intercept_ = max(intercept_, interceptProbability(info, time));
        *static_cast<UnitInfo *>(this) = info;  intercept = intercept_;
    }
};



namespace MoveFlag
{
    enum
    {
        FORW = 0, FIRST_LEFT  = 0, SECOND_LEFT  = 0,
        BACK = 1, FIRST_RIGHT = 2, SECOND_RIGHT = 4,
    };
}

constexpr double stepValue(double val)
{
    return val < 0 ? 0 : (val > 1 ? 1 : val);
}

template<typename T, typename Info> struct Mapper
{
    void addMapPoint(Info &info, int time, int flags, double turnTime)
    {
        static_cast<T *>(this)->addMapPoint(info, time, flags, turnTime);
    }

    void closePath(Info &info, int flags, double turnTime)
    {
        static_cast<T *>(this)->closePath(info, flags, turnTime);
    }

    void fillMapTail(const Info &info, int flags, double accel, int time)
    {
        Info cur = info;
        for(int i = time + 1; i <= maxLookahead; i++)
        {
            cur.nextStep(accel, 0, i);  addMapPoint(cur, i, flags, time);
        }
        closePath(cur, flags, time);
    }

    void fillMapTail(const Info &info, int flags, double accel, double turn, int time)
    {
        double turnTime = time + maxTurnSteps;
        Info cur = info;  flags ^= MoveFlag::BACK;
        for(int i = time + 1; i <= maxLookahead; i++)
        {
            cur.nextStep(accel, turn * stepValue(turnTime - i + 1), i);
            addMapPoint(cur, i, flags, turnTime);
        }
        closePath(cur, flags, turnTime);
    }

    void fillMap(const Info& info, int flags, double accel1, double accel2, double turn)
    {
        constexpr int step = 4;

        Info cur = info;
        int n = min(maxLookahead, info.knockdown + int(maxTurnSteps));
        for(int i = info.knockdown + 1; i <= n; i++)
        {
            cur.nextStep(accel1, turn, i);  addMapPoint(cur, i, flags, i);  if(i % step)continue;  // TODO: globalTick ?
            fillMapTail(cur, flags, accel1, i);  fillMapTail(cur, flags, accel2, turn, i);
        }
    }

    void fillMap(const Info& info)
    {
        Info cur = info;  addMapPoint(cur, 0, 0, 0);
        for(int i = 1; i <= info.knockdown; i++)
        {
            cur.nextStep(0, 0, i);  addMapPoint(cur, i, 0, 0);
        }
        fillMapTail(cur, MoveFlag::BACK, accelMin, info.knockdown);
        fillMapTail(cur, MoveFlag::FORW, accelMax, info.knockdown);
        fillMap(cur, MoveFlag::BACK | MoveFlag::FIRST_LEFT,  accelMin, accelMax, +turnAngle);
        fillMap(cur, MoveFlag::BACK | MoveFlag::FIRST_RIGHT, accelMin, accelMax, -turnAngle);
        fillMap(cur, MoveFlag::FORW | MoveFlag::FIRST_LEFT,  accelMax, accelMin, +turnAngle);
        fillMap(cur, MoveFlag::FORW | MoveFlag::FIRST_RIGHT, accelMax, accelMin, -turnAngle);
    }
};

struct ReachabilityMap : public Mapper<ReachabilityMap, HockeyistInfo>
{
    vector<int> body, stick, filteredStick;

    void reset()
    {
        body.resize(gridSize);  stick.resize(gridSize);
        for(auto &flag :  body)flag = maxLookahead;
        for(auto &flag : stick)flag = maxLookahead;
    }

    void addMapPoint(const HockeyistInfo &info, int time, int flags, double turnTime)
    {
        auto &cell1 = body[gridPos(info.pos)];  cell1 = min<int>(cell1, time);
        auto &cell2 = stick[gridPos(info.pos + 0.5 * stickLength * info.dir)];
        cell2 = min<int>(cell2, max(time, info.cooldown));
    }

    void closePath(HockeyistInfo &info, int flags, double turnTime)
    {
    }

    void filter()
    {
        int x, y, cell = 0;
        filteredStick.resize(gridSize);
        for(y = 0; y < gridBorder; y++)
            for(x = 0; x < gridLine; x++, cell++)
                filteredStick[cell] = maxLookahead;
        for(; y < gridHeight - gridBorder; y++)
        {
            for(x = 0; x < gridBorder; x++, cell++)
                filteredStick[cell] = maxLookahead;
            for(; x < gridLine - gridBorder; x++, cell++)
            {
                int val = stick[cell];
                val = min<int>(val, stick[cell - gridLine - 1]);
                val = min<int>(val, stick[cell - gridLine + 1]);
                val = min<int>(val, stick[cell + gridLine - 1]);
                val = min<int>(val, stick[cell + gridLine + 1]);
                val = min<int>(val, stick[cell - gridLine]);
                val = min<int>(val, stick[cell + gridLine]);
                val = min<int>(val, stick[cell - 1]);
                val = min<int>(val, stick[cell + 1]);
                filteredStick[cell] = val;
            }
            for(; x < gridLine; x++, cell++)
                filteredStick[cell] = maxLookahead;
        }
        for(; y < gridHeight; y++)
            for(x = 0; x < gridLine; x++, cell++)
                filteredStick[cell] = maxLookahead;
        assert(cell == gridSize);
    }
};


struct PassTarget : public Vec2D
{
    int id, time, flags;
    double turnTime, score;

    PassTarget() = default;

    PassTarget(const Vec2D &pos, int time_) : Vec2D(pos), time(time_)
    {
    }

    bool operator < (const PassTarget &cmp) const
    {
        return score < cmp.score;
    }
};


ReachabilityMap enemyMap;
PuckState puckPath[maxLookahead + 1];
vector<PassTarget> targets;
int puckPathLen, goalFlag;


inline double gridWeight(double x, double y)
{
    constexpr double C[] = {16.0 / 27, 24.0 / 27, -15.0 / 27, 2.0 / 27};
    double r2 = x * x + y * y;  r2 = max(1.0, min(2.0, r2));
    return C[0] + r2 * (C[1] + r2 * (C[2] + r2 * C[3]));
}

double checkField(const vector<int> &grid, const Vec2D &pos, int time)
{
    Vec2D offs = (pos - rinkCenter) / cellSize;
    int x = lround(offs.x), y = lround(offs.y);  offs -= Vec2D(x, y);
    assert(abs(x) <= gridHalfWidth && abs(y) <= gridHalfHeight);

    double res = 0;
    int index = gridCenter + (y - gridBorder) * gridLine + x;
    for(int y = -gridBorder; y <= gridBorder; y++, index += gridLine)
        for(int x = -gridBorder; x <= gridBorder; x++)if(grid[index + x] <= time)
            res = max(res, gridWeight(offs.x + x, offs.y + y));  return res;
}

bool checkGoalie(const Vec2D &pos, double rad)
{
    double x = min(pos.y - rinkLeft - hockeyistRad, rinkRight + hockeyistRad - pos.y);
    double y = max(0.0, abs(pos.y - goalCenter) - goalHalf);
    return x * x + y * y < sqr(hockeyistRad + rad);
}

bool validPuckPos(const Vec2D &pos)
{
    if(pos.x < rinkLeft + puckRad || pos.x > rinkRight  - puckRad)return false;
    if(pos.y < rinkTop  + puckRad || pos.y > rinkBottom - puckRad)return false;
    return !checkGoalie(pos, puckRad);
}

double checkSafety(const HockeyistInfo &info, int time)
{
    if(checkGoalie(info.pos, hockeyistRad))return 0;
    double res = 0.5 * checkField(enemyMap.body, info.pos, time);
    if(!puckPathLen)
    {
        Vec2D puck = info.pos + hockeyistFrict * info.spd + holdDist * sincos(info.angle - info.angSpd);
        if(!validPuckPos(puck))return 0;

        //res = max(res, knockdownChance * checkField(enemyMap.stick, info.pos, time));
        //res = max(res, strikeChance * checkField(enemyMap.stick, puck, time));
        if(enemyMap.filteredStick[gridPos(info.pos)] <= time)res = max(res, knockdownChance);
        if(enemyMap.filteredStick[gridPos(puck)] <= time)res = max(res, strikeChance );
    }
    return 1 - res / (1 + sqr(time / 50.0));
}

struct AllyState : public HockeyistInfo
{
    double timeFactor, safety;

    AllyState(const HockeyistInfo &info) : HockeyistInfo(info), timeFactor(1), safety(1)
    {
    }

    void nextStep(double accel, double turn, int time)
    {
        HockeyistInfo::nextStep(accel, turn);  timeFactor *= timeGamma;
        safety = min(safety, checkSafety(*this, time));
    }
};


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

double interceptProbability(const UnitInfo &info, int time)
{
    double chance = min(maxChance, strikeChance - chanceDrop * info.spd.len());
    return chance * checkField(enemyMap.stick, info.pos, time) / (1 + sqr(time / 50.0));
}

double checkInterception(const Vec2D &pos, const Vec2D &spd, int time, int duration)
{
    PuckInfo info(pos, spd);  double res = 0;
    for(int i = 1; i <= duration; i++)
    {
        if(info.nextStep())break;  // TODO: check self goal
        res = max(res, interceptProbability(info, time + i));
    }
    return res;
}

inline Vec2D puckPos(const HockeyistInfo &info)
{
    return info.pos + hockeyistFrict * info.spd + holdDist * sincos(info.angle - info.angSpd);
}

inline bool inStrikeSector(const HockeyistInfo &info, const Vec2D &puck)
{
    Vec2D delta = puck - info.pos;
    if(!(delta.sqr() < sqr(stickLength)))return false;
    double dot = delta * info.dir, cross = delta % info.dir;
    return abs(cross) < dot * stickSectorTan;
}

struct StrikeInfo
{
    static constexpr double goalMultiplier = 2;
    static constexpr double dangerMultiplier = 10;

    int strikeTime, swingTime, targetIndex;
    double score, passAngle, passPower;

    void reset()
    {
        strikeTime = maxLookahead;  swingTime = maxSwing;  targetIndex = -1;
        score = -numeric_limits<double>::infinity();  passAngle = 0;  passPower = 1;
    }

    void tryStrikeFlyby(const AllyState &info, int time, double val)
    {
        double drop = chanceDrop * puckPath[time].spd.len();  int swing = -1;
        if(goalFlag)  // TODO: check strike dir
        {
            val *= min(maxChance, strikeChance - drop) * dangerMultiplier;
            if(pickChance - drop < maxChance)swing = 0;
        }
        else
        {
            val *= min(maxChance, pickChance - drop);
            int deltaTime = time - enemyMap.filteredStick[gridPos(puckPath[time].pos)];
            val *= (1 - puckPath[time].intercept) * (1 + deltaTime / maxLookahead);
        }
        if(!(val > score))return;

        strikeTime = time;  swingTime = swing;  targetIndex = -1;  score = val;
    }

    void evaluateStrike(const HockeyistInfo &info, const Vec2D &puck,
        double power, double beta, double sector, int time, double val)
    {
        Sector res = estimateGoalAngle(puck, info.spd, power, leftPlayer);
        if(!(res.halfSpan > 0))return;

        double angle = relAngle(res.centerAngle - info.angle), offs = max(0.0, abs(angle) - sector);
        val *= erf(beta * (offs + res.halfSpan)) - erf(beta * (offs - res.halfSpan));
        if(!((val *= 0.5 * goalMultiplier) > score))return;

        strikeTime = time;  swingTime = (sector > 0 ? -1 : maxSwing);  targetIndex = -1;
        score = val;  passAngle = angle;  passPower = 1;
    }

    void evaluatePass(const HockeyistInfo &info, const Vec2D &puck,
        double x, double y, int time, int catchTime, double val, int target)
    {
        Vec2D delta = Vec2D(x, y) - puck;
        double len2 = delta.sqr();  if(len2 < sqr(256))return;
        double angle = relAngle(atan2(delta.y, delta.x) - info.angle);
        if(!(abs(angle) < passSector))return;

        double len = sqrt(len2);  delta /= len;
        double passEndSpd = (pickChance - maxChance) / chanceDrop;  // TODO: to initConsts
        double relSpd = info.spd * delta, puckSpd = min(passBase + relSpd, passEndSpd + puckBeta * len);
        puckSpd = min(puckSpd, puckBeta * len / (1 - exp(-puckBeta * (catchTime - time))));
        double duration = 1 - puckBeta * len / puckSpd;  if(!(duration > 0))return;
        duration = -log(duration) / puckBeta;

        /*
        val *= 1 - checkInterception(pos, puckSpd * delta, time, lround(duration));
        if(!(val > score))return;
        */

        strikeTime = time;  swingTime = -1;  targetIndex = target;
        score = val;  passAngle = angle;  passPower = (puckSpd - relSpd) / passBase;
    }


    void tryStrike(const AllyState &info, int time)
    {
        double mul = info.safety * info.timeFactor;
        double val = mul * timeGammaSwing;  Vec2D puck;
        if(val * goalMultiplier > score)
        {
            HockeyistInfo cur = info;  for(int i = 0; i < maxSwing; i++)cur.nextStep(0, 0);
            do
            {
                if(!puckPathLen)puck = puckPos(cur);
                else if(!inStrikeSector(cur, puck = puckPath[time + maxSwing].pos))break;
                evaluateStrike(cur, puck, strikeBase + strikeGrowth * maxSwing, strikeBeta, 0, time, val);
            }
            while(false);
        }

        if(!puckPathLen)puck = puckPos(info);
        else if(!inStrikeSector(info, puck = puckPath[time].pos))return;
        if(mul * goalMultiplier > score)evaluateStrike(info, puck, passBase, passBeta, passSector, time, mul);
        if(time < puckPathLen)tryStrikeFlyby(info, time, mul);

        for(size_t i = 0; i < targets.size(); i++)
        {
            if(!((val = mul * targets[i].score) > score))continue;
            evaluatePass(info, puck, targets[i].x, targets[i].y, time, targets[i].time, val, i);
        }
    }

    bool operator < (const StrikeInfo &info) const
    {
        return score < info.score;
    }
};

struct MovePlan : public StrikeInfo
{
    static constexpr int strikeDelta = 8;

    int flags;
    double firstTurnEnd, secondTurnStart;

    struct Helper
    {
        double accBase, accDelta, turn1, turn2;
        double accTime1, accTime2, turnTime1, turnTime2;

        Helper(const MovePlan &plan)
        {
            accBase = accelMin;  accDelta = accelMax;
            if(plan.flags & MoveFlag::BACK)swap(accBase, accDelta);  accDelta -= accBase;
            turn1 = (plan.flags & MoveFlag::FIRST_RIGHT  ? -turnAngle : turnAngle);
            turn2 = (plan.flags & MoveFlag::SECOND_RIGHT ? -turnAngle : turnAngle);

            accTime2 = plan.secondTurnStart + maxTurnSteps;
            accTime1 = min(plan.firstTurnEnd - maxTurnSteps, accTime2);
            turnTime2 = plan.secondTurnStart;  turnTime1 = min(plan.firstTurnEnd, turnTime2);
        }

        double accel(int time) const
        {
            return accBase + accDelta * (stepValue(accTime2 - time) - stepValue(accTime1 - time));
        }

        double turn(int time) const
        {
            return turn1 * stepValue(turnTime1 - time) + turn2 * (1 - stepValue(turnTime2 - time));
        }
    };

    MovePlan() : flags(0), firstTurnEnd(0), secondTurnStart(maxLookahead)
    {
        reset();
    }

    MovePlan(const StrikeInfo &info, int fl, double turnTime) : StrikeInfo(info),
        flags(fl), firstTurnEnd(turnTime), secondTurnStart(info.strikeTime)
    {
    }

    void set(const PassTarget &target)
    {
        /*
        cout << "Target score: " << target.score << ", point: " << target.x << ':' << target.y;
        cout << ", flags: " << target.flags << ", time: " << target.turnTime << endl;
        */

        reset();  flags = target.flags;  firstTurnEnd = target.turnTime;
        secondTurnStart = strikeTime = target.time + 1;
    }

    void evaluate(const HockeyistInfo &info)
    {
        int minStrike = max(strikeTime - strikeDelta, info.cooldown);
        int maxStrike = min(strikeTime + strikeDelta, maxLookahead);
        reset();

        AllyState cur = info;  Helper helper(*this);
        for(int i = 0;; i++)
        {
            if(i >= minStrike)tryStrike(cur, i);
            if(i >= maxStrike)break;

            cur.nextStep(helper.accel(i), helper.turn(i), i + 1);
        }
    }

    MovePlan(const HockeyistInfo &info, const MovePlan &old, double step) :
        flags(old.flags), firstTurnEnd(old.firstTurnEnd), secondTurnStart(old.secondTurnStart)
    {
        strikeTime = old.strikeTime;

        firstTurnEnd += step * frand();
        if(firstTurnEnd < 0)
        {
            flags ^= MoveFlag::FIRST_RIGHT;  firstTurnEnd = -firstTurnEnd;
        }
        if(firstTurnEnd > 2 * maxTurnSteps)firstTurnEnd = 4 * maxTurnSteps - firstTurnEnd;

        secondTurnStart += 4 * step * frand();
        int minTurn2 = max(0, strikeTime - (int)ceil(maxTurnSteps));
        if(secondTurnStart < minTurn2)secondTurnStart = 2 * minTurn2 - secondTurnStart;
        if(secondTurnStart > strikeTime)
        {
            flags ^= MoveFlag::SECOND_RIGHT;  secondTurnStart = 2 * strikeTime - secondTurnStart;
        }

        evaluate(info);
    }

    void execute(Move &move)
    {
        cout << "Best score: " << score << ", flags: " << flags;
        cout << ", times: " << firstTurnEnd << ':' << secondTurnStart << ':' << strikeTime;
        if(swingTime >= 0)cout << ", swing: " << swingTime;
        else if(puckPathLen)cout << ", pickup";
        else cout << ", pass: " << passAngle << ", power: " << passPower;
        if(targetIndex >= 0)cout << ", target: " << targetIndex;
        cout << endl;

        if(!strikeTime)
        {
            move.setSpeedUp(0);  move.setTurn(0);
            if(swingTime > 0)
            {
                move.setAction(SWING);  swingTime--;  return;
            }
            else if(!swingTime)move.setAction(STRIKE);
            else if(puckPathLen)move.setAction(TAKE_PUCK);
            else
            {
                move.setPassPower(passPower);  move.setPassAngle(passAngle);  move.setAction(PASS);
            }
            *this = MovePlan();  return;
        }

        Helper helper(*this);  double accel = helper.accel(0);
        move.setSpeedUp(accel > 0 ? accel / accelMax : -accel / accelMin);
        move.setTurn(helper.turn(0));  move.setAction(NONE);

        strikeTime--;
        firstTurnEnd = max(0.0, firstTurnEnd - 1);
        secondTurnStart = max(0.0, secondTurnStart - 1);
    }
};

struct OptimizerState : public AllyState, public StrikeInfo
{
    OptimizerState(const HockeyistInfo &info) : AllyState(info)
    {
        reset();
    }
};

template<typename T> size_t chooseBest(vector<T> &array, size_t max)
{
    auto rend = array.rend();  make_heap(array.rbegin(), rend);  size_t n = 0;
    for(; n < max && rend != array.rbegin(); n++, --rend)pop_heap(array.rbegin(), rend);
    array.resize(n);  return n;
}

struct Optimizer : public Mapper<Optimizer, OptimizerState>
{
    vector<MovePlan> moves;

    void addMapPoint(OptimizerState &info, int time, int flags, double turnTime)
    {
        constexpr int step = 4;

        if(time < info.cooldown)return;
        if(time < MovePlan::strikeDelta || !(time % step))info.tryStrike(info, time);  // TODO: globalTick ?
    }

    void closePath(OptimizerState &info, int flags, double turnTime)
    {
        moves.emplace_back(info, flags, turnTime);
    }

    void optimizeMove(const HockeyistInfo &info, double step)
    {
        constexpr int survive = 4, offspring = 4;

        size_t n = chooseBest(moves, survive);
        for(size_t i = 0; i < n; i++)for(int j = 0; j < offspring; j++)
            moves.emplace_back(info, moves[i], step);
    }

    MovePlan findBestMove(const HockeyistInfo &info, const MovePlan &old)
    {
        moves.clear();  moves.push_back(old);  fillMap(info);

        double delta = 8;
        for(int i = 0; i < 8; i++, delta /= 3)optimizeMove(info, delta);

        const MovePlan *res = &old;
        double best = -numeric_limits<double>::infinity();
        for(auto &move : moves)if(move.score > best)
        {
            res = &move;  best = move.score;
        }
        return *res;
    }
};


struct AllyInfo : public HockeyistInfo, public Mapper<AllyInfo, AllyState>
{
    long long id;
    MovePlan plan;
    vector<pair<double, int>> stick;
    PassTarget defend, attack;
    int swinging;

    AllyInfo(const Hockeyist& hockeyist) : id(hockeyist.getId())
    {
        set(hockeyist);
    }

    void set(const Hockeyist& hockeyist)
    {
        HockeyistInfo::set(hockeyist);  swinging = hockeyist.getSwingTicks();
    }

    void addMapPoint(AllyState &info, int time, int flags, double turnTime)
    {
        Vec2D puck = info.pos + 0.5 * stickLength * info.dir;
        if(!validPuckPos(puck))return;  int cell = gridPos(puck);
        double score = info.safety * info.timeFactor;
        score *= 1 - sqr(2 * (info.pos.x - rinkCenter.x) / rinkWidth);
        score *= 1 - sqr(2 * (info.pos.y - rinkCenter.y) / rinkHeight);
        if(!(score > stick[cell].first))return;
        stick[cell].first = score;

        PassTarget target(puck, time);  target.id = id;
        target.flags = flags;  target.turnTime = turnTime;  target.score = score;
        if(cell == defendCell && target.score > defend.score)defend = target;

        int delta = enemyMap.filteredStick[cell] - max(time, info.cooldown);
        target.score *= (1 + double(delta) / maxLookahead);
        if((cell == attackCell[0] || cell == attackCell[1]) &&
            target.score > attack.score)attack = target;

        if(stick[cell].second < 0)
        {
            stick[cell].second = targets.size();  targets.push_back(target);
        }
        else targets[stick[cell].second] = target;
    }

    void closePath(AllyState &info, int flags, double turnTime)
    {
    }

    void fillMap()
    {
        stick.resize(gridSize);
        defend.flags = attack.flags = 0;
        defend.turnTime = attack.turnTime = 0;
        defend.score = -numeric_limits<double>::infinity();
        attack.score = -numeric_limits<double>::infinity();
        for(auto &flag : stick)flag = {0, -1};
        Mapper::fillMap(*this);
    }

    const MovePlan &chooseMove()
    {
        constexpr double oldBonus = 0.01;

        if(swinging)return plan;  plan.evaluate(*this);  plan.score += oldBonus;
        plan = Optimizer().findBestMove(*this, plan);  return plan;
    }

    void execute(Move &move)
    {
        if(swinging && plan.strikeTime)
        {
            move.setSpeedUp(0);  move.setTurn(0);  move.setAction(CANCEL_STRIKE);
        }
        else plan.execute(move);

        switch(move.getAction())
        {
        case TAKE_PUCK:      cout << ">>>>>>>>>>>>> TAKE_PUCK"     << endl;  break;
        case SWING:          cout << ">>>>>>>>>>>>> SWING"         << endl;  break;
        case STRIKE:         cout << ">>>>>>>>>>>>> STRIKE"        << endl;  break;
        case CANCEL_STRIKE:  cout << ">>>>>>>>>>>>> CANCEL_STRIKE" << endl;  break;
        case PASS:           cout << ">>>>>>>>>>>>> PASS"          << endl;  break;
        case SUBSTITUTE:     cout << ">>>>>>>>>>>>> SUBSTITUTE"    << endl;  break;
        default:  break;
        }
    }
};

struct EnemyInfo : public HockeyistInfo
{
    long long id;

    EnemyInfo(const Hockeyist& hockeyist) : id(hockeyist.getId())
    {
        set(hockeyist);
    }
};

template<typename T> void synchronize(vector<T> &list, vector<T *> &ref, const Hockeyist &hockeyist)
{
    for(auto iter = list.begin();; ++iter)
        if(iter == list.end())
        {
            list.emplace_back(hockeyist);  ref.push_back(&*list.rbegin());  break;
        }
        else if(iter->id == hockeyist.getId())
        {
            iter->set(hockeyist);  ref.push_back(&*iter);  break;
        }
}

vector<AllyInfo> allyInfo;
vector<AllyInfo *> allies;
vector<EnemyInfo> enemyInfo;
vector<EnemyInfo *> enemies;


uint8_t gradient(float val)
{
    static const uint8_t palette[26] =
        {0, 232, 233, 234, 235, 236, 237, 238, 239, 240, 241, 242, 243, 244, 245, 246, 247, 248, 249, 250, 251, 252, 253, 254, 255, 15};

    //return palette[min(max(0, int(val * 26)), 25)];
    return palette[1 + min(max(0, int(val * 25)), 24)];
}

inline void printSmall(uint8_t up, uint8_t dn)
{
    cout << "\x1B[48;5;" << unsigned(up) << ";38;5;" << unsigned(dn) << "m▄";
}

inline void printLarge(uint8_t cell)
{
    cout << "\x1B[48;5;" << unsigned(cell) << "m  ";
}

void drawField(const vector<int> &grid, int target)
{
    static bool first = true;

    if(first)
    {
        cout << "\x1B[2J";  first = false;
    }

    cout << "\x1B[0;0H";
    for(int cell = 0; cell + gridLine <= gridSize;)
    {
        for(int x = 0; x < gridLine; x++, cell++)
            printLarge(cell == target ? 10 : gradient(grid[cell] / double(maxLookahead)));
        cout << "\x1B[0m\n";
    }
    cout.flush();
}


void MyStrategy::move(const Hockeyist& self, const World& world, const Game& game, Move& move)
{
    const auto &player = world.getMyPlayer();
    if(player.isJustMissedGoal() || player.isJustScoredGoal())
    {
        move.setTurn(self.getTeammateIndex() & 1 ? -pi : pi);
        move.setSpeedUp(0);  move.setAction(NONE);  return;
    }

    if(globalTick != world.getTick())
    {
        if(!(globalTick = world.getTick()))
        {
            initConsts(game, world);  srand(game.getRandomSeed());
            allyInfo.reserve(6);  allies.reserve(3);
            enemyInfo.reserve(6);  enemies.reserve(3);
        }

        long long owner = world.getPuck().getOwnerHockeyistId();
        PuckInfo puck;  puck.set(world.getPuck());  int totalGoals = 0;
        for(auto &player : world.getPlayers())totalGoals += player.getGoalCount();
        puck.goalieTime = (totalGoals ? 999999 : world.getTickCount());

        allies.clear();  enemies.clear();
        for(auto &hockeyist : world.getHockeyists())
        {
            if(hockeyist.getType() == GOALIE)
            {
                puck.goalie[hockeyist.isTeammate() ? 0 : 1] = Vec2D(hockeyist.getX(), hockeyist.getY());  continue;
            }
            if(hockeyist.isTeammate())synchronize(allyInfo, allies, hockeyist);
            else synchronize(enemyInfo, enemies, hockeyist);
        }

        enemyMap.reset();
        for(auto enemy : enemies)enemyMap.fillMap(*enemy);
        enemyMap.filter();

        targets.clear();  AllyInfo *allyPuck = nullptr;
        for(auto ally : allies)
            if(ally->id == owner)allyPuck = &*ally;
            else ally->fillMap();
        chooseBest(targets, 8);

        puckPathLen = 0;
        if(owner < 0)
        {
            static_cast<UnitInfo &>(puckPath[0]) = puck;
            double intercept = puckPath[0].intercept = 0;  goalFlag = 0;
            for(puckPathLen = 1; puckPathLen <= maxLookahead; puckPathLen++)
            {
                if((goalFlag = puck.nextStep()))break;
                puckPath[puckPathLen].set(puck, intercept, puckPathLen);
            }
            if(!leftPlayer)goalFlag = -goalFlag;
            if(goalFlag <= 0)
            {
                double best = -numeric_limits<double>::infinity();
                for(auto ally : allies)
                {
                    double score = ally->chooseMove().score;  if(!(score > best))continue;
                    best = score;  owner = ally->id;
                }
            }
        }

        int tg = -1;
        if(allyPuck)tg = allyPuck->chooseMove().targetIndex;
        long long pass = (tg < 0 ? -1 : targets[tg].id);
        for(auto ally : allies)
        {
            if(ally->id == pass)ally->plan.set(targets[tg]);
            else if(ally->id != owner)ally->plan.set(ally->defend.score > 0 ? ally->defend : ally->attack);
        }

        /*
        if(tg >= 0)
        {
            cout << "Target: " << targets[tg].x << ':' << targets[tg].y;
            cout << ", times:" << targets[tg].time << ':' << allyPuck->plan.strikeTime << endl;
        }
        */

        //drawField(enemyMap.body, tg < 0 ? -1 : gridPos(targets[tg]));
    }

    for(auto ally : allies)if(ally->id == self.getId())
    {
        ally->execute(move);  return;
    }
}

MyStrategy::MyStrategy()
{
    cout << setprecision(16);
}
