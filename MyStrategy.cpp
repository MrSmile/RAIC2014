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

inline constexpr Vec2D rotate(const Vec2D &v1, const Vec2D &v2)
{
    return Vec2D(v1.x * v2.x - v1.y * v2.y, v1.x * v2.y + v1.y * v2.x);
}

inline constexpr Vec2D conj(const Vec2D &v)
{
    return Vec2D(v.x, -v.y);
}




constexpr int maxLookahead = 256;
constexpr int turnStep = 4, strikeStep = 4, strikeDelta = 8;
constexpr int optSurvive = 4, optOffspring = 4, optStepCount = 4, optFinalCount = 4;
constexpr double optStepBase = 8, optStepMul = 0.5;

constexpr double timeGamma = 1 - 1.0 / maxLookahead;
constexpr double cellSize = 32;

constexpr double hockeyistFrict = 0.02, puckFrict = 0.001, puckBeta = -log(1 - puckFrict);
constexpr double wallBounce = 0.25, hockeyistBounce = 0.325, goalieFrict = 0.1;
constexpr double minDepth = 0.01, depthFactor = 0.8;

double rinkLeft, rinkRight, rinkTop, rinkBottom;
double rinkWidth, rinkHeight;  Vec2D rinkCenter;
int gridHalfWidth, gridHalfHeight, gridLine, gridHeight, gridCenter, gridSize;
double goalCenter, goalHalf, hockeyistRad, puckRad, goalieSpd, goalieRange;
double accelMin, accelMax, turnAngle, maxTurnSteps;
int maxTurnCount;  Vec2D turnRot, endTurnRot;
double stickLength, stickSectorTan, holdDist, timeGammaSwing;  int maxSwing;
double strikeBase, strikeGrowth, passBase, strikeBeta, passBeta;  Vec2D passSector;
double minChance, maxChance, pickChance, strikeChance, chanceDrop;
double knockdownChance, takeAwayChance;

int globalTick = -1, goalieTime;
int defendCell, attackCell[2];
bool leftPlayer;


void initConsts(const Game& game, const World& world)
{
    rinkLeft = game.getRinkLeft();  rinkRight = game.getRinkRight();
    rinkTop = game.getRinkTop();  rinkBottom = game.getRinkBottom();
    rinkWidth = rinkRight - rinkLeft;  rinkHeight = rinkBottom - rinkTop;
    rinkCenter.x = (rinkLeft + rinkRight) / 2;
    rinkCenter.y = (rinkTop + rinkBottom) / 2;

    gridHalfWidth  = lround(rinkWidth  / (2 * cellSize)) + 3;
    gridHalfHeight = lround(rinkHeight / (2 * cellSize)) + 3;
    gridLine = 2 * gridHalfWidth + 1;  gridHeight = 2 * gridHalfHeight + 1;
    gridCenter = gridHalfHeight * gridLine + gridHalfWidth;
    gridSize = gridHeight * gridLine;

    goalHalf = game.getGoalNetHeight() / 2;  goalCenter = game.getGoalNetTop() + goalHalf;
    hockeyistRad = world.getHockeyists()[0].getRadius();  puckRad = world.getPuck().getRadius();
    goalieSpd = game.getGoalieMaxSpeed();  goalieRange = goalHalf - hockeyistRad;
    accelMin = -game.getHockeyistSpeedDownFactor();  accelMax = game.getHockeyistSpeedUpFactor();
    turnAngle = game.getHockeyistTurnAngleFactor();  maxTurnSteps = (pi / 2) / turnAngle;  maxTurnCount = int(maxTurnSteps);
    turnRot = sincos(turnAngle);  endTurnRot = sincos(turnAngle * (maxTurnSteps - maxTurnCount));

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
    passSector = sincos(game.getPassSector() / 2);

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
    Vec2D offs = (pos - rinkCenter) / cellSize + Vec2D(gridHalfWidth + 0.5, gridHalfHeight + 0.5);
    return int(offs.y) * gridLine + int(offs.x);
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
    Vec2D dir;  Vec2D *rot;
    int knockdown, cooldown;

    void set(const Hockeyist& hockeyist, Vec2D *rot_)
    {
        double angSpd = hockeyist.getAngularSpeed();
        rot = rot_;  rot[0] = sincos(angSpd);
        for(int i = 1; i <= maxLookahead; i++)
        {
            angSpd -= 0.0270190131 * angSpd;  rot[i] = sincos(angSpd);
        }

        UnitInfo::set(hockeyist);
        dir = sincos(hockeyist.getAngle());
        knockdown = hockeyist.getRemainingKnockdownTicks();
        cooldown = max(knockdown, hockeyist.getRemainingCooldownTicks());
    }

    void nextStep(double accel, const Vec2D &turn, int time)
    {
        spd += accel * dir;  spd -= hockeyistFrict * spd;

        double delta;
        if((delta = rinkLeft + hockeyistRad - pos.x) > 0)
        {
            double sign = pos.y - goalCenter, offs = abs(sign) - goalHalf;
            if(offs < 0)
            {
                if(delta > +1)pos.x = rinkLeft - hockeyistRad - 1;
            }
            else if(!bounceCircle(Vec2D(rinkLeft, goalCenter + copysign(goalHalf, sign)), hockeyistRad, hockeyistBounce) &&
                !bounceCircle(Vec2D(rinkLeft, goalCenter + copysign(goalHalf + 25, sign)), hockeyistRad, hockeyistBounce / 10))
            {
                if(offs > 25)
                {
                    if(spd.x < 0)spd.x *= -hockeyistBounce;
                    if(spd.x < delta)pos.x += depthFactor * (delta + minDepth);
                }
                else
                {
                    if(spd.x < 0)spd.x *= -hockeyistBounce / 10;
                    if(spd.x < delta)pos.x += depthFactor * (delta + minDepth);
                }
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

        pos += spd;  dir = rotate(dir, rot[time]);  dir = rotate(dir, turn);
    }
};

struct PuckInfo : public UnitInfo
{
    Vec2D goalie[2];

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
        constexpr double eps = 10;
        if((delta = rinkLeft + puckRad - pos.x) > 0)
        {
            if(abs(pos.y - goalCenter) < goalHalf + eps)
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
            if(abs(pos.y - goalCenter) < goalHalf + eps)
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

inline Vec2D puckPos(const HockeyistInfo &info, int time)
{
    return info.pos + hockeyistFrict * info.spd + holdDist * rotate(info.dir, conj(info.rot[time]));
}

struct PuckState
{
    Vec2D pos;  double intercept, pick, strike;

    void set(const PuckInfo &info, int time, double intercept_ = 0)
    {
        pos = info.pos;  intercept = intercept_;
        double drop = chanceDrop * info.spd.len();
        pick = min(maxChance, pickChance - drop);
        strike = min(maxChance, strikeChance - drop);
    }

    void set(const HockeyistInfo &info, int time)
    {
        pos = puckPos(info, time);  intercept = 0;
        pick = 0.25;  strike = 0.75;  // TODO: initConsts
    }
};



namespace MoveFlag
{
    enum
    {
        FORW = 0, FIRST_LEFT  = 0, SECOND_LEFT  = 0,
        BACK = 1, FIRST_RIGHT = 2, SECOND_RIGHT = 4,
        STOP = 8
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
            cur.nextStep(accel, Vec2D(1, 0), i);  addMapPoint(cur, i, flags, time);
        }
        closePath(cur, flags, time);
    }

    void fillMapTail(const Info &info, int flags, double accel, const Vec2D &turn, const Vec2D &endTurn, int time)
    {
        double turnTime = time + maxTurnSteps;
        Info cur = info;  flags ^= MoveFlag::BACK;
        int i = time, end = time + maxTurnCount;
        for(i++; i <= maxLookahead; i++)
        {
            if(i > end)
            {
                cur.nextStep(accel, endTurn, i);  addMapPoint(cur, i, flags, turnTime);  break;
            }
            cur.nextStep(accel, turn, i);  addMapPoint(cur, i, flags, i);
        }
        for(i++; i <= maxLookahead; i++)
        {
            cur.nextStep(accel, Vec2D(1, 0), i);  addMapPoint(cur, i, flags, turnTime);
        }
        closePath(cur, flags, turnTime);
    }

    void fillMap(const Info& info, int flags)
    {
        double accel1 = accelMax, accel2 = accelMin;
        if(flags & MoveFlag::BACK)swap(accel1, accel2);
        Vec2D turn = turnRot, endTurn = endTurnRot;
        if(flags & MoveFlag::FIRST_RIGHT)
        {
            turn.y = -turn.y;  endTurn.y = -endTurn.y;
        }

        Info cur = info;
        int n = min(maxLookahead, info.knockdown + maxTurnCount);
        for(int i = info.knockdown + 1; i <= n; i++)
        {
            cur.nextStep(accel1, turn, i);  addMapPoint(cur, i, flags, i);  if(i % turnStep)continue;  // TODO: globalTick ?
            fillMapTail(cur, flags, accel1, i);  fillMapTail(cur, flags, accel2, turn, endTurn, i);
        }
    }

    void fillMap(const Info& info)
    {
        Info cur = info;  addMapPoint(cur, 0, 0, 0);
        for(int i = 1; i <= info.knockdown; i++)
        {
            cur.nextStep(0, Vec2D(1, 0), i);  addMapPoint(cur, i, 0, 0);
        }
        fillMapTail(cur, MoveFlag::BACK, accelMin, info.knockdown);
        fillMapTail(cur, MoveFlag::FORW, accelMax, info.knockdown);
        fillMap(cur, MoveFlag::BACK | MoveFlag::FIRST_LEFT);
        fillMap(cur, MoveFlag::BACK | MoveFlag::FIRST_RIGHT);
        fillMap(cur, MoveFlag::FORW | MoveFlag::FIRST_LEFT);
        fillMap(cur, MoveFlag::FORW | MoveFlag::FIRST_RIGHT);
    }
};

struct ReachabilityMap : public Mapper<ReachabilityMap, HockeyistInfo>
{
    vector<int> body, stick;

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

    static void filterField(vector<int> &field)
    {
        constexpr int border = 1;

        vector<int> old(gridSize);  swap(old, field);

        int x, y, cell = 0;
        for(y = 0; y < border; y++)
            for(x = 0; x < gridLine; x++, cell++)
                field[cell] = maxLookahead;
        for(; y < gridHeight - border; y++)
        {
            for(x = 0; x < border; x++, cell++)
                field[cell] = maxLookahead;
            for(; x < gridLine - border; x++, cell++)
            {
                int val = old[cell];
                val = min<int>(val, old[cell - gridLine - 1]);
                val = min<int>(val, old[cell - gridLine + 1]);
                val = min<int>(val, old[cell + gridLine - 1]);
                val = min<int>(val, old[cell + gridLine + 1]);
                val = min<int>(val, old[cell - gridLine]);
                val = min<int>(val, old[cell + gridLine]);
                val = min<int>(val, old[cell - 1]);
                val = min<int>(val, old[cell + 1]);
                field[cell] = val;
            }
            for(; x < gridLine; x++, cell++)
                field[cell] = maxLookahead;
        }
        for(; y < gridHeight; y++)
            for(x = 0; x < gridLine; x++, cell++)
                field[cell] = maxLookahead;
        assert(cell == gridSize);
    }

    void filter()
    {
        filterField(body);  filterField(stick);
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
HockeyistInfo *enemyPuck;


inline double gridWeight(double x, double y)
{
    constexpr double C[] = {16.0 / 27, 24.0 / 27, -15.0 / 27, 2.0 / 27};
    double r2 = x * x + y * y;  r2 = max(1.0, min(2.0, r2));
    return C[0] + r2 * (C[1] + r2 * (C[2] + r2 * C[3]));
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

    double res = 0;
    if(!puckPathLen)
    {
        Vec2D puck = puckPos(info, time);  if(!validPuckPos(puck))return 0;
        if(enemyMap.stick[gridPos(info.pos)] <= time)res = max(res, knockdownChance);
        if(enemyMap.stick[gridPos(puck)] <= time)res = max(res, strikeChance );
    }
    return 1 - res;
}

template<bool ally> struct StateScore : public HockeyistInfo
{
    double timeFactor;

    StateScore(const HockeyistInfo &info) : HockeyistInfo(info), timeFactor(1)
    {
    }

    void nextStep(double accel, const Vec2D &turn, int time)
    {
        HockeyistInfo::nextStep(accel, turn, time);  timeFactor *= timeGamma;
    }

    double multiplier() const
    {
        return timeFactor;
    }
};

template<> struct StateScore<true> : public HockeyistInfo
{
    double timeFactor, safety;

    StateScore(const HockeyistInfo &info) : HockeyistInfo(info), timeFactor(1), safety(1)
    {
    }

    void nextStep(double accel, const Vec2D &turn, int time)
    {
        HockeyistInfo::nextStep(accel, turn, time);  timeFactor *= timeGamma;
        safety = min(safety, checkSafety(*this, time));
    }

    double multiplier() const
    {
        return safety * timeFactor;
    }
};


struct Sector
{
    Vec2D dir;
    double span, time;

    Sector(const Vec2D &dir_ = {0, 0}, double span_ = 0, double time_ = 0) :
        dir(dir_), span(span_), time(time_)
    {
    }
};

struct GoalHelper
{
    Vec2D dir, bonus;
    double rad, power, dist, offs, end;

    GoalHelper(const Vec2D &pos, const Vec2D &spd, double power_) :
        dir(pos), bonus(spd), rad(hockeyistRad + puckRad), power(power_)
    {
    }

    void flipHorz()
    {
        dir.x = 2 * rinkCenter.x - dir.x;  bonus.x = -bonus.x;
    }

    void flipVert()
    {
        dir.y = 2 * goalCenter - dir.y;  bonus.y = -bonus.y;
    }

    bool init()
    {
        dir -= Vec2D(rinkLeft, goalCenter - goalHalf);
        dist = dir.x - hockeyistRad;  if(!(dist > rad))return false;
        offs = max(0.0, dir.y - 2 * goalHalf);  dir = normalize(dir);
        return true;
    }

    bool iterate(bool first)
    {
        double invSpd = 1 / (power - bonus * dir), cmp = 0;
        end = (dist - rad * dir.y) * invSpd;
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
};

Sector estimateGoalAngle(const Vec2D &pos, const Vec2D &spd, double power, bool right)
{
    GoalHelper helper(pos, spd, power);
    if(right)helper.flipHorz();  if(pos.y < goalCenter)helper.flipVert();

    if(!helper.init())return Sector();  Vec2D dir = helper.dir;
    if(!helper.iterate(true))return Sector();  double time = helper.end;
    for(int i = 0; i < 3; i++)helper.iterate(false);

    double dot = helper.dir * dir, cross = helper.dir % dir;
    double norm = 1 / sqrt(2 + 2 * dot);  dir += helper.dir;
    if(!right)dir.x = -dir.x;  if(!(pos.y < goalCenter))dir.y = -dir.y;
    return Sector(dir * norm, cross * norm, time);
}

double interceptProbability(const UnitInfo &info, int time)
{
    if(enemyMap.stick[gridPos(info.pos)] > time)return 0;
    return min(maxChance, strikeChance - chanceDrop * info.spd.len());
}

double checkInterception(const Vec2D &pos, const Vec2D &spd, int time, int duration)
{
    UnitInfo info;  double res = 0;
    info.pos = pos;  info.spd = spd;
    for(int i = 1; i <= duration; i++)
    {
        info.spd -= puckFrict * info.spd;  info.pos += info.spd;
        if(abs(info.pos.x - rinkCenter.x) > rinkWidth  / 2 + cellSize)break;

        double delta;
        if((delta = rinkTop + puckRad - info.pos.y) > 0)
        {
            if(info.spd.y < 0)info.spd.y *= -wallBounce;
            if(info.spd.y < delta)info.pos.y += depthFactor * (delta - minDepth);
        }
        if((delta = rinkBottom - puckRad - info.pos.y) < 0)
        {
            if(info.spd.y > 0)info.spd.y *= -wallBounce;
            if(info.spd.y > delta)info.pos.y += depthFactor * (delta + minDepth);
        }
        res = max(res, interceptProbability(info, time + i));
    }
    return res;
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
    static constexpr double goalMultiplier = 10;
    static constexpr double dangerMultiplier = 10;

    int strikeTime, swingTime, targetIndex;
    double score, passPower;  Vec2D passDir;

    void reset()
    {
        strikeTime = maxLookahead;  swingTime = -1;  targetIndex = -1;
        score = 0;  passDir = {0, 0};  passPower = 1;
    }

    void tryStrikeFlyby(int time, double val)
    {
        int swing = -1;
        if(goalFlag)  // TODO: check strike dir
        {
            val *= puckPath[time].strike * dangerMultiplier;
            if(puckPath[time].pick < maxChance)swing = 0;
        }
        else
        {
            val *= puckPath[time].pick;
            int deltaTime = time - enemyMap.stick[gridPos(puckPath[time].pos)];
            val *= (1 - puckPath[time].intercept) * (1 + deltaTime / maxLookahead);
        }
        if(!(val > score))return;

        strikeTime = time;  swingTime = swing;  targetIndex = -1;  score = val;
    }

    template<bool ally> void evaluateStrike(const HockeyistInfo &info, const Vec2D &puck,
        double power, double beta, bool sector, int time, double val)
    {
        Sector res = estimateGoalAngle(puck, info.spd, power, ally == leftPlayer);
        if(!(res.span > 0))return;

        Vec2D delta = rotate(res.dir, conj(info.dir)), offs(delta.x, abs(delta.y));
        if(sector)offs = (offs.x > passSector.x ? Vec2D(1, 0) : rotate(offs, conj(passSector)));
        if(!(offs.x > 0))return;

        val *= erf(beta * (offs.y + res.span)) - erf(beta * (offs.y - res.span));
        if(!((val *= 0.5 * goalMultiplier) > score))return;

        if(ally)
        {
            Vec2D dir = res.dir;  // TODO: better approx
            val *= 1 - checkInterception(puck, (power + info.spd * dir) * dir, time, lround(res.time));
            if(!(val > score))return;
        }

        strikeTime = time;  swingTime = (sector ? -1 : maxSwing);  targetIndex = -1;
        score = val;  passPower = 1;  passDir = delta;
    }

    void evaluatePass(const HockeyistInfo &info, const Vec2D &puck,
        double x, double y, int time, int catchTime, double val, int target)
    {
        constexpr double minPassDist = 256;

        Vec2D delta = Vec2D(x, y) - puck;
        double len2 = delta.sqr();  if(len2 < sqr(minPassDist))return;
        double len = sqrt(len2);  delta /= len;
        Vec2D offs = rotate(delta, conj(info.dir));
        if(!(offs.x > passSector.x))return;

        double passEndSpd = (pickChance - maxChance) / chanceDrop;  // TODO: to initConsts
        double relSpd = info.spd * delta, puckSpd = min(passBase + relSpd, passEndSpd + puckBeta * len);
        puckSpd = min(puckSpd, puckBeta * len / (1 - exp(-puckBeta * (catchTime - time))));
        double duration = 1 - puckBeta * len / puckSpd;  if(!(duration > 0))return;
        duration = -log(duration) / puckBeta;

        if(min(x, puck.x) < rinkLeft + (2 * hockeyistRad + puckRad))
        {
            double cross = delta % (Vec2D(rinkLeft + hockeyistRad, goalCenter) - puck);
            double cmp = abs(delta.x) * goalieRange + (hockeyistRad + puckRad);
            if(cross > -cmp && cross < cmp)return;
        }
        if(max(x, puck.x) > rinkRight - (2 * hockeyistRad + puckRad))
        {
            double cross = delta % (Vec2D(rinkRight - hockeyistRad, goalCenter) - puck);
            double cmp = abs(delta.x) * goalieRange + (hockeyistRad + puckRad);
            if(cross > -cmp && cross < cmp)return;
        }
        val *= 1 - checkInterception(puck, puckSpd * delta, time, lround(duration));
        if(!(val > score))return;

        strikeTime = time;  swingTime = -1;  targetIndex = target;
        score = val;  passPower = (puckSpd - relSpd) / passBase;  passDir = offs;
    }


    template<bool ally> void tryStrike(const StateScore<ally> &info, int time)
    {
        double mul = info.multiplier();  Vec2D puck;
        if(time + maxSwing <= maxLookahead && mul * timeGammaSwing * goalMultiplier > score)do
        {
            StateScore<ally> cur = info;
            for(int i = 1; i <= maxSwing; i++)cur.nextStep(0, Vec2D(1, 0), time + i);

            double val = cur.multiplier();
            if(!puckPathLen)puck = puckPos(cur, time + maxSwing);
            else
            {
                if(!inStrikeSector(cur, puck = puckPath[time + maxSwing].pos))break;
                val *= puckPath[time + maxSwing].strike * (1 - puckPath[time + maxSwing].intercept);
            }
            evaluateStrike<ally>(cur, puck, strikeBase + strikeGrowth * maxSwing, strikeBeta, false, time, val);
        }
        while(false);

        if(!puckPathLen)puck = puckPos(info, time);
        else
        {
            if(!inStrikeSector(info, puck = puckPath[time].pos))return;
            mul *= puckPath[time].strike * (1 - puckPath[time].intercept);
        }
        if(mul * goalMultiplier > score)evaluateStrike<ally>(info, puck, passBase, passBeta, true, time, mul);
        if(!ally)return;

        if(puckPathLen)
        {
            if(time < puckPathLen)tryStrikeFlyby(time, mul);  return;
        }

        double y1 = rinkTop + puckRad, y2 = rinkBottom + puckRad;
        for(size_t i = 0; i < targets.size(); i++)
        {
            double val = mul * targets[i].score;  if(!(val > score))continue;
            evaluatePass(info, puck, targets[i].x, targets[i].y, time, targets[i].time, val, i);
            evaluatePass(info, puck, targets[i].x, y1 - (targets[i].y - y1) / wallBounce, time, targets[i].time, val, i);
            evaluatePass(info, puck, targets[i].x, y2 - (targets[i].y - y2) / wallBounce, time, targets[i].time, val, i);
        }
    }

    bool operator < (const StrikeInfo &info) const
    {
        return score < info.score;
    }
};

struct MovePlan : public StrikeInfo
{
    int flags;
    double firstTurnEnd, secondTurnStart;

    struct Helper
    {
        double accBase, accDelta, accTime1, accTime2;
        int turnStep[3], pos;  Vec2D rot[5];

        Helper(const MovePlan &plan) : pos(0)
        {
            accBase = accelMin;  accDelta = accelMax;
            if(plan.flags & MoveFlag::BACK)swap(accBase, accDelta);  accDelta -= accBase;
            accTime2 = plan.secondTurnStart + maxTurnSteps;
            accTime1 = min(plan.firstTurnEnd - maxTurnSteps, accTime2);

            double turnTime2 = plan.secondTurnStart;
            double turnTime1 = min(plan.firstTurnEnd, turnTime2);
            turnStep[0] = (int)floor(turnTime1);  turnTime1 -= turnStep[0];
            turnStep[1] = (int)floor(turnTime2);  turnTime2 -= turnStep[1];
            turnStep[2] = maxLookahead;

            rot[0] = turnRot;  rot[4] = turnRot;
            double angle1 = turnAngle * turnTime1;
            double angle2 = turnAngle * turnTime2;
            if(plan.flags & MoveFlag::FIRST_RIGHT)
            {
                rot[0].y = -rot[0].y;  angle1 = -angle1;
            }
            if(plan.flags & MoveFlag::SECOND_RIGHT)
            {
                rot[4].y = -rot[4].y;  angle2 = -angle2;
            }
            if(turnStep[0] == turnStep[1])
            {
                turnStep[1] = maxLookahead;  rot[1] = sincos(angle1 + angle2);  rot[2] = rot[4];
            }
            else
            {
                rot[1] = sincos(angle1);  rot[2] = Vec2D(1, 0);  rot[3] = sincos(angle2);
            }
        }

        double accel(int time) const
        {
            return accBase + accDelta * (stepValue(accTime2 - time) - stepValue(accTime1 - time));
        }

        Vec2D turn(int time)
        {
            return time < turnStep[pos] ? rot[2 * pos] : rot[2 * (++pos) - 1];
        }
    };

    MovePlan() : flags(MoveFlag::STOP), firstTurnEnd(0), secondTurnStart(maxLookahead)
    {
        reset();
    }

    MovePlan(const StrikeInfo &info, int flags_, double turnTime) : StrikeInfo(info),
        flags(flags_), firstTurnEnd(turnTime), secondTurnStart(info.strikeTime)
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

    template<bool ally> void evaluate(const HockeyistInfo &info)
    {
        int minStrike = max(strikeTime - strikeDelta, info.cooldown);
        int maxStrike = min(strikeTime + strikeDelta, maxLookahead);
        reset();

        StateScore<ally> cur = info;  Helper helper(*this);
        for(int i = 0;; i++)
        {
            if(i >= minStrike)tryStrike<ally>(cur, i);  if(i >= maxStrike)break;
            cur.nextStep(helper.accel(i), helper.turn(i), i + 1);
        }
    }

    MovePlan(const HockeyistInfo &info, const MovePlan &old, double step, bool ally) :
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
            flags ^= MoveFlag::SECOND_RIGHT;
            secondTurnStart = 2 * strikeTime - secondTurnStart;
        }

        if(ally)evaluate<true>(info);  else evaluate<false>(info);
    }

    void skipTurn()
    {
        if(flags & MoveFlag::STOP)return;
        if(!strikeTime)
        {
            if(swingTime > 0)swingTime--;
            else flags = MoveFlag::STOP;  return;
        }
        firstTurnEnd = max(0.0, firstTurnEnd - 1);
        secondTurnStart = max(0.0, secondTurnStart - 1);
        strikeTime--;
    }

    void execute(Move &move)
    {
        /*
        cout << "Best score: " << score << ", flags: " << flags;
        if(!(flags & MoveFlag::STOP))
        {
            cout << ", times: " << firstTurnEnd << ':' << secondTurnStart << ':' << strikeTime;
            if(swingTime >= 0)cout << ", swing: " << swingTime;
            else if(puckPathLen)cout << ", pickup";
            else cout << ", pass: " << passDir.x << ':' << passDir.y <<
                ", power: " << passPower;
            if(targetIndex >= 0)cout << ", target: " << targetIndex;
        }
        cout << endl;
        */

        if(flags & MoveFlag::STOP)
        {
            move.setSpeedUp(0);  move.setTurn(0);  move.setAction(NONE);  return;
        }

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
                move.setPassPower(passPower);
                move.setPassAngle(atan2(passDir.y, passDir.x));  move.setAction(PASS);
            }
            flags = MoveFlag::STOP;  return;
        }

        Helper helper(*this);
        double accel = helper.accel(0);  Vec2D rot = helper.turn(0);
        move.setSpeedUp(accel > 0 ? accel / accelMax : -accel / accelMin);
        move.setTurn(atan2(rot.y, rot.x));  move.setAction(NONE);

        firstTurnEnd = max(0.0, firstTurnEnd - 1);
        secondTurnStart = max(0.0, secondTurnStart - 1);
        strikeTime--;
    }
};

template<bool ally> struct OptimizerState : public StateScore<ally>, public StrikeInfo
{
    OptimizerState(const HockeyistInfo &info) : StateScore<ally>(info)
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

template<bool ally> struct Optimizer : public Mapper<Optimizer<ally>, OptimizerState<ally>>
{
    vector<MovePlan> moves;

    void addMapPoint(OptimizerState<ally> &info, int time, int flags, double turnTime)
    {
        if(time < info.cooldown)return;
        if(time < strikeDelta || !(time % strikeStep))info.StrikeInfo::tryStrike<ally>(info, time);  // TODO: globalTick ?
    }

    void closePath(OptimizerState<ally> &info, int flags, double turnTime)
    {
        moves.emplace_back(info, flags, turnTime);
    }

    void optimizeMove(const HockeyistInfo &info, double step, int survive, int offspring)
    {
        size_t n = chooseBest(moves, survive);
        for(size_t i = 0; i < n; i++)for(int j = 0; j < offspring; j++)
            moves.emplace_back(info, moves[i], step, ally);
    }

    MovePlan findBestMove(const HockeyistInfo &info, const MovePlan &old)
    {
        moves.clear();  moves.push_back(old);
        Mapper<Optimizer<ally>, OptimizerState<ally>>::fillMap(info);

        double delta = optStepBase;
        for(int i = 0; i < optStepCount; i++, delta *= optStepMul)
            optimizeMove(info, delta, optSurvive, optOffspring);
        if(ally)for(int i = 0; i < optFinalCount; i++, delta *= optStepMul)
            optimizeMove(info, delta, 1, optOffspring);

        const MovePlan *res = &old;
        double best = -numeric_limits<double>::infinity();
        for(auto &move : moves)if(move.score > best)
        {
            res = &move;  best = move.score;
        }
        return *res;
    }
};


struct AllyInfo : public HockeyistInfo, public Mapper<AllyInfo, StateScore<true>>
{
    long long id;
    vector<Vec2D> rot;
    MovePlan plan;  bool activePlan;
    vector<pair<double, int>> stick;
    PassTarget defend, attack;
    Vec2D defendPoint;
    int swinging;

    AllyInfo(const Hockeyist& hockeyist) : id(hockeyist.getId()), rot(maxLookahead + 1)
    {
        set(hockeyist);  double offs = rinkWidth / 2 - 2 * goalHalf;
        defendPoint = rinkCenter + Vec2D(leftPlayer ? -offs : offs, 0);
    }

    void set(const Hockeyist& hockeyist)
    {
        HockeyistInfo::set(hockeyist, rot.data());  swinging = hockeyist.getSwingTicks();
    }

    void addMapPoint(StateScore<true> &info, int time, int flags, double turnTime)
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

        int delta = enemyMap.stick[cell] - max(time, info.cooldown);
        target.score *= (1 + double(delta) / maxLookahead);
        if((cell == attackCell[0] || cell == attackCell[1]) &&
            target.score > attack.score)attack = target;

        if(stick[cell].second < 0)
        {
            stick[cell].second = targets.size();  targets.push_back(target);
        }
        else targets[stick[cell].second] = target;
    }

    void closePath(StateScore<true> &info, int flags, double turnTime)
    {
    }

    bool tryKnockdown()
    {
        if(!enemyPuck || cooldown)return false;

        /*
        HockeyistInfo info = *this, enemy = *enemyPuck;
        for(int i = 0; i < maxSwing; i++)
        {
            info.nextStep(0, 0);  enemy.nextStep(0, 0);
        }
        Vec2D delta = info.pos + 0.5 * stickLength * info.dir - enemy.pos;
        if(delta.sqr() < sqr(0.3 * stickLength))
        {
            plan.strikeTime = 0;  plan.swingTime = maxSwing;
            plan.targetIndex = -1;  plan.score = 0;  return true;
        }
        */
        if(inStrikeSector(*this, puckPath[0].pos))
        {
            plan.strikeTime = 0;  plan.swingTime = 0;
            plan.targetIndex = -1;  plan.score = 0;
            return activePlan = true;
        }
        return false;
    }

    void follow(const Vec2D &target, int &flags, double &turnEnd)
    {
        Vec2D delta = rotate(target - pos, conj(dir));
        double turn = atan2(delta.y, delta.x);
        flags = (turn > 0 ? MoveFlag::FIRST_LEFT : MoveFlag::FIRST_RIGHT);
        turnEnd = abs(turn) / turnAngle;
    }

    void tryAttack()
    {
        /*
        Vec2D target = puckPath[0].pos - defendPoint;
        double len = target.len();  target /= len;
        target = defendPoint + target * max(0.0, len - 5 * time);
        */
        Vec2D target = puckPath[0].pos + 0.5 * enemyPuck->spd / hockeyistFrict;
        plan.reset();  plan.secondTurnStart = maxLookahead;
        follow(target, plan.flags, plan.firstTurnEnd);
    }

    void fillMap()
    {
        if(swinging || knockdown)return;

        defend.score = 0;
        follow(defendPoint, defend.flags, defend.turnTime);
        attack.score = -numeric_limits<double>::infinity();

        stick.resize(gridSize);
        for(auto &flag : stick)flag = {0, -1};  activePlan = false;
        Mapper::fillMap(*this);
    }

    const MovePlan &chooseMove()
    {
        if(swinging || knockdown)return plan;  plan.evaluate<true>(*this);
        plan = Optimizer<true>().findBestMove(*this, plan);  return plan;
    }

    void execute(Move &move)
    {
        if(swinging && ((plan.flags & MoveFlag::STOP) || plan.strikeTime || plan.swingTime < 0))
        {
            move.setSpeedUp(0);  move.setTurn(0);  move.setAction(CANCEL_STRIKE);
        }
        else plan.execute(move);

        /*
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
        */
    }
};

struct EnemyInfo : public HockeyistInfo
{
    long long id;
    vector<Vec2D> rot;
    MovePlan plan;
    int swinging;

    EnemyInfo(const Hockeyist& hockeyist) : id(hockeyist.getId()), rot(maxLookahead + 1)
    {
        set(hockeyist);
    }

    void set(const Hockeyist& hockeyist)
    {
        HockeyistInfo::set(hockeyist, rot.data());  swinging = hockeyist.getSwingTicks();
    }

    const MovePlan &chooseMove()
    {
        if(swinging || knockdown)return plan;  plan.evaluate<false>(*this);
        plan = Optimizer<false>().findBestMove(*this, plan);  return plan;
    }

    void execute()
    {
        puckPath[0].set(*this, 0);  int time = 1;
        HockeyistInfo info = *this;  MovePlan::Helper helper(plan);
        for(; time <= plan.strikeTime; time++)
        {
            info.nextStep(helper.accel(time - 1), helper.turn(time - 1), time);
            puckPath[time].set(info, time);
        }
        for(; time <= plan.strikeTime + plan.swingTime; time++)
        {
            info.nextStep(0, Vec2D(1, 0), time);  puckPath[time].set(info, time);
        }

        Vec2D dir = info.dir;  double power;
        if(plan.swingTime < 0)
        {
            dir = rotate(dir, plan.passDir);  power = passBase * plan.passPower;
        }
        else power = strikeBase + strikeGrowth * plan.swingTime;
        power += info.spd * dir;

        PuckInfo puck(puckPath[time - 1].pos, power * dir);  int flag = 0;
        for(; time <= maxLookahead; time++)
        {
            if((flag = puck.nextStep()))break;  puckPath[time].set(puck, time);
        }
        puckPathLen = time;  goalFlag = leftPlayer ? flag : -flag;  plan.skipTurn();
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
        goalieTime = (totalGoals ? numeric_limits<int>::max() : world.getTickCount() - globalTick);
        puckPathLen = 0;  goalFlag = 0;

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

        enemyMap.reset();  EnemyInfo *enemyPuck = nullptr;
        for(auto enemy : enemies)
            if(enemy->id == owner)enemyPuck = &*enemy;
            else enemyMap.fillMap(*enemy);
        ::enemyPuck = enemyPuck;  enemyMap.filter();

        if(owner < 0)
        {
            puckPath[0].set(puck, 0);
            double intercept = 0;  int time = 1, flag = 0;
            for(; time <= maxLookahead; time++)
            {
                if((flag = puck.nextStep()))break;
                intercept = max(intercept, interceptProbability(puck, time));
                puckPath[time].set(puck, time, intercept);
            }
            puckPathLen = time;  goalFlag = leftPlayer ? flag : -flag;
        }
        if(enemyPuck)
        {
            enemyPuck->chooseMove();  enemyPuck->execute();
        }

        targets.clear();  AllyInfo *allyPuck = nullptr;
        for(auto ally : allies)
            if(ally->id == owner)
            {
                allyPuck = &*ally;  ally->activePlan = true;
            }
            else ally->fillMap();
        chooseBest(targets, 8);

        if(!allyPuck && goalFlag <= 0)
        {
            AllyInfo *follower = nullptr;
            double best = enemyPuck ? 1 : 0;
            for(auto ally : allies)
            {
                if(ally->tryKnockdown())continue;  double score = ally->chooseMove().score;
                if(!(score > best))continue;  best = score;  follower = &*ally;
            }
            if(follower)follower->activePlan = true;
        }

        long long defender = -1;
        double best = -numeric_limits<double>::infinity();
        for(auto ally : allies)if(!ally->activePlan)
        {
            if(!(ally->defend.score > best))continue;
            best = ally->defend.score;  defender = ally->id;
        }

        int tg = -1;
        if(allyPuck)tg = allyPuck->chooseMove().targetIndex;
        long long pass = (tg < 0 ? -1 : targets[tg].id);
        for(auto ally : allies)if(!ally->activePlan)
        {
            if(ally->id == pass)ally->plan.set(targets[tg]);
            else if(!enemyPuck || ally->id == defender)ally->plan.set(ally->defend);
            else ally->tryAttack();
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
