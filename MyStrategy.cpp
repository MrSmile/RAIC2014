#include "MyStrategy.h"

#include <cmath>
#include <cstdlib>
#include <cassert>
#include <algorithm>
#include <limits>
#include <map>

//#define PRINT_LOG
//#define CHECK_PREDICTION

#include <iostream>  // DEBUG
#include <iomanip>  // DEBUG

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



enum PuckStatus
{
    HOLD, FLY, PUCK_STATUS_COUNT
};

constexpr int maxLookahead = 256;
constexpr int turnStep = 4, strikeStep = 4, strikeDelta = 8, targetCount = 2;
constexpr int optSurvive = 4, optOffspring = 4, optStepCount = 4, optFinalCount = 4;
constexpr double optStepBase = 8, optStepMul = 0.5;

constexpr double cellSize = 32;
constexpr double dangerMultiplier = 10;
constexpr double minPassDist = 256, passMultiplier = 0.2;
constexpr double timeGamma = 1 - 1.0 / maxLookahead;

constexpr double hockeyistFrict = 0.02, angularFrict = 0.0270190131;
constexpr double puckFrict = 0.001, puckBeta = -log(1 - puckFrict);
constexpr double wallBounce = 0.25, hockeyistBounce = 0.325, goalieFrict = 0.1;
constexpr double minDepth = 0.01, depthFactor = 0.8;

double rinkLeft, rinkRight, rinkTop, rinkBottom;
double rinkWidth, rinkHeight;  Vec2D rinkCenter;
int gridHalfWidth, gridHalfHeight, gridLine, gridHeight, gridCenter, gridSize;
double goalCenter, goalHalf, hockeyistRad, puckRad, goalieSpd, goalieRange;
double staminaBase, staminaGrowth, accelMin, accelMax, turnAngle;
double stickLength, stickSectorTan, holdDist, timeGammaSwing;  int maxSwing;
double strikeBase, strikeGrowth, passBase, strikeBeta, passBeta;  Vec2D passSector;
double minChance, maxChance, knockdownChance, chanceDrop;
double strikeChance, pickChanceDelta[PUCK_STATUS_COUNT];

int globalTick = -1, goalieTime;
int defendCell, attackCell[2];
bool leftPlayer;


void initConsts(const model::Game& game, const model::World& world)
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

    staminaBase = 0.01 * game.getZeroStaminaHockeyistEffectivenessFactor();
    staminaGrowth = (0.01 - staminaBase) / game.getHockeyistMaxStamina();
    accelMin = -game.getHockeyistSpeedDownFactor();  accelMax = game.getHockeyistSpeedUpFactor();
    turnAngle = game.getHockeyistTurnAngleFactor();

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
    knockdownChance = game.getKnockdownChanceFactor();
    chanceDrop = 1 / game.getStruckPuckInitialSpeedFactor();
    strikeChance = game.getStrikePuckBaseChance();
    pickChanceDelta[HOLD] = strikeChance - game.getTakePuckAwayBaseChance();
    pickChanceDelta[FLY]  = strikeChance - game.getPickUpPuckBaseChance();


    leftPlayer = (2 * world.getMyPlayer().getNetBack() < rinkLeft + rinkRight);
    int defHOffs = lround(0.40 * rinkWidth  / cellSize);
    int attHOffs = lround(0.00 * rinkWidth  / cellSize);
    int attVOffs = lround(0.05 * rinkHeight / cellSize);
    if(leftPlayer)
    {
        defHOffs = -defHOffs;  attHOffs = -attHOffs;
    }
    defendCell = gridCenter + defHOffs;
    attackCell[0] = gridCenter - attHOffs - attVOffs * gridLine;
    attackCell[1] = gridCenter - attHOffs + attVOffs * gridLine;
}

inline double toChance(double chance)
{
    return max(minChance, min(maxChance, chance));
}

inline int gridPos(const Vec2D &pos)
{
    Vec2D offs = (pos - rinkCenter) / cellSize + Vec2D(gridHalfWidth + 0.5, gridHalfHeight + 0.5);
    return int(offs.y) * gridLine + int(offs.x);
}


struct UnitState
{
    Vec2D pos, spd;

    void set(const model::Unit& unit)
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

struct HockeyistInfo
{
    long long id;
    double accelMin, accelMax, turnAngle, maxTurnSteps;
    int maxTurnCount;  Vec2D turnRot, endTurnRot;
    double strikeBase, strikeGrowth, passBase, strikeBeta, passBeta;
    double attack[PUCK_STATUS_COUNT], defence;
    double knockupChance, knockdownChance;

    bool havePuck;
    int knockdown, cooldown;
    vector<Vec2D> rot;

    HockeyistInfo(const model::Hockeyist& hockeyist) : id(hockeyist.getId())
    {
    }

    void set(const model::Hockeyist& hockeyist, long long puckOwner)
    {
        double mul = staminaBase + staminaGrowth * hockeyist.getStamina();
        double str = mul * hockeyist.getStrength();
        double dex = mul * hockeyist.getDexterity();
        double end = mul * hockeyist.getEndurance();
        double agi = mul * hockeyist.getAgility();

        accelMin = ::accelMin * agi;  accelMax = ::accelMax * agi;
        turnAngle = ::turnAngle * agi;  maxTurnSteps = (pi / 2) / turnAngle;  maxTurnCount = int(maxTurnSteps);
        turnRot = sincos(turnAngle);  endTurnRot = sincos(turnAngle * (maxTurnSteps - maxTurnCount));

        strikeBase = ::strikeBase * str;  strikeGrowth = ::strikeGrowth * str;  passBase = ::passBase * str;
        strikeBeta = ::strikeBeta * dex;  passBeta = ::passBeta * dex;

        attack[HOLD] = ::strikeChance + max(str, dex);
        attack[FLY] = ::strikeChance + max(dex, agi);  defence = max(end, agi);
        knockupChance = toChance(::strikeChance - max(end, agi) + 1);
        knockdownChance = toChance(::knockdownChance - end + 1);


        havePuck = (id == puckOwner);
        knockdown = hockeyist.getRemainingKnockdownTicks();
        cooldown = max(knockdown, hockeyist.getRemainingCooldownTicks());

        double angSpd = hockeyist.getAngularSpeed();
        rot.resize(maxLookahead + 1);  rot[0] = sincos(angSpd);
        for(int i = 1; i <= maxLookahead; i++)
        {
            angSpd -= angularFrict * angSpd;  rot[i] = sincos(angSpd);
        }
    }
};

struct HockeyistState : public UnitState
{
    Vec2D dir;

    void set(const model::Hockeyist& hockeyist)
    {
        UnitState::set(hockeyist);
        dir = sincos(hockeyist.getAngle());
    }

    void nextStep(const HockeyistInfo &info, double accel, const Vec2D &turn, int time)
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

        pos += spd;  dir = rotate(dir, info.rot[time]);  dir = rotate(dir, turn);
    }
};

struct SimplePuckState : public UnitState
{
    double spdLen;

    SimplePuckState(const Vec2D &pos_, const Vec2D &spd_)
    {
        pos = pos_;  spd = spd_;  spdLen = spd.len();
    }

    bool nextStep()
    {
        spd -= puckFrict * spd;  spdLen -= puckFrict * spdLen;  pos += spd;
        if(abs(pos.x - rinkCenter.x) > rinkWidth / 2 + cellSize)return false;

        double delta;
        if((delta = rinkTop + puckRad - pos.y) > 0)
        {
            if(spd.y < 0)
            {
                spd.y *= -wallBounce;  spdLen = spd.len();
            }
            if(spd.y < delta)pos.y += depthFactor * (delta - minDepth);
        }
        if((delta = rinkBottom - puckRad - pos.y) < 0)
        {
            if(spd.y > 0)
            {
                spd.y *= -wallBounce;  spdLen = spd.len();
            }
            if(spd.y > delta)pos.y += depthFactor * (delta + minDepth);
        }
        return true;
    }
};

struct PuckState : public UnitState
{
    Vec2D goalie[2];

    PuckState() = default;

    PuckState(const Vec2D &pos_, const Vec2D &spd_)
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

inline Vec2D puckPos(const HockeyistInfo &info, const HockeyistState &state, int time)
{
    return state.pos + hockeyistFrict * state.spd + holdDist * rotate(state.dir, conj(info.rot[time]));
}

struct PathPoint
{
    Vec2D pos;  PuckStatus status;
    double intercept, chanceDrop;

    void set(const PuckState &state, int time, double intercept_ = 0)
    {
        pos = state.pos;  status = FLY;  intercept = intercept_;
        chanceDrop = ::chanceDrop * state.spd.len();
    }

    void set(const HockeyistInfo &info, const HockeyistState &state, int time)
    {
        pos = puckPos(info, state, time);  status = HOLD;  intercept = 0;
        chanceDrop = info.defence;
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

template<typename T, typename State> using AddPointFunc =
    void (T::*)(const HockeyistInfo &info, State &state, int time, int flags, double turnTime);
template<typename T, typename State> using ClosePathFunc =
    void (T::*)(const HockeyistInfo &info, State &state, int flags, double turnTime);

template<typename T, typename State, AddPointFunc<T, State> addPoint, ClosePathFunc<T, State> closePath> struct Mapper
{
    T &obj;
    const HockeyistInfo &info;

    Mapper(T &obj_, const HockeyistInfo &info_) : obj(obj_), info(info_)
    {
    }

    void fillMapTail(const State &state, int flags, double accel, int time)
    {
        State cur = state;
        for(int i = time + 1; i <= maxLookahead; i++)
        {
            cur.nextStep(info, accel, Vec2D(1, 0), i);  (obj.*addPoint)(info, cur, i, flags, time);
        }
        (obj.*closePath)(info, cur, flags, time);
    }

    void fillMapTail(const State &state, int flags, double accel, const Vec2D &turn, const Vec2D &endTurn, int time)
    {
        double turnTime = time + info.maxTurnSteps;
        State cur = state;  flags ^= MoveFlag::BACK;
        int i = time, end = time + info.maxTurnCount;
        for(i++; i <= maxLookahead; i++)
        {
            if(i > end)
            {
                cur.nextStep(info, accel, endTurn, i);  (obj.*addPoint)(info, cur, i, flags, turnTime);  break;
            }
            cur.nextStep(info, accel, turn, i);  (obj.*addPoint)(info, cur, i, flags, i);
        }
        for(i++; i <= maxLookahead; i++)
        {
            cur.nextStep(info, accel, Vec2D(1, 0), i);  (obj.*addPoint)(info, cur, i, flags, turnTime);
        }
        (obj.*closePath)(info, cur, flags, turnTime);
    }

    void fillMap(const State& state, int flags)
    {
        double accel1 = info.accelMax, accel2 = info.accelMin;
        if(flags & MoveFlag::BACK)swap(accel1, accel2);
        Vec2D turn = info.turnRot, endTurn = info.endTurnRot;
        if(flags & MoveFlag::FIRST_RIGHT)
        {
            turn.y = -turn.y;  endTurn.y = -endTurn.y;
        }

        State cur = state;
        int n = min(maxLookahead, info.knockdown + info.maxTurnCount);
        for(int i = info.knockdown + 1; i <= n; i++)
        {
            cur.nextStep(info, accel1, turn, i);  (obj.*addPoint)(info, cur, i, flags, i);  if(i % turnStep)continue;  // TODO: globalTick ?
            fillMapTail(cur, flags, accel1, i);  fillMapTail(cur, flags, accel2, turn, endTurn, i);
        }
    }

    void fillMap(const State& state)
    {
        State cur = state;  (obj.*addPoint)(info, cur, 0, 0, 0);
        for(int i = 1; i <= info.knockdown; i++)
        {
            cur.nextStep(info, 0, Vec2D(1, 0), i);  (obj.*addPoint)(info, cur, i, 0, 0);
        }
        fillMapTail(cur, MoveFlag::BACK, info.accelMin, info.knockdown);
        fillMapTail(cur, MoveFlag::FORW, info.accelMax, info.knockdown);
        fillMap(cur, MoveFlag::BACK | MoveFlag::FIRST_LEFT);
        fillMap(cur, MoveFlag::BACK | MoveFlag::FIRST_RIGHT);
        fillMap(cur, MoveFlag::FORW | MoveFlag::FIRST_LEFT);
        fillMap(cur, MoveFlag::FORW | MoveFlag::FIRST_RIGHT);
    }
};

struct ReachabilityMap
{
    vector<int> body, stick;

    void reset()
    {
        body.resize(gridSize);  stick.resize(gridSize);
        for(auto &flag :  body)flag = maxLookahead;
        for(auto &flag : stick)flag = maxLookahead;
    }

    void addMapPoint(const HockeyistInfo &info, HockeyistState &state, int time, int flags, double turnTime)
    {
        auto &cell1 = body[gridPos(state.pos)];  cell1 = min<int>(cell1, time);
        auto &cell2 = stick[gridPos(state.pos + 0.5 * stickLength * state.dir)];
        cell2 = min<int>(cell2, max(time, info.cooldown));
    }

    void closePath(const HockeyistInfo &info, HockeyistState &state, int flags, double turnTime)
    {
    }

    void fillMap(const HockeyistInfo &info, const HockeyistState &state)
    {
        Mapper<ReachabilityMap, HockeyistState, &ReachabilityMap::addMapPoint, &ReachabilityMap::closePath>(*this, info).fillMap(state);
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


struct PassTarget;
struct EnemyInfo;
struct AllyInfo;

ReachabilityMap enemyMap;
PathPoint puckPath[maxLookahead + 1];
vector<PassTarget> targets, nextTargets;
int puckPathLen, goalFlag;
EnemyInfo *enemyPuck;
AllyInfo *allyPuck;


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

double checkSafety(const HockeyistInfo &info, const HockeyistState &state, int time)
{
    if(checkGoalie(state.pos, hockeyistRad))return 0;

    double res = 0;
    if(!puckPathLen)
    {
        Vec2D puck = puckPos(info, state, time);  if(!validPuckPos(puck))return 0;
        if(enemyMap.stick[gridPos(state.pos)] <= time)res = max(res, info.knockdownChance);
        if(enemyMap.stick[gridPos(puck)] <= time)res = max(res, info.knockupChance);
    }
    return 1 - res;
}

template<bool ally> struct StateScore : public HockeyistState
{
    double timeFactor;

    StateScore(const HockeyistState &state) : HockeyistState(state), timeFactor(1)
    {
    }

    void nextStep(const HockeyistInfo &info, double accel, const Vec2D &turn, int time)
    {
        HockeyistState::nextStep(info, accel, turn, time);  timeFactor *= timeGamma;
    }

    double multiplier() const
    {
        return timeFactor;
    }
};

template<> struct StateScore<true> : public HockeyistState
{
    double timeFactor, safety;

    StateScore(const HockeyistState &state) : HockeyistState(state), timeFactor(1), safety(1)
    {
    }

    void nextStep(const HockeyistInfo &info, double accel, const Vec2D &turn, int time)
    {
        HockeyistState::nextStep(info, accel, turn, time);  timeFactor *= timeGamma;
        safety = min(safety, checkSafety(info, *this, time));
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

double interceptProbability(const Vec2D &pos, double spd, int time)
{
    if(enemyMap.stick[gridPos(pos)] > time)return 0;
    return min(maxChance, 1 + strikeChance - chanceDrop * spd);
}

double checkInterception(const Vec2D &pos, const Vec2D &spd, int time, int duration)
{
    SimplePuckState puck(pos, spd);  double res = 0;
    for(int i = 1; i <= duration; i++)
    {
        if(!puck.nextStep())break;
        res = max(res, interceptProbability(puck.pos, puck.spdLen, time + i));
    }
    return res;
}

inline bool inStrikeSector(const HockeyistState &state, const Vec2D &puck)
{
    Vec2D delta = puck - state.pos;
    if(!(delta.sqr() < sqr(stickLength)))return false;
    double dot = delta * state.dir, cross = delta % state.dir;
    return abs(cross) < dot * stickSectorTan;
}

struct StrikeInfo
{
    int strikeTime, swingTime, targetIndex;
    double score, passPower;  Vec2D passDir, strikePos;

    StrikeInfo() : strikeTime(maxLookahead), swingTime(0), targetIndex(-1), score(0), passPower(1), passDir{0, 0}
    {
    }

    StrikeInfo(int strikeTime_, double score_) :
        strikeTime(strikeTime_), swingTime(0), targetIndex(-1), score(score_), passPower(1), passDir{0, 0}
    {
    }

    void reset()
    {
        *this = StrikeInfo();
    }

    void tryStrikeFlyby(double strike, int time, double val)
    {
        double pick = strike - pickChanceDelta[puckPath[time].status];

        int swing = -1;
        if(goalFlag)  // TODO: check strike dir
        {
            val *= toChance(strike) * dangerMultiplier;
            if(pick < maxChance)swing = 0;
        }
        else
        {
            val *= toChance(pick) * passMultiplier;
            int deltaTime = time - enemyMap.stick[gridPos(puckPath[time].pos)];
            val *= (1 - puckPath[time].intercept) * (1 + deltaTime / maxLookahead);
        }
        if(!(val > score))return;

        strikeTime = time;  swingTime = swing;  targetIndex = -1;  score = val;
    }

    template<bool ally> void evaluateStrike(const HockeyistState &state, const Vec2D &puck,
        double power, double beta, int swing, int time, double val)
    {
        Sector res = estimateGoalAngle(puck, state.spd, power, ally == leftPlayer);
        if(!(res.span > 0))return;

        Vec2D delta = rotate(res.dir, conj(state.dir)), offs(delta.x, abs(delta.y));
        if(swing < 0)offs = (offs.x > passSector.x ? Vec2D(1, 0) : rotate(offs, conj(passSector)));
        if(!(offs.x > 0))return;

        val *= erf(beta * (offs.y + res.span)) - erf(beta * (offs.y - res.span));
        if(!((val *= 0.5) > score))return;

        if(ally)
        {
            Vec2D dir = res.dir;  // TODO: better approx
            val *= 1 - checkInterception(puck, (power + state.spd * dir) * dir, time, lround(res.time));
            if(!(val > score))return;
        }

        strikeTime = time;  swingTime = swing;  targetIndex = -1;
        score = val;  passPower = 1;  passDir = delta;  strikePos = puck;
    }

    void evaluatePass(const HockeyistState &state, const Vec2D &puck, double power,
        double x, double y, double spd, int catchTime, int time, double val, int target)
    {
        Vec2D delta = Vec2D(x, y) - puck;
        double len2 = delta.sqr();  if(len2 < sqr(minPassDist))return;
        double len = sqrt(len2);  delta /= len;
        Vec2D offs = rotate(delta, conj(state.dir));
        if(!(offs.x > passSector.x))return;

        double puckSpd = puckBeta * len / (1 - exp(-puckBeta * (catchTime - time)));
        double relSpd = state.spd * delta;  if(!(puckSpd < power + relSpd))return;
        val *= maxChance - max(0.0, puckSpd - spd - puckBeta * len) * chanceDrop;
        if(!(val > score))return;

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
        score = val;  passPower = (puckSpd - relSpd) / power;  passDir = offs;
    }

    void tryPass(const HockeyistState &state, const Vec2D &puck, double power, int time, double mul);

    template<bool ally> void tryStrike(const HockeyistInfo &info, const StateScore<ally> &state, int time)
    {
        double mul = state.multiplier();  Vec2D puck;
        if(time + maxSwing <= maxLookahead && mul * timeGammaSwing > score)do
        {
            StateScore<ally> cur = state;
            for(int i = 1; i <= maxSwing; i++)cur.nextStep(info, 0, Vec2D(1, 0), time + i);

            double val = cur.multiplier();
            if(!puckPathLen)puck = puckPos(info, cur, time + maxSwing);
            else
            {
                if(!inStrikeSector(cur, puck = puckPath[time + maxSwing].pos))break;
                double strike = info.attack[puckPath[time + maxSwing].status] - puckPath[time + maxSwing].chanceDrop;
                val *= toChance(strike) * (1 - puckPath[time + maxSwing].intercept);
            }
            if(val > score)evaluateStrike<ally>(cur, puck,
                info.strikeBase + info.strikeGrowth * maxSwing, info.strikeBeta, maxSwing, time, val);
        }
        while(false);

        double strike = 0;
        if(!puckPathLen)puck = puckPos(info, state, time);
        else
        {
            if(!inStrikeSector(state, puck = puckPath[time].pos))return;
            strike = info.attack[puckPath[time].status] - puckPath[time].chanceDrop;
            mul *= toChance(strike) * (1 - puckPath[time].intercept);
        }
        if(mul > score)
        {
            if(info.havePuck)
                 evaluateStrike<ally>(state, puck, info.passBase,   info.passBeta,  -1, time, mul);
            else evaluateStrike<ally>(state, puck, info.strikeBase, info.strikeBeta, 0, time, mul);
        }
        if(!ally)return;

        if(!puckPathLen)tryPass(state, puck, info.passBase, time, mul);
        else if(time < puckPathLen)tryStrikeFlyby(strike, time, mul);
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

        Helper(const HockeyistInfo &info, const MovePlan &plan) : pos(0)
        {
            accBase = info.accelMin;  accDelta = info.accelMax;
            if(plan.flags & MoveFlag::BACK)swap(accBase, accDelta);  accDelta -= accBase;
            accTime2 = plan.secondTurnStart + info.maxTurnSteps;
            accTime1 = min(plan.firstTurnEnd - info.maxTurnSteps, accTime2);

            double turnTime2 = plan.secondTurnStart;
            double turnTime1 = min(plan.firstTurnEnd, turnTime2);
            turnStep[0] = (int)floor(turnTime1);  turnTime1 -= turnStep[0];
            turnStep[1] = (int)floor(turnTime2);  turnTime2 -= turnStep[1];
            turnStep[2] = maxLookahead;

            rot[0] = info.turnRot;  rot[4] = info.turnRot;
            double angle1 = info.turnAngle * turnTime1;
            double angle2 = info.turnAngle * turnTime2;
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

    MovePlan() : StrikeInfo(0, 0), flags(MoveFlag::STOP), firstTurnEnd(0), secondTurnStart(0)
    {
    }

    MovePlan(int flags_, double turnTime, int strikeTime, double score) : StrikeInfo(strikeTime, score),
        flags(flags_ | MoveFlag::STOP), firstTurnEnd(turnTime), secondTurnStart(strikeTime)
    {
    }

    MovePlan(const StrikeInfo &info, int flags_, double turnTime) : StrikeInfo(info),
        flags(flags_), firstTurnEnd(turnTime), secondTurnStart(info.strikeTime)
    {
    }

    template<bool ally> void evaluate(const HockeyistInfo &info, const HockeyistState &state)
    {
        int minStrike = max(strikeTime - strikeDelta, info.cooldown);
        int maxStrike = min(strikeTime + strikeDelta, maxLookahead);
        reset();

        StateScore<ally> cur = state;  Helper helper(info, *this);
        for(int i = 0;; i++)
        {
            if(i >= minStrike)tryStrike<ally>(info, cur, i);  if(i >= maxStrike)break;
            cur.nextStep(info, helper.accel(i), helper.turn(i), i + 1);
        }
    }

    static int mutate(double &val, double step, double minVal, double maxVal)
    {
        val += step * frand();
        if(val < minVal)
        {
            val = min(maxVal, 2 * minVal - val);  return -1;
        }
        if(val > maxVal)
        {
            val = max(minVal, 2 * maxVal - val);  return +1;
        }
        return 0;
    }

    MovePlan(const HockeyistInfo &info, const HockeyistState &state, const MovePlan &old, double step, bool ally) :
        flags(old.flags), firstTurnEnd(old.firstTurnEnd), secondTurnStart(old.secondTurnStart)
    {
        strikeTime = old.strikeTime;  int minTurn2 = max(0, strikeTime - info.maxTurnCount);
        if(mutate(firstTurnEnd, step, 0, 2 * info.maxTurnSteps) < 0)flags ^= MoveFlag::FIRST_RIGHT;
        if(mutate(secondTurnStart, 4 * step, minTurn2, strikeTime) > 0)flags ^= MoveFlag::SECOND_RIGHT;
        if(ally)evaluate<true>(info, state);  else evaluate<false>(info, state);
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

    void print(const char *prefix, long long id)
    {
        cout << prefix << id << ": " << score << ", flags: " << flags;
        cout << ", times: " << firstTurnEnd << '\'' << secondTurnStart << '\'' << strikeTime;
        if(!(flags & MoveFlag::STOP))
        {
            if(swingTime >= 0)cout << ", swing: " << swingTime;
            else if(puckPathLen)cout << ", pickup";
            else cout << ", pass: " << passDir.x << '\'' << passDir.y << ", power: " << passPower;
            if(targetIndex >= 0)cout << ", target: " << targetIndex;
        }
        cout << endl;
    }

    void execute(const HockeyistInfo &info, model::Move &move)
    {
        if(!strikeTime)
        {
            move.setSpeedUp(0);  move.setTurn(0);
            if(flags & MoveFlag::STOP)move.setAction(model::NONE);
            else if(swingTime > 0)
            {
                move.setAction(model::SWING);  swingTime--;  return;
            }
            else if(!swingTime)move.setAction(model::STRIKE);
            else if(puckPathLen)move.setAction(model::TAKE_PUCK);
            else
            {
                move.setPassPower(passPower);
                move.setPassAngle(atan2(passDir.y, passDir.x));  move.setAction(model::PASS);
            }
            flags = MoveFlag::STOP;  return;
        }

        Helper helper(info, *this);
        double accel = helper.accel(0);  Vec2D rot = helper.turn(0);
        move.setSpeedUp(accel > 0 ? accel / info.accelMax : -accel / info.accelMin);
        move.setTurn(atan2(rot.y, rot.x));  move.setAction(model::NONE);

        firstTurnEnd = max(0.0, firstTurnEnd - 1);
        secondTurnStart = max(0.0, secondTurnStart - 1);
        strikeTime--;
    }
};

struct PassTarget : public Vec2D, public MovePlan
{
    int id;  double spd;

    PassTarget() = default;

    PassTarget(const HockeyistInfo &info, const MovePlan &plan) : Vec2D(plan.strikePos), MovePlan(plan), id(info.id)
    {
        spd = (info.attack[FLY] - pickChanceDelta[FLY] - maxChance) / chanceDrop;
    }

    PassTarget(const HockeyistInfo &info, const Vec2D &puck, int flags, double turnTime, int strikeTime, double score) :
        Vec2D(puck), MovePlan(flags, turnTime, strikeTime, score), id(info.id)
    {
        spd = (info.attack[FLY] - pickChanceDelta[FLY] - maxChance) / chanceDrop;
    }

    bool operator < (const PassTarget &cmp) const
    {
        return score < cmp.score;
    }
};

void StrikeInfo::tryPass(const HockeyistState &state, const Vec2D &puck, double power, int time, double mul)
{
    double y1 = rinkTop + puckRad, y2 = rinkBottom - puckRad;
    for(size_t i = 0; i < targets.size(); i++)
    {
        int hit = targets[i].strikeTime + targets[i].swingTime;  if(hit <= time)continue;
        double val = mul * targets[i].score;  if(!(val > score))continue;

        evaluatePass(state, puck, power, targets[i].x, targets[i].y, targets[i].spd, hit, time, val, i);
        evaluatePass(state, puck, power, targets[i].x, y1 - (targets[i].y - y1) / wallBounce, targets[i].spd, hit, time, val, i);
        evaluatePass(state, puck, power, targets[i].x, y2 - (targets[i].y - y2) / wallBounce, targets[i].spd, hit, time, val, i);
    }
}


template<bool ally> struct OptimizerState : public StateScore<ally>, public StrikeInfo
{
    OptimizerState(const HockeyistState &state) : StateScore<ally>(state)
    {
    }
};

template<typename T> size_t chooseBest(vector<T> &array, size_t max)
{
    auto rend = array.rend();  make_heap(array.rbegin(), rend);  size_t n = 0;
    for(; n < max && rend != array.rbegin(); n++, --rend)pop_heap(array.rbegin(), rend);
    array.resize(n);  return n;
}

template<bool ally> struct Optimizer
{
    typedef OptimizerState<ally> State;

    vector<MovePlan> moves;

    void reset()
    {
        moves.clear();
    }

    void addPoint(const HockeyistInfo &info, State &state, int time, int flags, double turnTime)
    {
        if(time < info.cooldown)return;
        if(time < strikeDelta || !(time % strikeStep))state.StrikeInfo::tryStrike<ally>(info, state, time);  // TODO: globalTick ?
    }

    void closePath(const HockeyistInfo &info, State &state, int flags, double turnTime)
    {
        moves.emplace_back(state, flags, turnTime);
    }

    void optimizeMove(const HockeyistInfo &info, const HockeyistState &state, double step, int survive, int offspring)
    {
        size_t n = chooseBest(moves, survive);
        for(size_t i = 0; i < n; i++)for(int j = 0; j < offspring; j++)
            moves.emplace_back(info, state, moves[i], step, ally);
    }

    void fillMap(const HockeyistInfo &info, const HockeyistState &state)
    {
        typedef Optimizer<ally> Self;  assert(!moves.size());
        Mapper<Self, State, &Self::addPoint, &Self::closePath>(*this, info).fillMap(state);
    }

    void findBestMove(const HockeyistInfo &info, const HockeyistState &state, MovePlan &plan)
    {
        assert(moves.size());  plan.evaluate<ally>(info, state);  moves.push_back(plan);

        double delta = optStepBase;
        for(int i = 0; i < optStepCount; i++, delta *= optStepMul)
            optimizeMove(info, state, delta, optSurvive, optOffspring);
        if(ally)for(int i = 0; i < optFinalCount; i++, delta *= optStepMul)
            optimizeMove(info, state, delta, 1, optOffspring);

        const MovePlan *res = &plan;
        double best = -numeric_limits<double>::infinity();
        for(auto &move : moves)if(move.score > best)
        {
            res = &move;  best = move.score;
        }
        plan = *res;
    }
};

struct EnemyInfo : public HockeyistInfo, public HockeyistState, public Optimizer<false>
{
    MovePlan plan;
    int swinging;

    EnemyInfo(const model::Hockeyist& hockeyist, long long puckOwner) : HockeyistInfo(hockeyist)
    {
        set(hockeyist, puckOwner);
    }

    void set(const model::Hockeyist& hockeyist, long long puckOwner)
    {
        HockeyistInfo::set(hockeyist, puckOwner);  HockeyistState::set(hockeyist);
        reset();  swinging = hockeyist.getSwingTicks();
    }

    const MovePlan &choosePuckMove()
    {
        if(swinging)return plan;  Optimizer<false>::fillMap(*this, *this);
        findBestMove(*this, *this, plan);  return plan;
    }

    void execute()
    {
        puckPath[0].set(*this, *this, 0);  int time = 1;
        HockeyistState state = *this;  MovePlan::Helper helper(*this, plan);
        for(; time <= plan.strikeTime; time++)
        {
            state.nextStep(*this, helper.accel(time - 1), helper.turn(time - 1), time);
            puckPath[time].set(*this, state, time);
        }
        for(; time <= plan.strikeTime + plan.swingTime; time++)
        {
            state.nextStep(*this, 0, Vec2D(1, 0), time);  puckPath[time].set(*this, state, time);
        }

        Vec2D dir = state.dir;  double power;
        if(plan.swingTime < 0)
        {
            dir = rotate(dir, plan.passDir);  power = passBase * plan.passPower;
        }
        else power = strikeBase + strikeGrowth * plan.swingTime;
        power += state.spd * dir;

        PuckState puck(puckPath[time - 1].pos, power * dir);  int flag = 0;
        for(; time <= maxLookahead; time++)
        {
            if((flag = puck.nextStep()))break;  puckPath[time].set(puck, time);
        }
        puckPathLen = time;  goalFlag = leftPlayer ? flag : -flag;  plan.skipTurn();
    }
};

struct AllyInfo : public HockeyistInfo, public HockeyistState, public Optimizer<true>
{
    MovePlan plan;  bool activePlan;
    vector<pair<double, int>> stick;
    PassTarget defend, attack;
    Vec2D defendPoint;
    int swinging;

    AllyInfo(const model::Hockeyist& hockeyist, long long puckOwner) : HockeyistInfo(hockeyist)
    {
        set(hockeyist, puckOwner);  double offs = rinkWidth / 2 - 2 * goalHalf;
        defendPoint = rinkCenter + Vec2D(leftPlayer ? -offs : offs, 0);
    }

    void set(const model::Hockeyist& hockeyist, long long puckOwner)
    {
#ifdef CHECK_PREDICTION
        Vec2D errPos = pos, errSpd = spd, errDir = dir;
#endif

        HockeyistInfo::set(hockeyist, puckOwner);  HockeyistState::set(hockeyist);
        reset();  swinging = hockeyist.getSwingTicks();

#ifdef CHECK_PREDICTION
        errPos -= pos;  errSpd -= spd;  errDir -= dir;
        if(abs(errPos.x) > 1e-3 || abs(errPos.y) > 1e-3 ||
           abs(errSpd.x) > 1e-4 || abs(errSpd.y) > 1e-4 ||
           abs(errDir.x) > 1e-5 || abs(errDir.y) > 1e-5)
        {
            cout << "Error: " << id << ' ';
            cout << errPos.x << ' ' << errPos.y << ' ';
            cout << errSpd.x << ' ' << errSpd.y << ' ';
            cout << errDir.x << ' ' << errDir.y << endl;
        }
#endif
    }

    void addPoint(const HockeyistInfo &info, State &state, int time, int flags, double turnTime)
    {
        Optimizer<true>::addPoint(info, state, time, flags, turnTime);

        Vec2D puck = state.pos + 0.5 * stickLength * state.dir;
        if(!validPuckPos(puck))return;  int cell = gridPos(puck);
        double score = state.safety * state.timeFactor * passMultiplier;
        score *= 1 - sqr(2 * (state.pos.x - rinkCenter.x) / rinkWidth);
        score *= 1 - sqr(2 * (state.pos.y - rinkCenter.y) / rinkHeight);
        int delta = enemyMap.stick[cell] - max(time, info.cooldown);
        score *= 0.5 + double(delta) * (0.5 / maxLookahead);
        if(!(score > stick[cell].first))return;
        stick[cell].first = score;

        PassTarget target(info, puck, flags, turnTime, time, score);
        if(cell == defendCell && score > defend.score)defend = target;
        if((cell == attackCell[0] || cell == attackCell[1]) && score > attack.score)attack = target;

        if(stick[cell].second < 0)
        {
            stick[cell].second = nextTargets.size();  nextTargets.push_back(target);
        }
        else nextTargets[stick[cell].second] = target;
    }

    void closePath(const HockeyistInfo &info, State &state, int flags, double turnTime)
    {
        Optimizer<true>::closePath(info, state, flags, turnTime);
    }

    bool tryKnockdown()
    {
        if(!enemyPuck || cooldown)return false;

        /*
        HockeyistState state = *this, enemy = *enemyPuck;
        for(int i = 1; i <= maxSwing; i++)
        {
            state.nextStep(*this, 0, Vec2D(1, 0), i);
            enemy.nextStep(*enemyPuck, 0, Vec2D(1, 0), i);
        }
        Vec2D delta = state.pos + 0.5 * stickLength * state.dir - enemy.pos;
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

    const MovePlan &chooseInterceptMove()
    {
        if(swinging || knockdown)return plan;

        attack.id = defend.id = -1;  attack.spd = defend.spd = 0;
        defend.reset();  defend.secondTurnStart = maxLookahead;
        follow(defendPoint, defend.flags, defend.firstTurnEnd);
        static_cast<MovePlan &>(attack) = MovePlan();

        stick.resize(gridSize);  for(auto &flag : stick)flag = {0, -1};  activePlan = false;
        Mapper<AllyInfo, State, &AllyInfo::addPoint, &AllyInfo::closePath>(*this, *this).fillMap(*this);
        findBestMove(*this, *this, plan);  return plan;
    }

    const MovePlan &choosePuckMove()
    {
        if(swinging)return plan;  Optimizer<true>::fillMap(*this, *this);
        findBestMove(*this, *this, plan);  return plan;
    }

    void execute(model::Move &move)
    {
        if(swinging && ((plan.flags & MoveFlag::STOP) || plan.strikeTime || plan.swingTime < 0))
        {
            move.setSpeedUp(0);  move.setTurn(0);  move.setAction(model::CANCEL_STRIKE);
        }
        else plan.execute(*this, move);

#ifdef PRINT_LOG
        switch(move.getAction())
        {
        case model::TAKE_PUCK:      cout << ">>>>>>>>>>>>> TAKE_PUCK"     << endl;  break;
        case model::SWING:          cout << ">>>>>>>>>>>>> SWING"         << endl;  break;
        case model::STRIKE:         cout << ">>>>>>>>>>>>> STRIKE"        << endl;  break;
        case model::CANCEL_STRIKE:  cout << ">>>>>>>>>>>>> CANCEL_STRIKE" << endl;  break;
        case model::PASS:           cout << ">>>>>>>>>>>>> PASS"          << endl;  break;
        case model::SUBSTITUTE:     cout << ">>>>>>>>>>>>> SUBSTITUTE"    << endl;  break;
        default:  break;
        }
        if(move.getAction() == model::STRIKE && !havePuck && !inStrikeSector(*this, puckPath[0].pos))
            cout << "----------------------- BAD STRIKE -----------------------" << endl;
#endif

#ifdef CHECK_PREDICTION
        double accel = (move.getSpeedUp() > 0 ? accelMax : -accelMin) * move.getSpeedUp();
        nextStep(*this, accel, sincos(move.getTurn()), 1);
#endif
    }
};

template<typename T> void synchronize(vector<T> &list, vector<T *> &ref, const model::Hockeyist &hockeyist, long long puckOwner)
{
    for(auto iter = list.begin();; ++iter)
        if(iter == list.end())
        {
            list.emplace_back(hockeyist, puckOwner);  ref.push_back(&*list.rbegin());  break;
        }
        else if(iter->id == hockeyist.getId())
        {
            iter->set(hockeyist, puckOwner);  ref.push_back(&*iter);  break;
        }
}

vector<EnemyInfo> enemyInfo;
vector<EnemyInfo *> enemies;
vector<AllyInfo> allyInfo;
vector<AllyInfo *> allies;


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


void MyStrategy::move(const model::Hockeyist& self, const model::World& world, const model::Game& game, model::Move& move)
{
    const auto &player = world.getMyPlayer();
    if(player.isJustMissedGoal() || player.isJustScoredGoal())
    {
        move.setTurn(self.getTeammateIndex() & 1 ? -pi : pi);
        move.setSpeedUp(0);  move.setAction(model::NONE);  return;
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
        PuckState puck;  puck.set(world.getPuck());  int totalGoals = 0;
        for(auto &player : world.getPlayers())totalGoals += player.getGoalCount();
        goalieTime = (totalGoals ? numeric_limits<int>::max() : world.getTickCount() - globalTick);
        puckPathLen = 0;  goalFlag = 0;

        allies.clear();  enemies.clear();
        for(auto &hockeyist : world.getHockeyists())if(hockeyist.getState() != model::RESTING)
        {
            if(hockeyist.getType() == model::GOALIE)
            {
                puck.goalie[hockeyist.isTeammate() ? 0 : 1] = Vec2D(hockeyist.getX(), hockeyist.getY());  continue;
            }
            if(hockeyist.isTeammate())
                 synchronize(allyInfo,  allies,  hockeyist, owner);
            else synchronize(enemyInfo, enemies, hockeyist, owner);
        }

        enemyMap.reset();  EnemyInfo *enemyPuck = nullptr;
        for(auto enemy : enemies)
            if(enemy->id == owner)enemyPuck = &*enemy;
            else enemyMap.fillMap(*enemy, *enemy);
        ::enemyPuck = enemyPuck;  enemyMap.filter();

        if(owner < 0)
        {
            puckPath[0].set(puck, 0);
            double intercept = 0;  int time = 1, flag = 0;
            for(; time <= maxLookahead; time++)
            {
                if((flag = puck.nextStep()))break;
                intercept = max(intercept, interceptProbability(puck.pos, puck.spd.len(), time));
                puckPath[time].set(puck, time, intercept);
            }
            puckPathLen = time;  goalFlag = leftPlayer ? flag : -flag;
        }
        if(enemyPuck)
        {
            enemyPuck->choosePuckMove();  enemyPuck->execute();
        }

        AllyInfo *allyPuck = nullptr, *follower = nullptr;
        targets.clear();  double best = enemyPuck ? 1 : 0;
        for(auto ally : allies)
            if(ally->id == owner)
            {
                allyPuck = &*ally;  ally->activePlan = true;
            }
            else
            {
                const MovePlan &plan = ally->chooseInterceptMove();  nextTargets.emplace_back(*ally, plan);
                if(!(plan.score > best))continue;  best = plan.score;  follower = &*ally;
            }
        if(!allyPuck && goalFlag <= 0 && follower)follower->activePlan = true;
        chooseBest(nextTargets, targetCount);  swap(targets, nextTargets);

        int tg = -1;
        if(allyPuck)tg = allyPuck->choosePuckMove().targetIndex;
        long long pass = (tg < 0 ? -1 : targets[tg].id);

        long long defender = -1;
        best = -numeric_limits<double>::infinity();
        for(auto ally : allies)if(!ally->activePlan)
        {
            if(!(ally->defend.score > best))continue;
            best = ally->defend.score;  defender = ally->id;
        }
        for(auto ally : allies)if(!ally->activePlan && !ally->tryKnockdown())
        {
            if(ally->id == pass)ally->plan = targets[tg];
            else if(ally->id == defender)ally->plan = ally->defend;
            else if(!enemyPuck)ally->plan = ally->attack;
            else ally->tryAttack();
        }

#ifdef PRINT_LOG
        if(follower && follower->activePlan)follower->plan.print("Follower ", follower->id);
        if(allyPuck)allyPuck->plan.print("Holder ", allyPuck->id);
        if(tg >= 0)targets[tg].print("Target ", pass);
#endif

        //drawField(enemyMap.body, tg < 0 ? -1 : gridPos(targets[tg]));
    }

    for(auto ally : allies)if(ally->id == self.getId())
    {
        ally->execute(move);  return;
    }
}

MyStrategy::MyStrategy()
{
    cout << fixed << setprecision(5);
}
