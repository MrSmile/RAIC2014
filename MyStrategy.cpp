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

constexpr int MAX_ACTIVE = 3;

constexpr int maxLookahead = 160, stickLookahead = 20;
constexpr int turnStep = 4, strikeStep = 4, strikeDelta = 8, targetCount = 2;
constexpr int optSurvive = 4, optOffspring = 4, optStepCount = 4, optFinalCount = 4;
constexpr double optStepBase = 8, optStepMul = 0.5;
constexpr double stickThreshold = 20;

constexpr double cellSize = 32;
constexpr double timeGamma = 1 - 1.0 / maxLookahead;
constexpr double dangerMultiplier = 4, passBonus = 2, passMultiplier = 0.1;
constexpr double dangerThreshold = 0.1;

constexpr double hockeyistFrict = 0.02, angularFrict = 0.0270190131;
constexpr double puckFrict = 0.001, puckBeta = -log(1 - puckFrict);
constexpr double wallBounce = 0.25, hockeyistBounce = 0.325, goalieFrict = 0.1;
constexpr double minDepth = 0.01, depthFactor = 0.8;

bool leftPlayer;
double rinkLeft, rinkRight, rinkTop, rinkBottom;
double rinkWidth, rinkHeight;  Vec2D rinkCenter;
int gridHalfWidth, gridHalfHeight, gridLine, gridHeight, gridCenter, gridSize;
double goalCenter, goalHalf, hockeyistRad, puckRad, goalieSpd, goalieRange;
double staminaBase, staminaGrowth, accelMin, accelMax, turnAngle;
double stickLength, effStickLength, stickSector, stickSectorTan, holdDist, timeGammaSwing;  int maxSwing;
double strikeBase, strikeGrowth, passBase, strikeBeta, passBeta, timeSpread;  Vec2D passSector;
double minChance, maxChance, knockdownChance, chanceDrop;
double strikeChance, pickChanceDelta[PUCK_STATUS_COUNT];

int globalTick = -1, goalieTime, activeCount = 0;


void initConsts(const model::Game& game, const model::World& world)
{
    rinkLeft = game.getRinkLeft();  rinkRight = game.getRinkRight();
    rinkTop = game.getRinkTop();  rinkBottom = game.getRinkBottom();
    rinkWidth = rinkRight - rinkLeft;  rinkHeight = rinkBottom - rinkTop;
    rinkCenter.x = (rinkLeft + rinkRight) / 2;
    rinkCenter.y = (rinkTop + rinkBottom) / 2;

    leftPlayer = (world.getMyPlayer().getNetBack() < rinkCenter.x);

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

    stickLength = game.getStickLength();  effStickLength = stickLength / sqrt(3.0);
    stickSector = game.getStickSector() / 2;  stickSectorTan = tan(stickSector);
    holdDist = game.getPuckBindingRange();  timeGammaSwing = pow(timeGamma, maxSwing);
    maxSwing = game.getMaxEffectiveSwingTicks();

    strikeBase   = game.getStruckPuckInitialSpeedFactor() * game.getStrikePowerBaseFactor();
    strikeGrowth = game.getStruckPuckInitialSpeedFactor() * game.getStrikePowerGrowthFactor();
    passBase     = game.getStruckPuckInitialSpeedFactor() * game.getPassPowerFactor();
    strikeBeta = 1 / (game.getStrikeAngleDeviation() * sqrt(2.0));
    passBeta   = 1 / (game.getPassAngleDeviation()   * sqrt(2.0));
    timeSpread = strikeBase / strikeBeta;
    passSector = sincos(game.getPassSector() / 2);

    minChance = game.getMinActionChance();
    maxChance = game.getMaxActionChance();
    knockdownChance = game.getKnockdownChanceFactor();
    chanceDrop = 1 / game.getStruckPuckInitialSpeedFactor();
    strikeChance = game.getStrikePuckBaseChance();
    pickChanceDelta[HOLD] = strikeChance - game.getTakePuckAwayBaseChance();
    pickChanceDelta[FLY]  = strikeChance - game.getPickUpPuckBaseChance();
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

struct StickSector
{
    Vec2D pos, dir, side;
    double back, forw, base, coeff;

    void set(const Vec2D &pos_, const Vec2D &dir_, double back_, const Vec2D &forw_, double angle)
    {
        pos = pos_;  dir = dir_;  side = sincos(angle + stickSector);  back = back_;
        double phi = angle + stickSector / 2, dr = forw_.y / sin(phi);  forw = forw_.x - dr * cos(phi);
        base = 0.5 * (stickLength + dr);  coeff = 0.5 / (stickLength + dr);
    }

    double distance(const Vec2D &pt) const
    {
        double dot = (pt - pos) * dir, cross = (pt - pos) % dir;
        return max((sqr(dot - forw) + sqr(cross)) * coeff - base, side.x * abs(cross) - side.y * (dot + back));
    }
};

struct HockeyistInfo
{
    long long id;  int index;  double stamina;
    double accelMin, accelMax, turnAngle, maxTurnSteps;
    int maxTurnCount;  Vec2D turnRot, endTurnRot;
    double strikeBase, strikeGrowth, passBase, strikeBeta, passBeta;
    double attack[PUCK_STATUS_COUNT], defence;
    double knockoutChance, knockdownChance;
    double followMoveCoeff, followRotCoeff;
    double staminaCoeff, staminaBias, defenceBias;

    bool havePuck;
    int knockdown, cooldown;
    vector<StickSector> stick;
    vector<Vec2D> rot;

    HockeyistInfo(const model::Hockeyist& hockeyist) : id(hockeyist.getId())
    {
    }

    void set(const model::Hockeyist& hockeyist, long long puckOwner, int index_)
    {
        index = index_;  stamina = hockeyist.getStamina();
        double mul = staminaBase + staminaGrowth * stamina;
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
        knockoutChance = toChance(::strikeChance - defence + 1);
        knockdownChance = toChance(::knockdownChance - end + 1);

        double timeCoeff = 0.01 / (1 - timeGamma);
        followMoveCoeff = timeCoeff * sqr(hockeyistFrict / accelMax);
        followRotCoeff = 0.02 * timeCoeff / sqr(turnAngle);
        staminaCoeff = 1;  staminaBias = 500 * staminaCoeff;
        defenceBias = 500 * defence;


        havePuck = (id == puckOwner);
        knockdown = hockeyist.getRemainingKnockdownTicks();
        cooldown = max(knockdown, hockeyist.getRemainingCooldownTicks());

        double angSpd = hockeyist.getAngularSpeed();
        rot.resize(maxLookahead + 1);  rot[0] = sincos(angSpd);
        for(int i = 1; i <= maxLookahead; i++)
        {
            angSpd -= angularFrict * angSpd;  rot[i] = sincos(angSpd);
        }

        UnitState state;  state.set(hockeyist);  Vec2D dir = sincos(hockeyist.getAngle());
        stick.resize(stickLookahead + 1);  stick[0].set(state.pos, dir, 0, Vec2D(0, 0), 0);
        double back = 0, backSpd = 0;  Vec2D forw(0, 0), forwSpd(0, 0), forwDir(1, 0);
        for(int i = 1; i <= stickLookahead; i++)
        {
            state.spd -= hockeyistFrict * state.spd;
            state.pos += state.spd;  dir = rotate(dir, rot[i]);
            if(i > knockdown)
            {
                backSpd -= accelMin;  backSpd -= hockeyistFrict * backSpd;  back += backSpd;
                forwSpd += accelMax * forwDir;  forwSpd -= hockeyistFrict * forwSpd;
                forw += forwSpd;  forwDir = rotate(forwDir, turnRot);
            }
            stick[i].set(state.pos, dir, back, forw, (i - knockdown) * turnAngle);
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

struct SubstituteInfo
{
    int index;  double stamina;

    SubstituteInfo(const model::Hockeyist& hockeyist) :
        index(hockeyist.getTeammateIndex()), stamina(hockeyist.getStamina())
    {
    }

    bool operator < (const SubstituteInfo &info) const
    {
        return stamina > info.stamina;
    }
};


struct Safety
{
    double val[MAX_ACTIVE];

    Safety(double val_ = 1)
    {
        for(int i = 0; i < activeCount; i++)val[i] = val_;
    }

    Safety &update(const Safety &safety)
    {
        for(int i = 0; i < activeCount; i++)val[i] = min(val[i], safety.val[i]);  return *this;
    }

    double update(int index, double val_)
    {
        return val[index] = min(val[index], val_);
    }

    double multiply(double mul) const
    {
        for(int i = 0; i < activeCount; i++)mul *= val[i];  return mul;
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

    int nextStep(int time)
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
        if(time <= goalieTime)for(int i = 0; i < 2; i++)
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
    Safety intact;  double chanceDrop;

    void set(const PuckState &state, int time)
    {
        pos = state.pos;  status = FLY;  intact = 1;
        chanceDrop = ::chanceDrop * state.spd.len();
    }

    void set(const HockeyistInfo &info, const HockeyistState &state, int time)
    {
        pos = puckPos(info, state, time);  status = HOLD;  intact = 1;
        chanceDrop = info.defence;
    }
};



namespace MoveFlag
{
    enum
    {
        FORW = 0, FIRST_LEFT  = 0, SECOND_LEFT  = 0,
        BACK = 1, FIRST_RIGHT = 2, SECOND_RIGHT = 4,
        ROTATE_ONLY = 8, STOP = 16, SUBSTITUTE = 32
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

struct EnemyMap
{
    vector<double> danger;
    vector<int> stick[MAX_ACTIVE];

    void init();

    void reset()
    {
        for(int i = 0; i < activeCount; i++)
        {
            stick[i].resize(gridSize);
            for(auto &flag : stick[i])flag = maxLookahead;
        }
    }

    void addMapPoint(const HockeyistInfo &info, HockeyistState &state, int time, int flags, double turnTime)
    {
        auto &cell = stick[info.index][gridPos(state.pos + effStickLength * state.dir)];
        cell = min<int>(cell, max(time, info.cooldown));
    }

    void closePath(const HockeyistInfo &info, HockeyistState &state, int flags, double turnTime)
    {
    }

    void fillMap(const HockeyistInfo &info, const HockeyistState &state)
    {
        Mapper<EnemyMap, HockeyistState, &EnemyMap::addMapPoint, &EnemyMap::closePath>(*this, info).fillMap(state);
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
        for(int i = 0; i < activeCount; i++)filterField(stick[i]);
    }

    double getDanger(const Vec2D &pos)
    {
        Vec2D offs = (pos - rinkCenter) / cellSize + Vec2D(gridHalfWidth, gridHalfHeight);
        int x = int(offs.x), y = int(offs.y), cell = y * gridLine + x;  offs -= Vec2D(x, y);
        double val1 = danger[cell] + (danger[cell + 1] - danger[cell]) * offs.x;  cell += gridLine;
        double val2 = danger[cell] + (danger[cell + 1] - danger[cell]) * offs.x;
        return val1 + (val2 - val1) * offs.y;
    }

    void updateSafety(Safety &res, const Vec2D &puck, double spd, int time);
    void updateSafety(Safety &res, const HockeyistInfo &info, const HockeyistState &state, int time, bool havePuck);

    int enemyTime(int cell) const
    {
        int res = maxLookahead;
        for(int i = 0; i < activeCount; i++)res = min(res, stick[i][cell]);
        return res;
    }
};


struct TargetPlace
{
    Vec2D pos, view;
};

struct PassTarget;
struct EnemyInfo;
struct AllyInfo;

EnemyMap enemyMap;
double currentDanger;
vector<SubstituteInfo> substitutes;
TargetPlace attackPlace, defencePlace;
vector<PassTarget> targets, nextTargets;
PathPoint puckPath[maxLookahead + 1];
int puckPathLen, goalFlag;
EnemyInfo *enemyPuck;
AllyInfo *allyPuck;


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

    Safety intact() const
    {
        return 1;
    }
};

template<> struct StateScore<true> : public HockeyistState
{
    double timeFactor;
    Safety safety;

    StateScore(const HockeyistState &state) : HockeyistState(state), timeFactor(1), safety(1)
    {
    }

    void nextStep(const HockeyistInfo &info, double accel, const Vec2D &turn, int time)
    {
        enemyMap.updateSafety(safety, info, *this, time - 1, !puckPathLen);
        HockeyistState::nextStep(info, accel, turn, time);  timeFactor *= timeGamma;
    }

    const Safety &intact() const
    {
        return safety;
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
    Vec2D start, bonus, dir;
    double rad, power, dist, offs, end, duration;

    GoalHelper(const Vec2D &pos, const Vec2D &spd, double power_, bool right) :
        start(pos), bonus(spd), rad(hockeyistRad + puckRad), power(power_)
    {
        if(right)
        {
            start.x = 2 * rinkCenter.x - start.x;  bonus.x = -bonus.x;
        }
        if(pos.y < goalCenter)
        {
            start.y = 2 * goalCenter - start.y;  bonus.y = -bonus.y;
        }
    }

    bool init()
    {
        start -= Vec2D(rinkLeft, goalCenter - goalHalf);
        dist = start.x - hockeyistRad;  if(!(dist > rad))return false;
        offs = max(0.0, start.y - goalHalf - goalieRange);
        dir = start = normalize(start);  return true;
    }

    void withoutGoalie()
    {
        start -= Vec2D(rinkLeft, goalCenter - goalHalf);
        dir = normalize(Vec2D(start.x, start.y - goalHalf - goalieRange));
        double len = start.len();  start /= len;

        double end = len / (power - bonus * start);
        duration = -log(1 - puckBeta * end) / puckBeta;
    }

    template<bool first> bool iterate()
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
        if(first)duration = end;  return true;
    }

    Sector result(const Vec2D &pos, bool right) const
    {
        double dot = dir * start, cross = dir % start;
        double norm = 1 / sqrt(2 + 2 * dot);  Vec2D res = dir + start;
        if(!right)res.x = -res.x;  if(!(pos.y < goalCenter))res.y = -res.y;
        return Sector(res * norm, cross * norm, duration);
    }
};

Sector estimateGoalAngle(const Vec2D &pos, const Vec2D &spd, double power, bool right)
{
    GoalHelper helper(pos, spd, power, right);
    if(goalieTime > 0)
    {
        if(!helper.init())return Sector();
        if(!helper.iterate<true>())return Sector();
        for(int i = 0; i < 3; i++)helper.iterate<false>();
    }
    else helper.withoutGoalie();  return helper.result(pos, right);
}


void interception(Safety &res, const Vec2D &pos, const Vec2D &spd, int time, int duration)
{
    SimplePuckState puck(pos, spd);
    for(int i = 0; i < duration; i++)
    {
        enemyMap.updateSafety(res, puck.pos, puck.spdLen, time + i);  if(!puck.nextStep())break;
    }
}

inline bool inStrikeSector(const HockeyistState &state, const Vec2D &puck)
{
    Vec2D delta = puck - state.pos;
    if(!(delta.sqr() < sqr(stickLength)))return false;
    double dot = delta * state.dir, cross = delta % state.dir;
    return abs(cross) < dot * stickSectorTan;
}

inline double probFunc(double val)
{
    constexpr double eps = 1e-3;
    double val1 = erf(val), val2 = val / (1 + abs(val));
    return (val1 + eps * (val2 - val1)) / 2;
}

inline double probFactor(double beta, double center, double span)
{
    return probFunc(beta * (center + span)) - probFunc(beta * (center - span));
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

    void tryPickup(double strike, int time, double val, const Safety &safety)
    {
        int deltaTime = enemyMap.enemyTime(gridPos(puckPath[time].pos)) - time;
        val *= safety.multiply(currentDanger + passMultiplier * (1 + double(deltaTime) / maxLookahead));
        val *= toChance(strike - pickChanceDelta[puckPath[time].status]);  if(!(val > score))return;

        strikeTime = time;  swingTime = -1;  targetIndex = -1;  score = val;
    }

    template<bool ally> void evaluateStrike(const HockeyistState &state, const Vec2D &puck,
        double power, double beta, int swing, int time, double mul, const Safety &safety, double bonus = 0)
    {
        double val = 0;  Vec2D delta(1, 0);
        Sector res = estimateGoalAngle(puck, state.spd, power, ally == leftPlayer);
        if(res.span > 0)
        {
            delta = rotate(res.dir, conj(state.dir));  Vec2D offs(delta.x, abs(delta.y));
            if(swing < 0)offs = (offs.x > passSector.x ? Vec2D(1, 0) : rotate(offs, conj(passSector)));
            if(offs.x > 0)
            {
                val = mul * probFactor(beta, offs.y, res.span);
                if(ally)
                {
                    if(!(safety.multiply(val) + bonus > score))return;  Safety intact = safety;
                    interception(intact, puck, (power + state.spd * res.dir) * res.dir, time + swing, lround(res.time));
                    val = intact.multiply(val);
                }
                else val = safety.multiply(val);
            }
        }
        if(!(val + bonus > score))return;

        strikeTime = time;  swingTime = swing;  targetIndex = -1;
        score = val + bonus;  passPower = 1;  passDir = delta;  strikePos = puck;
    }

    void evaluatePass(const HockeyistState &state, const Vec2D &puck, double power,
        double x, double y, double spd, int catchTime, int time, double val, const Safety &safety, int target)
    {
        Vec2D delta = Vec2D(x, y) - puck;
        double len = delta.len();  delta /= len;
        Vec2D offs = rotate(delta, conj(state.dir));
        if(!(offs.x > passSector.x))return;

        double puckSpd = puckBeta * len / (1 - exp(-puckBeta * (catchTime - time)));
        double relSpd = state.spd * delta;  if(!(puckSpd < power + relSpd))return;
        val *= maxChance - max(0.0, puckSpd - spd - puckBeta * len) * chanceDrop;
        if(!(safety.multiply(val) > score))return;

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
        Safety intact = safety;
        interception(intact, puck, puckSpd * delta, time, lround(duration));
        if(!((val = intact.multiply(val)) > score))return;

        strikeTime = time;  swingTime = -1;  targetIndex = target;
        score = val;  passPower = (puckSpd - relSpd) / power;  passDir = offs;
    }

    void tryPass(const HockeyistState &state, const Vec2D &puck, double power, int time, double mul, const Safety &safety);

    template<bool ally> void tryStrike(const HockeyistInfo &info, const StateScore<ally> &state, int time)
    {
        double mul = state.timeFactor;  Vec2D puck;  Safety safety = state.intact();
        if(time + maxSwing <= maxLookahead && safety.multiply(mul * timeGammaSwing) > score)do
        {
            StateScore<ally> cur = state;
            for(int i = 1; i <= maxSwing; i++)cur.nextStep(info, 0, Vec2D(1, 0), time + i);

            double val = cur.timeFactor;  Safety safety = cur.intact();
            if(!puckPathLen)puck = puckPos(info, cur, time + maxSwing);
            else
            {
                if(time >= puckPathLen || !inStrikeSector(cur, puck = puckPath[time + maxSwing].pos))break;
                double strike = info.attack[puckPath[time + maxSwing].status] - puckPath[time + maxSwing].chanceDrop;
                safety.update(puckPath[time + maxSwing].intact);  val *= toChance(strike);
            }
            if(safety.multiply(val) > score)evaluateStrike<ally>(cur, puck,
                info.strikeBase + info.strikeGrowth * maxSwing, info.strikeBeta, maxSwing, time, val, safety);
        }
        while(false);

        if(!puckPathLen)
        {
            puck = puckPos(info, state, time);
            if(ally)tryPass(state, puck, info.passBase, time, mul, safety);
        }
        else
        {
            if(time >= puckPathLen || !inStrikeSector(state, puck = puckPath[time].pos))return;
            double strike = info.attack[puckPath[time].status] - puckPath[time].chanceDrop;
            safety.update(puckPath[time].intact);  if(ally)tryPickup(strike, time, mul, safety);
            mul *= toChance(strike);
        }
        if(info.havePuck)
        {
            if(!(safety.multiply(mul) > score))return;
            evaluateStrike<ally>(state, puck, info.passBase, info.passBeta, -1, time, mul, safety);
        }
        else
        {
            double bonus = 0;
            if(ally && puckPathLen)
            {
                double autoGoal = 0;
                Sector res = estimateGoalAngle(puck, state.spd, info.strikeBase, !leftPlayer);
                if(res.span > 0)
                {
                    Vec2D offs = rotate(res.dir, conj(state.dir));
                    if(offs.x > 0)autoGoal = probFactor(info.strikeBeta, abs(offs.y), res.span);
                }
                bonus = safety.multiply(mul * currentDanger) * (1 - autoGoal);
            }
            if(!(safety.multiply(mul) + bonus > score))return;
            evaluateStrike<ally>(state, puck, info.strikeBase, info.strikeBeta, 0, time, mul, safety, bonus);
        }
    }

    template<bool ally> void trySwing(const HockeyistInfo &info, const StateScore<ally> &state, int swing, int time)
    {
        double val = state.timeFactor;  Vec2D puck;  Safety safety = state.intact();
        if(!puckPathLen)puck = puckPos(info, state, time);
        else
        {
            if(!inStrikeSector(state, puck = puckPath[time].pos))return;
            double strike = info.attack[puckPath[time].status] - puckPath[time].chanceDrop;
            safety.update(puckPath[time].intact);  val *= toChance(strike);
        }
        if(safety.multiply(val) > score)evaluateStrike<ally>(state, puck,
            info.strikeBase + info.strikeGrowth * swing, info.strikeBeta, time, 0, val, safety);
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
            if(plan.flags & MoveFlag::ROTATE_ONLY)accBase = accDelta = 0;
            else
            {
                accBase = info.accelMin;  accDelta = info.accelMax;
                if(plan.flags & MoveFlag::BACK)swap(accBase, accDelta);  accDelta -= accBase;
            }
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

    void followPlace(const HockeyistInfo &info, const HockeyistState &state, const Vec2D &target)
    {
        Vec2D delta = rotate(target - state.pos, conj(state.dir));  double turn = atan2(delta.y, delta.x);
        *static_cast<StrikeInfo *>(this) = StrikeInfo(maxLookahead, info.staminaBias);
        flags = (turn > 0 ? MoveFlag::FIRST_LEFT : MoveFlag::FIRST_RIGHT) | MoveFlag::STOP;
        firstTurnEnd = abs(turn) / info.turnAngle;  secondTurnStart = maxLookahead;
    }

    static double placeScore(const HockeyistInfo &info, const Vec2D &pos, const TargetPlace &place, double dirFactor)
    {
        return info.followRotCoeff * dirFactor - info.followMoveCoeff * (pos - place.pos).sqr();
    }

    void tryPlace(const HockeyistInfo &info, const Vec2D &pos, const Vec2D &end, const Vec2D &dir,
        const TargetPlace &place, int flags_, double turnTime, int endTime, double val)
    {
        Vec2D view = normalize(place.view - pos);
        val += placeScore(info, end, place, dir * view) - endTime;  if(!(val > score))return;

        *static_cast<StrikeInfo *>(this) = StrikeInfo(endTime, val);
        flags = flags_ | MoveFlag::STOP;  firstTurnEnd = turnTime;  secondTurnStart = endTime;
    }

    static double substituteScore(const HockeyistInfo &info, const Vec2D &pos)
    {
        double dx = max(0.0, (leftPlayer ? pos.x - rinkCenter.x : rinkCenter.x - pos.x) + 2 * hockeyistRad);
        double dy = max(0.0, pos.y - rinkTop), val = info.followMoveCoeff * (dx * dx + dy * dy);
        return info.staminaCoeff * (substitutes[0].stamina - info.stamina) - val;
    }

    void trySubstitute(const HockeyistInfo &info, const Vec2D &pos, int index, int flags_, double turnTime, int endTime)
    {
        double val = substituteScore(info, pos) - endTime;  if(!(val > score))return;

        *static_cast<StrikeInfo *>(this) = StrikeInfo(endTime, val);  targetIndex = index;
        flags = flags_ | MoveFlag::STOP | MoveFlag::SUBSTITUTE;  firstTurnEnd = turnTime;  secondTurnStart = endTime;
    }

    void updatePlace(const HockeyistInfo &info, const HockeyistState &state, const TargetPlace &place, double val)
    {
        HockeyistState cur = state;  Helper helper(info, *this);
        for(int i = 0; i < strikeTime; i++)cur.nextStep(info, helper.accel(i), helper.turn(i), i + 1);
        Vec2D end = cur.pos + cur.spd / hockeyistFrict;
        if(flags & MoveFlag::SUBSTITUTE)
        {
            score = substituteScore(info, end) - strikeTime;  return;
        }

        Vec2D view = normalize(place.view - cur.pos);
        if((flags & MoveFlag::ROTATE_ONLY) || !strikeTime)
        {
            Vec2D delta = rotate(view, conj(state.dir));  double turn = atan2(delta.y, delta.x);
            flags = (turn > 0 ? MoveFlag::FIRST_LEFT : MoveFlag::FIRST_RIGHT) | MoveFlag::ROTATE_ONLY | MoveFlag::STOP;
            firstTurnEnd = abs(turn) / info.turnAngle;  secondTurnStart = strikeTime = int(firstTurnEnd + 1);
            score = val + placeScore(info, end, place, 1);
        }
        else score = val + placeScore(info, end, place, cur.dir * view) - strikeTime;
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

    template<bool ally> void evaluateSwing(const HockeyistInfo &info, const HockeyistState &state, int swing)
    {
        reset();  StateScore<ally> cur = state;  trySwing<ally>(info, cur, swing, 0);
        for(int i = swing + 1; i <= maxSwing; i++)
        {
            cur.nextStep(info, 0, Vec2D(1, 0), i);  trySwing<ally>(info, cur, i, i - swing);
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

    void execute(const HockeyistInfo &info, model::Move &move, int swinging) const
    {
        if(swinging)
        {
            move.setAction(strikeTime ? model::CANCEL_STRIKE : (swingTime ? model::SWING : model::STRIKE));
            move.setSpeedUp(0);  move.setTurn(0);  return;
        }

        if(strikeTime)
        {
            Helper helper(info, *this);  double accel = helper.accel(0);
            move.setSpeedUp(accel > 0 ? accel / info.accelMax : -accel / info.accelMin);
            Vec2D rot = helper.turn(0);  move.setTurn(atan2(rot.y, rot.x));
            if(flags & MoveFlag::SUBSTITUTE)
            {
                move.setTeammateIndex(targetIndex);  move.setAction(model::SUBSTITUTE);
            }
            else move.setAction(model::NONE);  return;
        }
        move.setSpeedUp(0);  move.setTurn(0);
        if(flags & MoveFlag::STOP)move.setAction(model::NONE);
        else if(swingTime > 0)move.setAction(model::SWING);
        else if(!swingTime)move.setAction(model::STRIKE);
        else if(puckPathLen)move.setAction(model::TAKE_PUCK);
        else
        {
            move.setAction(model::PASS);  move.setPassPower(passPower);
            move.setPassAngle(atan2(passDir.y, passDir.x));
        }
    }

    void skipTurn()
    {
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
};

struct PassTarget : public Vec2D, public MovePlan
{
    AllyInfo *receiver;  double spd;

    PassTarget() = default;

    PassTarget(AllyInfo *rcv, const HockeyistInfo &info, const MovePlan &plan) : Vec2D(plan.strikePos), MovePlan(plan), receiver(rcv)
    {
        spd = (info.attack[FLY] - pickChanceDelta[FLY] - maxChance) / chanceDrop;  score *= passBonus;
    }

    PassTarget(AllyInfo *rcv, const HockeyistInfo &info, const Vec2D &puck, int flags, double turnTime, int strikeTime, double score) :
        Vec2D(puck), MovePlan(flags, turnTime, strikeTime, score), receiver(rcv)
    {
        spd = (info.attack[FLY] - pickChanceDelta[FLY] - maxChance) / chanceDrop;
    }

    bool operator < (const PassTarget &cmp) const
    {
        return score < cmp.score;
    }
};

void StrikeInfo::tryPass(const HockeyistState &state, const Vec2D &puck, double power, int time, double mul, const Safety &safety)
{
    double y1 = rinkTop + puckRad, y2 = rinkBottom - puckRad;
    for(size_t i = 0; i < targets.size(); i++)
    {
        int hit = targets[i].strikeTime + targets[i].swingTime;  if(hit <= time)continue;
        double val = mul * targets[i].score;  if(!(safety.multiply(val) > score))continue;

        evaluatePass(state, puck, power, targets[i].x, targets[i].y, targets[i].spd, hit, time, val, safety, i);
        evaluatePass(state, puck, power, targets[i].x, y1 - (targets[i].y - y1) / wallBounce, targets[i].spd, hit, time, val, safety, i);
        evaluatePass(state, puck, power, targets[i].x, y2 - (targets[i].y - y2) / wallBounce, targets[i].spd, hit, time, val, safety, i);
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
        if(state.score > 0)moves.emplace_back(state, flags, turnTime);
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
        plan.evaluate<ally>(info, state);  moves.push_back(plan);

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

    EnemyInfo(const model::Hockeyist& hockeyist, long long puckOwner, int index) : HockeyistInfo(hockeyist)
    {
        set(hockeyist, puckOwner, index);
    }

    void set(const model::Hockeyist& hockeyist, long long puckOwner, int index)
    {
        HockeyistInfo::set(hockeyist, puckOwner, index);  HockeyistState::set(hockeyist);
        reset();  swinging = hockeyist.getSwingTicks();
    }

    void addPoint(const HockeyistInfo &info, State &state, int time, int flags, double turnTime)
    {
        Optimizer<false>::addPoint(info, state, time, flags, turnTime);
        enemyMap.addMapPoint(info, state, time, flags, turnTime);
    }

    void closePath(const HockeyistInfo &info, State &state, int flags, double turnTime)
    {
        Optimizer<false>::closePath(info, state, flags, turnTime);
        enemyMap.closePath(info, state, flags, turnTime);
    }

    void process()
    {
        if(!puckPathLen)enemyMap.fillMap(*this, *this);
        else if(knockdown)plan = MovePlan();
        else if(swinging)plan.evaluateSwing<false>(*this, *this, swinging);
        else
        {
            Mapper<EnemyInfo, State, &EnemyInfo::addPoint, &EnemyInfo::closePath>(*this, *this).fillMap(*this);
            findBestMove(*this, *this, plan);
        }
    }

    const MovePlan &choosePuckMove()
    {
        if(!swinging)
        {
            Optimizer<false>::fillMap(*this, *this);  findBestMove(*this, *this, plan);
        }
        else plan.evaluateSwing<false>(*this, *this, swinging);  return plan;
    }

    bool execute()
    {
        if(havePuck)puckPath[0].set(*this, *this, 0);  int time = 1;
        HockeyistState state = *this;  MovePlan::Helper helper(*this, plan);
        for(; time <= plan.strikeTime; time++)
        {
            state.nextStep(*this, helper.accel(time - 1), helper.turn(time - 1), time);
            if(havePuck)puckPath[time].set(*this, state, time);
        }
        for(; time <= plan.strikeTime + plan.swingTime; time++)
        {
            state.nextStep(*this, 0, Vec2D(1, 0), time);
            if(havePuck)puckPath[time].set(*this, state, time);
        }

        Vec2D dir = state.dir;  double power;
        if(plan.swingTime < 0)
        {
            dir = rotate(dir, plan.passDir);  power = passBase * plan.passPower;
        }
        else power = strikeBase + strikeGrowth * plan.swingTime;
        power += state.spd * dir;

        PuckState puck(puckPath[time - 1].pos, power * dir);  int flag = 0;
        if(!havePuck && !inStrikeSector(state, puck.pos))return false;
        for(; time <= maxLookahead; time++)
        {
            if((flag = puck.nextStep(time)))break;  puckPath[time].set(puck, time);
        }
        puckPathLen = time;  goalFlag = leftPlayer ? flag : -flag;
        plan.skipTurn();  return true;
    }
};

struct AllyInfo : public HockeyistInfo, public HockeyistState, public Optimizer<true>
{
    enum PlanType
    {
        ATTACK, DEFENCE, MAIN, PLAN_TYPE_COUNT
    };

    PlanType activePlan;
    MovePlan plan[PLAN_TYPE_COUNT];
    vector<pair<double, int>> stick;
    int swinging;

    AllyInfo(const model::Hockeyist& hockeyist, long long puckOwner, int index) : HockeyistInfo(hockeyist)
    {
        set(hockeyist, puckOwner, index);
    }

    void set(const model::Hockeyist& hockeyist, long long puckOwner, int index)
    {
#ifdef CHECK_PREDICTION
        Vec2D errPos = pos, errSpd = spd, errDir = dir;
#endif

        HockeyistInfo::set(hockeyist, puckOwner, index);  HockeyistState::set(hockeyist);
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

        Vec2D end = state.pos + state.spd / hockeyistFrict;
        if(substitutes.size() && substitutes[0].stamina > stamina)
            plan[ATTACK].trySubstitute(info, end, substitutes[0].index, flags, turnTime, time);
        else if(!enemyPuck)plan[ATTACK].tryPlace(info, state.pos, end, state.dir, attackPlace, flags, turnTime, time, staminaBias);
        plan[DEFENCE].tryPlace(info, state.pos, end, state.dir, defencePlace, flags, turnTime, time, defenceBias);

        Vec2D puck = state.pos + effStickLength * state.dir;
        if(!validPuckPos(puck))return;  int cell = gridPos(puck);
        double score = state.safety.multiply(state.timeFactor * passMultiplier);
        score *= 1 - sqr(2 * (state.pos.x - rinkCenter.x) / rinkWidth);
        score *= 1 - sqr(2 * (state.pos.y - rinkCenter.y) / rinkHeight);
        int deltaTime = enemyMap.enemyTime(cell) - max(time, cooldown);
        score *= 0.5 + deltaTime * (0.5 / maxLookahead);
        if(!(score > stick[cell].first))return;

        stick[cell].first = score;
        PassTarget target(this, info, puck, flags, turnTime, time, score);
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

    bool tryKnockout()
    {
        if(!enemyPuck || cooldown || !inStrikeSector(*this, puckPath[0].pos))return false;

        double autoGoal = 0;
        Sector res = estimateGoalAngle(puckPath[0].pos, spd, strikeBase, !leftPlayer);
        if(res.span > 0)
        {
            Vec2D offs = rotate(res.dir, conj(dir));
            if(offs.x > 0)autoGoal = probFactor(strikeBeta, abs(offs.y), res.span);
        }
        if(autoGoal > 1e-3)return false;

        plan[MAIN] = MovePlan();  plan[MAIN].flags = 0;  activePlan = MAIN;  return true;
    }

    const MovePlan &chooseInterceptMove()
    {
        if(knockdown)return plan[activePlan = ATTACK] = plan[DEFENCE] = MovePlan();
        if(swinging)
        {
            plan[ATTACK] = plan[DEFENCE] = MovePlan();
            plan[MAIN].evaluateSwing<true>(*this, *this, swinging);  return plan[activePlan = MAIN];
        }

        activePlan = ATTACK;
        if(enemyPuck)plan[ATTACK].followPlace(*this, *this, puckPath[0].pos + enemyPuck->spd * (0.5 / hockeyistFrict));
        else plan[ATTACK].updatePlace(*this, *this, attackPlace, staminaBias);
        plan[DEFENCE].updatePlace(*this, *this, defencePlace, defenceBias);

        stick.resize(gridSize);  for(auto &flag : stick)flag = {0, -1};
        Mapper<AllyInfo, State, &AllyInfo::addPoint, &AllyInfo::closePath>(*this, *this).fillMap(*this);
        findBestMove(*this, *this, plan[MAIN]);  return plan[MAIN];
    }

    const MovePlan &choosePuckMove()
    {
        if(!swinging)
        {
            Optimizer<true>::fillMap(*this, *this);  findBestMove(*this, *this, plan[MAIN]);
        }
        else plan[MAIN].evaluateSwing<true>(*this, *this, swinging);  return plan[MAIN];
    }

    void substitute(const SubstituteInfo &sub)
    {
        activePlan = ATTACK;  double x;
        if(leftPlayer)x = min(pos.x, rinkCenter.x - 2 * hockeyistRad);
        else          x = max(pos.x, rinkCenter.x + 2 * hockeyistRad);
        plan[ATTACK].followPlace(*this, *this, Vec2D(x, rinkTop));
        plan[ATTACK].flags |= MoveFlag::SUBSTITUTE;
        plan[ATTACK].targetIndex = sub.index;
    }

    void doNothing()
    {
        plan[activePlan = ATTACK] = plan[DEFENCE] = MovePlan();
    }

    void execute(model::Move &move)
    {
        plan[activePlan].execute(*this, move, swinging);
        for(int i = 0; i < PLAN_TYPE_COUNT; i++)plan[i].skipTurn();

#ifdef PRINT_LOG
        switch(move.getAction())
        {
        case model::TAKE_PUCK:      cout << ">>>>>>>>>>>>> TAKE_PUCK"     << endl;  break;
        case model::SWING:          cout << ">>>>>>>>>>>>> SWING"         << endl;  break;
        case model::STRIKE:         cout << ">>>>>>>>>>>>> STRIKE"        << endl;  break;
        case model::CANCEL_STRIKE:  cout << ">>>>>>>>>>>>> CANCEL_STRIKE" << endl;  break;
        case model::PASS:           cout << ">>>>>>>>>>>>> PASS"          << endl;  break;
        //case model::SUBSTITUTE:     cout << ">>>>>>>>>>>>> SUBSTITUTE"    << endl;  break;
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
            list.emplace_back(hockeyist, puckOwner, ref.size());  ref.push_back(&*list.rbegin());  break;
        }
        else if(iter->id == hockeyist.getId())
        {
            iter->set(hockeyist, puckOwner, ref.size());  ref.push_back(&*iter);  break;
        }
}

vector<EnemyInfo> enemyInfo;
vector<EnemyInfo *> enemies;
vector<AllyInfo> allyInfo;
vector<AllyInfo *> allies;


inline double stickFactor(const StickSector &sector, const Vec2D &pt)
{
    return stepValue(1 - sector.distance(pt) / stickThreshold);
}

void EnemyMap::init()
{
    danger.resize(gridSize);  int cell = 0;
    double power = strikeBase + accelMax / hockeyistFrict;
    for(int y = -gridHalfHeight; y <= gridHalfHeight; y++)
        for(int x = -gridHalfWidth; x <= gridHalfWidth; x++, cell++)
        {
            Sector res = estimateGoalAngle(rinkCenter + Vec2D(x, y) * cellSize, Vec2D(0, 0), power, !leftPlayer);
            double border = (rinkHeight / 2 + hockeyistRad - abs(y) * cellSize) / abs(res.dir.y * stickLength);
            danger[cell] = stepValue(border) * probFactor(strikeBeta, 0, res.span);
        }
    assert(cell == gridSize);
}

void EnemyMap::updateSafety(Safety &res, const Vec2D &puck, double spd, int time)
{
    double val = (1 - min(maxChance, 1 + strikeChance - chanceDrop * spd)) * (1 - getDanger(puck));
    if(time > stickLookahead)
    {
        int cell = gridPos(puck);
        for(int i = 0; i < activeCount; i++)if(stick[i][cell] <= time)res.update(i, val);
    }
    else for(auto &enemy : enemies)if(time >= enemy->cooldown)
        res.update(enemy->index, 1 - (1 - val) * stickFactor(enemy->stick[time], puck));
}

void EnemyMap::updateSafety(Safety &res, const HockeyistInfo &info, const HockeyistState &state, int time, bool havePuck)
{
    if(checkGoalie(state.pos, hockeyistRad))
    {
        res = 0;  return;
    }
    if(!havePuck)return;

    Vec2D puck = puckPos(info, state, time);
    if(!validPuckPos(puck))
    {
        res = 0;  return;
    }

    double val1 = 1 - info.knockdownChance;
    double val2 = (1 - info.knockoutChance) * (1 - getDanger(puck));
    if(time > stickLookahead)
    {
        int cell1 = gridPos(state.pos), cell2 = gridPos(puck);
        for(int i = 0; i < activeCount; i++)
        {
            if(stick[i][cell1] <= time)res.update(i, val1);
            if(stick[i][cell2] <= time)res.update(i, val2);
        }
    }
    else for(auto &enemy : enemies)if(time >= enemy->cooldown)
    {
        res.update(enemy->index, 1 - (1 - val1) * stickFactor(enemy->stick[time], state.pos));
        res.update(enemy->index, 1 - (1 - val2) * stickFactor(enemy->stick[time], puck));
    }
}


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

void drawField(int target)
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
            printLarge(cell == target ? 10 : gradient(enemyMap.danger[cell]));
        cout << "\x1B[0m\n";
    }
    cout.flush();
}


AllyInfo *chooseHockeyist(AllyInfo::PlanType action, double score = -numeric_limits<double>::infinity())
{
    AllyInfo *best = nullptr;
    for(auto ally : allies)if(!ally->activePlan && ally->plan[action].score > score)
    {
        score = ally->plan[action].score;  best = ally;
    }
#ifdef PRINT_LOG
    if(action == AllyInfo::MAIN && best)best->plan[action].print("Follower ", best->id);
#endif
    if(best)best->activePlan = action;  return best;
}

void createPlan(PuckState &puck, bool freePuck)
{
    if(freePuck)
    {
        puckPath[0].set(puck, 0);  int time = 1, flag = 0;
        for(; time <= maxLookahead; time++)
        {
            if((flag = puck.nextStep(time)))break;  puckPath[time].set(puck, time);
        }
        puckPathLen = time;  goalFlag = leftPlayer ? flag : -flag;
    }

    enemyMap.reset();  enemyPuck = nullptr;  Vec2D avg(0, 0);
    for(auto enemy : enemies)
    {
        avg += enemy->pos + enemy->spd / hockeyistFrict;
        if(enemy->havePuck)enemyPuck = enemy;
        else enemy->process();
    }
    enemyMap.filter();  avg /= enemies.size();

    double offs = 1.5 * goalHalf;
    attackPlace.pos = rinkCenter + Vec2D(0, avg.y > rinkCenter.y ? -offs : offs);
    attackPlace.view = attackPlace.pos + Vec2D(leftPlayer ? effStickLength : -effStickLength, 0);

    offs = rinkWidth / 2 - (goalieTime > 0 ? goalHalf : hockeyistRad);
    defencePlace.pos = rinkCenter + Vec2D(leftPlayer ? -offs : offs, 0);
    defencePlace.view = defencePlace.pos + effStickLength * normalize(puckPath[0].pos - defencePlace.pos);

    if(goalFlag < 0)currentDanger = dangerMultiplier;
    else if(freePuck)
    {
        EnemyInfo *best = nullptr;  double score = dangerThreshold;
        for(auto enemy : enemies)if(enemy != enemyPuck && enemy->plan.score > score)
        {
            score = enemy->plan.score;  best = enemy;
        }
        if(!best || !best->execute())
        {
            Safety intact;
            for(int time = 0; time < puckPathLen - 1; time++)
            {
                enemyMap.updateSafety(intact, puckPath[time].pos, puckPath[time].chanceDrop / chanceDrop, time);
                puckPath[time + 1].intact = intact;
            }
        }
        else currentDanger = dangerMultiplier * best->plan.score;
    }
    else if(enemyPuck)
    {
        enemyPuck->choosePuckMove();  enemyPuck->execute();
        currentDanger = dangerMultiplier * enemyPuck->plan.score;
    }
    else currentDanger = 0;

    targets.clear();  allyPuck = nullptr;
    for(auto ally : allies)
        if(ally->havePuck)
        {
            allyPuck = ally;  ally->activePlan = AllyInfo::MAIN;
        }
        else
        {
            const MovePlan &plan = ally->chooseInterceptMove();
            nextTargets.emplace_back(ally, *ally, plan);
        }
    chooseBest(nextTargets, targetCount);  swap(targets, nextTargets);

    int tg = -1;
    if(allyPuck && (tg = allyPuck->choosePuckMove().targetIndex) >= 0)
    {
        targets[tg].receiver->plan[AllyInfo::MAIN] = targets[tg];
        targets[tg].receiver->activePlan = AllyInfo::MAIN;
#ifdef PRINT_LOG
        targets[tg].print("Target ", targets[tg].receiver->id);
#endif
    }
#ifdef PRINT_LOG
    if(allyPuck)allyPuck->plan[AllyInfo::MAIN].print("Holder ", allyPuck->id);
#endif

    if(enemyPuck)
    {
        chooseHockeyist(AllyInfo::DEFENCE);  chooseHockeyist(AllyInfo::MAIN, dangerThreshold);
    }
    else
    {
        if(!allyPuck && goalFlag <= 0)chooseHockeyist(AllyInfo::MAIN, 0);
        chooseHockeyist(AllyInfo::DEFENCE);
    }

    for(auto ally : allies)if(ally->activePlan < AllyInfo::MAIN)ally->tryKnockout();
}

void MyStrategy::move(const model::Hockeyist& self, const model::World& world, const model::Game& game, model::Move& move)
{
    if(globalTick != world.getTick())
    {
        int totalGoals = 0;
        for(auto &player : world.getPlayers())totalGoals += player.getGoalCount();
        goalieTime = (totalGoals ? numeric_limits<int>::max() : world.getTickCount() - globalTick);

        if(!(globalTick = world.getTick()))
        {
            initConsts(game, world);  srand(game.getRandomSeed());
            enemyMap.init();  substitutes.reserve(MAX_ACTIVE);
            enemyInfo.reserve(6);  enemies.reserve(MAX_ACTIVE);
            allyInfo.reserve(6);  allies.reserve(MAX_ACTIVE);

            //drawField(-1);
        }

        long long owner = world.getPuck().getOwnerHockeyistId();
        PuckState puck;  puck.set(world.getPuck());
        puckPathLen = 0;  goalFlag = 0;

        enemies.clear();  allies.clear();  substitutes.clear();
        for(auto &hockeyist : world.getHockeyists())
        {
            if(hockeyist.getType() == model::GOALIE)
                puck.goalie[hockeyist.isTeammate() ? 0 : 1] = Vec2D(hockeyist.getX(), hockeyist.getY());
            else if(hockeyist.getState() != model::RESTING)
            {
                if(hockeyist.isTeammate())
                     synchronize(allyInfo,  allies,  hockeyist, owner);
                else synchronize(enemyInfo, enemies, hockeyist, owner);
            }
            else if(hockeyist.isTeammate())substitutes.emplace_back(hockeyist);
        }
        assert(allies.size() <= MAX_ACTIVE && enemies.size() <= MAX_ACTIVE);
        activeCount = enemies.size();

        sort(allies.begin(), allies.end(),
            [] (AllyInfo *ally1, AllyInfo *ally2) { return ally1->stamina < ally2->stamina; } );
        sort(substitutes.begin(), substitutes.end());

        const auto &player = world.getMyPlayer();
        if(player.isJustMissedGoal() || player.isJustScoredGoal())
        {
            for(size_t i = 0; i < allies.size(); i++)
                if(i < substitutes.size() && allies[i]->stamina < substitutes[i].stamina)
                    allies[i]->substitute(substitutes[i]);
                else allies[i]->doNothing();
        }
        else createPlan(puck, owner < 0);
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
