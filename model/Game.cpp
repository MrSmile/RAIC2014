#include "Game.h"

using namespace model;

Game::Game()
        : randomSeed(-1), tickCount(-1), worldWidth(-1.0), worldHeight(-1.0), goalNetTop(-1.0), goalNetWidth(-1.0),
        goalNetHeight(-1.0), rinkTop(-1.0), rinkLeft(-1.0), rinkBottom(-1.0), rinkRight(-1.0),
        afterGoalStateTickCount(-1), overtimeTickCount(-1), defaultActionCooldownTicks(-1),
        swingActionCooldownTicks(-1), cancelStrikeActionCooldownTicks(-1), actionCooldownTicksAfterLosingPuck(-1),
        stickLength(-1.0), stickSector(-1.0), passSector(-1.0), hockeyistAttributeBaseValue(-1), minActionChance(-1.0),
        maxActionChance(-1.0), strikeAngleDeviation(-1.0), passAngleDeviation(-1.0), pickUpPuckBaseChance(-1.0),
        takePuckAwayBaseChance(-1.0), maxEffectiveSwingTicks(-1), strikePowerBaseFactor(-1.0),
        strikePowerGrowthFactor(-1.0), strikePuckBaseChance(-1.0), knockdownChanceFactor(-1.0),
        knockdownTicksFactor(-1.0), maxSpeedToAllowSubstitute(-1.0), substitutionAreaHeight(-1.0),
        passPowerFactor(-1.0), hockeyistMaxStamina(-1.0), activeHockeyistStaminaGrowthPerTick(-1.0),
        restingHockeyistStaminaGrowthPerTick(-1.0), zeroStaminaHockeyistEffectivenessFactor(-1.0),
        speedUpStaminaCostFactor(-1.0), turnStaminaCostFactor(-1.0), takePuckStaminaCost(-1.0), swingStaminaCost(-1.0),
        strikeStaminaBaseCost(-1.0), strikeStaminaCostGrowthFactor(-1.0), cancelStrikeStaminaCost(-1.0),
        passStaminaCost(-1.0), goalieMaxSpeed(-1.0), hockeyistMaxSpeed(-1.0), struckHockeyistInitialSpeedFactor(-1.0),
        hockeyistSpeedUpFactor(-1.0), hockeyistSpeedDownFactor(-1.0), hockeyistTurnAngleFactor(-1.0),
        versatileHockeyistStrength(-1), versatileHockeyistEndurance(-1), versatileHockeyistDexterity(-1),
        versatileHockeyistAgility(-1), forwardHockeyistStrength(-1), forwardHockeyistEndurance(-1),
        forwardHockeyistDexterity(-1), forwardHockeyistAgility(-1), defencemanHockeyistStrength(-1),
        defencemanHockeyistEndurance(-1), defencemanHockeyistDexterity(-1), defencemanHockeyistAgility(-1),
        minRandomHockeyistParameter(-1), maxRandomHockeyistParameter(-1), struckPuckInitialSpeedFactor(-1.0),
        puckBindingRange(-1.0) { }

Game::Game(long long randomSeed, int tickCount, double worldWidth, double worldHeight, double goalNetTop,
        double goalNetWidth, double goalNetHeight, double rinkTop, double rinkLeft, double rinkBottom, double rinkRight,
        int afterGoalStateTickCount, int overtimeTickCount, int defaultActionCooldownTicks,
        int swingActionCooldownTicks, int cancelStrikeActionCooldownTicks, int actionCooldownTicksAfterLosingPuck,
        double stickLength, double stickSector, double passSector, int hockeyistAttributeBaseValue,
        double minActionChance, double maxActionChance, double strikeAngleDeviation, double passAngleDeviation,
        double pickUpPuckBaseChance, double takePuckAwayBaseChance, int maxEffectiveSwingTicks,
        double strikePowerBaseFactor, double strikePowerGrowthFactor, double strikePuckBaseChance,
        double knockdownChanceFactor, double knockdownTicksFactor, double maxSpeedToAllowSubstitute,
        double substitutionAreaHeight, double passPowerFactor, double hockeyistMaxStamina,
        double activeHockeyistStaminaGrowthPerTick, double restingHockeyistStaminaGrowthPerTick,
        double zeroStaminaHockeyistEffectivenessFactor, double speedUpStaminaCostFactor, double turnStaminaCostFactor,
        double takePuckStaminaCost, double swingStaminaCost, double strikeStaminaBaseCost,
        double strikeStaminaCostGrowthFactor, double cancelStrikeStaminaCost, double passStaminaCost,
        double goalieMaxSpeed, double hockeyistMaxSpeed, double struckHockeyistInitialSpeedFactor,
        double hockeyistSpeedUpFactor, double hockeyistSpeedDownFactor, double hockeyistTurnAngleFactor,
        int versatileHockeyistStrength, int versatileHockeyistEndurance, int versatileHockeyistDexterity,
        int versatileHockeyistAgility, int forwardHockeyistStrength, int forwardHockeyistEndurance,
        int forwardHockeyistDexterity, int forwardHockeyistAgility, int defencemanHockeyistStrength,
        int defencemanHockeyistEndurance, int defencemanHockeyistDexterity, int defencemanHockeyistAgility,
        int minRandomHockeyistParameter, int maxRandomHockeyistParameter, double struckPuckInitialSpeedFactor,
        double puckBindingRange)
        : randomSeed(randomSeed), tickCount(tickCount), worldWidth(worldWidth), worldHeight(worldHeight),
        goalNetTop(goalNetTop), goalNetWidth(goalNetWidth), goalNetHeight(goalNetHeight), rinkTop(rinkTop),
        rinkLeft(rinkLeft), rinkBottom(rinkBottom), rinkRight(rinkRight),
        afterGoalStateTickCount(afterGoalStateTickCount), overtimeTickCount(overtimeTickCount),
        defaultActionCooldownTicks(defaultActionCooldownTicks), swingActionCooldownTicks(swingActionCooldownTicks),
        cancelStrikeActionCooldownTicks(cancelStrikeActionCooldownTicks),
        actionCooldownTicksAfterLosingPuck(actionCooldownTicksAfterLosingPuck), stickLength(stickLength),
        stickSector(stickSector), passSector(passSector), hockeyistAttributeBaseValue(hockeyistAttributeBaseValue),
        minActionChance(minActionChance), maxActionChance(maxActionChance), strikeAngleDeviation(strikeAngleDeviation),
        passAngleDeviation(passAngleDeviation), pickUpPuckBaseChance(pickUpPuckBaseChance),
        takePuckAwayBaseChance(takePuckAwayBaseChance), maxEffectiveSwingTicks(maxEffectiveSwingTicks),
        strikePowerBaseFactor(strikePowerBaseFactor), strikePowerGrowthFactor(strikePowerGrowthFactor),
        strikePuckBaseChance(strikePuckBaseChance), knockdownChanceFactor(knockdownChanceFactor),
        knockdownTicksFactor(knockdownTicksFactor), maxSpeedToAllowSubstitute(maxSpeedToAllowSubstitute),
        substitutionAreaHeight(substitutionAreaHeight), passPowerFactor(passPowerFactor),
        hockeyistMaxStamina(hockeyistMaxStamina),
        activeHockeyistStaminaGrowthPerTick(activeHockeyistStaminaGrowthPerTick),
        restingHockeyistStaminaGrowthPerTick(restingHockeyistStaminaGrowthPerTick),
        zeroStaminaHockeyistEffectivenessFactor(zeroStaminaHockeyistEffectivenessFactor),
        speedUpStaminaCostFactor(speedUpStaminaCostFactor), turnStaminaCostFactor(turnStaminaCostFactor),
        takePuckStaminaCost(takePuckStaminaCost), swingStaminaCost(swingStaminaCost),
        strikeStaminaBaseCost(strikeStaminaBaseCost), strikeStaminaCostGrowthFactor(strikeStaminaCostGrowthFactor),
        cancelStrikeStaminaCost(cancelStrikeStaminaCost), passStaminaCost(passStaminaCost),
        goalieMaxSpeed(goalieMaxSpeed), hockeyistMaxSpeed(hockeyistMaxSpeed),
        struckHockeyistInitialSpeedFactor(struckHockeyistInitialSpeedFactor),
        hockeyistSpeedUpFactor(hockeyistSpeedUpFactor), hockeyistSpeedDownFactor(hockeyistSpeedDownFactor),
        hockeyistTurnAngleFactor(hockeyistTurnAngleFactor), versatileHockeyistStrength(versatileHockeyistStrength),
        versatileHockeyistEndurance(versatileHockeyistEndurance),
        versatileHockeyistDexterity(versatileHockeyistDexterity), versatileHockeyistAgility(versatileHockeyistAgility),
        forwardHockeyistStrength(forwardHockeyistStrength), forwardHockeyistEndurance(forwardHockeyistEndurance),
        forwardHockeyistDexterity(forwardHockeyistDexterity), forwardHockeyistAgility(forwardHockeyistAgility),
        defencemanHockeyistStrength(defencemanHockeyistStrength),
        defencemanHockeyistEndurance(defencemanHockeyistEndurance),
        defencemanHockeyistDexterity(defencemanHockeyistDexterity),
        defencemanHockeyistAgility(defencemanHockeyistAgility),
        minRandomHockeyistParameter(minRandomHockeyistParameter),
        maxRandomHockeyistParameter(maxRandomHockeyistParameter),
        struckPuckInitialSpeedFactor(struckPuckInitialSpeedFactor), puckBindingRange(puckBindingRange) { }

long long Game::getRandomSeed() const {
    return randomSeed;
}

int Game::getTickCount() const {
    return tickCount;
}

double Game::getWorldWidth() const {
    return worldWidth;
}

double Game::getWorldHeight() const {
    return worldHeight;
}

double Game::getGoalNetTop() const {
    return goalNetTop;
}

double Game::getGoalNetWidth() const {
    return goalNetWidth;
}

double Game::getGoalNetHeight() const {
    return goalNetHeight;
}

double Game::getRinkTop() const {
    return rinkTop;
}

double Game::getRinkLeft() const {
    return rinkLeft;
}

double Game::getRinkBottom() const {
    return rinkBottom;
}

double Game::getRinkRight() const {
    return rinkRight;
}

int Game::getAfterGoalStateTickCount() const {
    return afterGoalStateTickCount;
}

int Game::getOvertimeTickCount() const {
    return overtimeTickCount;
}

int Game::getDefaultActionCooldownTicks() const {
    return defaultActionCooldownTicks;
}

int Game::getSwingActionCooldownTicks() const {
    return swingActionCooldownTicks;
}

int Game::getCancelStrikeActionCooldownTicks() const {
    return cancelStrikeActionCooldownTicks;
}

int Game::getActionCooldownTicksAfterLosingPuck() const {
    return actionCooldownTicksAfterLosingPuck;
}

double Game::getStickLength() const {
    return stickLength;
}

double Game::getStickSector() const {
    return stickSector;
}

double Game::getPassSector() const {
    return passSector;
}

int Game::getHockeyistAttributeBaseValue() const {
    return hockeyistAttributeBaseValue;
}

double Game::getMinActionChance() const {
    return minActionChance;
}

double Game::getMaxActionChance() const {
    return maxActionChance;
}

double Game::getStrikeAngleDeviation() const {
    return strikeAngleDeviation;
}

double Game::getPassAngleDeviation() const {
    return passAngleDeviation;
}

double Game::getPickUpPuckBaseChance() const {
    return pickUpPuckBaseChance;
}

double Game::getTakePuckAwayBaseChance() const {
    return takePuckAwayBaseChance;
}

int Game::getMaxEffectiveSwingTicks() const {
    return maxEffectiveSwingTicks;
}

double Game::getStrikePowerBaseFactor() const {
    return strikePowerBaseFactor;
}

double Game::getStrikePowerGrowthFactor() const {
    return strikePowerGrowthFactor;
}

double Game::getStrikePuckBaseChance() const {
    return strikePuckBaseChance;
}

double Game::getKnockdownChanceFactor() const {
    return knockdownChanceFactor;
}

double Game::getKnockdownTicksFactor() const {
    return knockdownTicksFactor;
}

double Game::getMaxSpeedToAllowSubstitute() const {
    return maxSpeedToAllowSubstitute;
}

double Game::getSubstitutionAreaHeight() const {
    return substitutionAreaHeight;
}

double Game::getPassPowerFactor() const {
    return passPowerFactor;
}

double Game::getHockeyistMaxStamina() const {
    return hockeyistMaxStamina;
}

double Game::getActiveHockeyistStaminaGrowthPerTick() const {
    return activeHockeyistStaminaGrowthPerTick;
}

double Game::getRestingHockeyistStaminaGrowthPerTick() const {
    return restingHockeyistStaminaGrowthPerTick;
}

double Game::getZeroStaminaHockeyistEffectivenessFactor() const {
    return zeroStaminaHockeyistEffectivenessFactor;
}

double Game::getSpeedUpStaminaCostFactor() const {
    return speedUpStaminaCostFactor;
}

double Game::getTurnStaminaCostFactor() const {
    return turnStaminaCostFactor;
}

double Game::getTakePuckStaminaCost() const {
    return takePuckStaminaCost;
}

double Game::getSwingStaminaCost() const {
    return swingStaminaCost;
}

double Game::getStrikeStaminaBaseCost() const {
    return strikeStaminaBaseCost;
}

double Game::getStrikeStaminaCostGrowthFactor() const {
    return strikeStaminaCostGrowthFactor;
}

double Game::getCancelStrikeStaminaCost() const {
    return cancelStrikeStaminaCost;
}

double Game::getPassStaminaCost() const {
    return passStaminaCost;
}

double Game::getGoalieMaxSpeed() const {
    return goalieMaxSpeed;
}

double Game::getHockeyistMaxSpeed() const {
    return hockeyistMaxSpeed;
}

double Game::getStruckHockeyistInitialSpeedFactor() const {
    return struckHockeyistInitialSpeedFactor;
}

double Game::getHockeyistSpeedUpFactor() const {
    return hockeyistSpeedUpFactor;
}

double Game::getHockeyistSpeedDownFactor() const {
    return hockeyistSpeedDownFactor;
}

double Game::getHockeyistTurnAngleFactor() const {
    return hockeyistTurnAngleFactor;
}

int Game::getVersatileHockeyistStrength() const {
    return versatileHockeyistStrength;
}

int Game::getVersatileHockeyistEndurance() const {
    return versatileHockeyistEndurance;
}

int Game::getVersatileHockeyistDexterity() const {
    return versatileHockeyistDexterity;
}

int Game::getVersatileHockeyistAgility() const {
    return versatileHockeyistAgility;
}

int Game::getForwardHockeyistStrength() const {
    return forwardHockeyistStrength;
}

int Game::getForwardHockeyistEndurance() const {
    return forwardHockeyistEndurance;
}

int Game::getForwardHockeyistDexterity() const {
    return forwardHockeyistDexterity;
}

int Game::getForwardHockeyistAgility() const {
    return forwardHockeyistAgility;
}

int Game::getDefencemanHockeyistStrength() const {
    return defencemanHockeyistStrength;
}

int Game::getDefencemanHockeyistEndurance() const {
    return defencemanHockeyistEndurance;
}

int Game::getDefencemanHockeyistDexterity() const {
    return defencemanHockeyistDexterity;
}

int Game::getDefencemanHockeyistAgility() const {
    return defencemanHockeyistAgility;
}

int Game::getMinRandomHockeyistParameter() const {
    return minRandomHockeyistParameter;
}

int Game::getMaxRandomHockeyistParameter() const {
    return maxRandomHockeyistParameter;
}

double Game::getStruckPuckInitialSpeedFactor() const {
    return struckPuckInitialSpeedFactor;
}

double Game::getPuckBindingRange() const {
    return puckBindingRange;
}
