#pragma once

#ifndef _GAME_H_
#define _GAME_H_

namespace model {
    class Game {
    private:
        long long randomSeed;
        int tickCount;
        double worldWidth;
        double worldHeight;
        double goalNetTop;
        double goalNetWidth;
        double goalNetHeight;
        double rinkTop;
        double rinkLeft;
        double rinkBottom;
        double rinkRight;
        int afterGoalStateTickCount;
        int overtimeTickCount;
        int defaultActionCooldownTicks;
        int swingActionCooldownTicks;
        int cancelStrikeActionCooldownTicks;
        int actionCooldownTicksAfterLosingPuck;
        double stickLength;
        double stickSector;
        double passSector;
        int hockeyistAttributeBaseValue;
        double minActionChance;
        double maxActionChance;
        double strikeAngleDeviation;
        double passAngleDeviation;
        double pickUpPuckBaseChance;
        double takePuckAwayBaseChance;
        int maxEffectiveSwingTicks;
        double strikePowerBaseFactor;
        double strikePowerGrowthFactor;
        double strikePuckBaseChance;
        double knockdownChanceFactor;
        double knockdownTicksFactor;
        double maxSpeedToAllowSubstitute;
        double substitutionAreaHeight;
        double passPowerFactor;
        double hockeyistMaxStamina;
        double activeHockeyistStaminaGrowthPerTick;
        double restingHockeyistStaminaGrowthPerTick;
        double zeroStaminaHockeyistEffectivenessFactor;
        double speedUpStaminaCostFactor;
        double turnStaminaCostFactor;
        double takePuckStaminaCost;
        double swingStaminaCost;
        double strikeStaminaBaseCost;
        double strikeStaminaCostGrowthFactor;
        double cancelStrikeStaminaCost;
        double passStaminaCost;
        double goalieMaxSpeed;
        double hockeyistMaxSpeed;
        double struckHockeyistInitialSpeedFactor;
        double hockeyistSpeedUpFactor;
        double hockeyistSpeedDownFactor;
        double hockeyistTurnAngleFactor;
        int versatileHockeyistStrength;
        int versatileHockeyistEndurance;
        int versatileHockeyistDexterity;
        int versatileHockeyistAgility;
        int forwardHockeyistStrength;
        int forwardHockeyistEndurance;
        int forwardHockeyistDexterity;
        int forwardHockeyistAgility;
        int defencemanHockeyistStrength;
        int defencemanHockeyistEndurance;
        int defencemanHockeyistDexterity;
        int defencemanHockeyistAgility;
        int minRandomHockeyistParameter;
        int maxRandomHockeyistParameter;
        double struckPuckInitialSpeedFactor;
        double puckBindingRange;
    public:
        Game();
        Game(long long randomSeed, int tickCount, double worldWidth, double worldHeight, double goalNetTop,
                double goalNetWidth, double goalNetHeight, double rinkTop, double rinkLeft, double rinkBottom,
                double rinkRight, int afterGoalStateTickCount, int overtimeTickCount, int defaultActionCooldownTicks,
                int swingActionCooldownTicks, int cancelStrikeActionCooldownTicks,
                int actionCooldownTicksAfterLosingPuck, double stickLength, double stickSector, double passSector,
                int hockeyistAttributeBaseValue, double minActionChance, double maxActionChance,
                double strikeAngleDeviation, double passAngleDeviation, double pickUpPuckBaseChance,
                double takePuckAwayBaseChance, int maxEffectiveSwingTicks, double strikePowerBaseFactor,
                double strikePowerGrowthFactor, double strikePuckBaseChance, double knockdownChanceFactor,
                double knockdownTicksFactor, double maxSpeedToAllowSubstitute, double substitutionAreaHeight,
                double passPowerFactor, double hockeyistMaxStamina, double activeHockeyistStaminaGrowthPerTick,
                double restingHockeyistStaminaGrowthPerTick, double zeroStaminaHockeyistEffectivenessFactor,
                double speedUpStaminaCostFactor, double turnStaminaCostFactor, double takePuckStaminaCost,
                double swingStaminaCost, double strikeStaminaBaseCost, double strikeStaminaCostGrowthFactor,
                double cancelStrikeStaminaCost, double passStaminaCost, double goalieMaxSpeed, double hockeyistMaxSpeed,
                double struckHockeyistInitialSpeedFactor, double hockeyistSpeedUpFactor,
                double hockeyistSpeedDownFactor, double hockeyistTurnAngleFactor, int versatileHockeyistStrength,
                int versatileHockeyistEndurance, int versatileHockeyistDexterity, int versatileHockeyistAgility,
                int forwardHockeyistStrength, int forwardHockeyistEndurance, int forwardHockeyistDexterity,
                int forwardHockeyistAgility, int defencemanHockeyistStrength, int defencemanHockeyistEndurance,
                int defencemanHockeyistDexterity, int defencemanHockeyistAgility, int minRandomHockeyistParameter,
                int maxRandomHockeyistParameter, double struckPuckInitialSpeedFactor, double puckBindingRange);

        long long getRandomSeed() const;
        int getTickCount() const;
        double getWorldWidth() const;
        double getWorldHeight() const;
        double getGoalNetTop() const;
        double getGoalNetWidth() const;
        double getGoalNetHeight() const;
        double getRinkTop() const;
        double getRinkLeft() const;
        double getRinkBottom() const;
        double getRinkRight() const;
        int getAfterGoalStateTickCount() const;
        int getOvertimeTickCount() const;
        int getDefaultActionCooldownTicks() const;
        int getSwingActionCooldownTicks() const;
        int getCancelStrikeActionCooldownTicks() const;
        int getActionCooldownTicksAfterLosingPuck() const;
        double getStickLength() const;
        double getStickSector() const;
        double getPassSector() const;
        int getHockeyistAttributeBaseValue() const;
        double getMinActionChance() const;
        double getMaxActionChance() const;
        double getStrikeAngleDeviation() const;
        double getPassAngleDeviation() const;
        double getPickUpPuckBaseChance() const;
        double getTakePuckAwayBaseChance() const;
        int getMaxEffectiveSwingTicks() const;
        double getStrikePowerBaseFactor() const;
        double getStrikePowerGrowthFactor() const;
        double getStrikePuckBaseChance() const;
        double getKnockdownChanceFactor() const;
        double getKnockdownTicksFactor() const;
        double getMaxSpeedToAllowSubstitute() const;
        double getSubstitutionAreaHeight() const;
        double getPassPowerFactor() const;
        double getHockeyistMaxStamina() const;
        double getActiveHockeyistStaminaGrowthPerTick() const;
        double getRestingHockeyistStaminaGrowthPerTick() const;
        double getZeroStaminaHockeyistEffectivenessFactor() const;
        double getSpeedUpStaminaCostFactor() const;
        double getTurnStaminaCostFactor() const;
        double getTakePuckStaminaCost() const;
        double getSwingStaminaCost() const;
        double getStrikeStaminaBaseCost() const;
        double getStrikeStaminaCostGrowthFactor() const;
        double getCancelStrikeStaminaCost() const;
        double getPassStaminaCost() const;
        double getGoalieMaxSpeed() const;
        double getHockeyistMaxSpeed() const;
        double getStruckHockeyistInitialSpeedFactor() const;
        double getHockeyistSpeedUpFactor() const;
        double getHockeyistSpeedDownFactor() const;
        double getHockeyistTurnAngleFactor() const;
        int getVersatileHockeyistStrength() const;
        int getVersatileHockeyistEndurance() const;
        int getVersatileHockeyistDexterity() const;
        int getVersatileHockeyistAgility() const;
        int getForwardHockeyistStrength() const;
        int getForwardHockeyistEndurance() const;
        int getForwardHockeyistDexterity() const;
        int getForwardHockeyistAgility() const;
        int getDefencemanHockeyistStrength() const;
        int getDefencemanHockeyistEndurance() const;
        int getDefencemanHockeyistDexterity() const;
        int getDefencemanHockeyistAgility() const;
        int getMinRandomHockeyistParameter() const;
        int getMaxRandomHockeyistParameter() const;
        double getStruckPuckInitialSpeedFactor() const;
        double getPuckBindingRange() const;
    };
}

#endif
