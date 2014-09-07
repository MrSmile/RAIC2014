#pragma once

#ifndef _HOCKEYIST_H_
#define _HOCKEYIST_H_

#include "ActionType.h"
#include "HockeyistState.h"
#include "HockeyistType.h"
#include "Unit.h"

namespace model {
    class Hockeyist : public Unit {
    private:
        long long playerId;
        int teammateIndex;
        bool teammate;
        HockeyistType type;
        int strength;
        int endurance;
        int dexterity;
        int agility;
        double stamina;
        HockeyistState state;
        int originalPositionIndex;
        int remainingKnockdownTicks;
        int remainingCooldownTicks;
        int swingTicks;
        ActionType lastAction;
        int lastActionTick;
    public:
        Hockeyist();
        Hockeyist(long long id, long long playerId, int teammateIndex, double mass, double radius, double x, double y,
                double speedX, double speedY, double angle, double angularSpeed, bool teammate, HockeyistType type,
                int strength, int endurance, int dexterity, int agility, double stamina, HockeyistState state,
                int originalPositionIndex, int remainingKnockdownTicks, int remainingCooldownTicks, int swingTicks,
                ActionType lastAction, int lastActionTick);

        long long getPlayerId() const;
        int getTeammateIndex() const;
        bool isTeammate() const;
        HockeyistType getType() const;
        int getStrength() const;
        int getEndurance() const;
        int getDexterity() const;
        int getAgility() const;
        double getStamina() const;
        HockeyistState getState() const;
        int getOriginalPositionIndex() const;
        int getRemainingKnockdownTicks() const;
        int getRemainingCooldownTicks() const;
        int getSwingTicks() const;
        ActionType getLastAction() const;
        int getLastActionTick() const;
    };
}

#endif
