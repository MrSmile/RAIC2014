#include "Hockeyist.h"

using namespace model;

Hockeyist::Hockeyist()
        : Unit(0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0), playerId(-1), teammateIndex(-1), teammate(false),
        type(_UNKNOWN_HOCKEYIST_TYPE_), strength(-1), endurance(-1), dexterity(-1), agility(-1), stamina(-1.0),
        state(_UNKNOWN_HOCKEYIST_STATE_), originalPositionIndex(-1), remainingKnockdownTicks(-1),
        remainingCooldownTicks(-1), swingTicks(-1), lastAction(_UNKNOWN_ACTION_TYPE_), lastActionTick(-1) { }

Hockeyist::Hockeyist(long long id, long long playerId, int teammateIndex, double mass, double radius, double x,
        double y, double speedX, double speedY, double angle, double angularSpeed, bool teammate, HockeyistType type,
        int strength, int endurance, int dexterity, int agility, double stamina, HockeyistState state,
        int originalPositionIndex, int remainingKnockdownTicks, int remainingCooldownTicks, int swingTicks,
        ActionType lastAction, int lastActionTick)
        : Unit(id, mass, radius, x, y, speedX, speedY, angle, angularSpeed), playerId(playerId),
        teammateIndex(teammateIndex), teammate(teammate), type(type), strength(strength), endurance(endurance),
        dexterity(dexterity), agility(agility), stamina(stamina), state(state),
        originalPositionIndex(originalPositionIndex), remainingKnockdownTicks(remainingKnockdownTicks),
        remainingCooldownTicks(remainingCooldownTicks), swingTicks(swingTicks), lastAction(lastAction),
        lastActionTick(lastActionTick) { }

long long Hockeyist::getPlayerId() const {
    return playerId;
}

int Hockeyist::getTeammateIndex() const {
    return teammateIndex;
}

bool Hockeyist::isTeammate() const {
    return teammate;
}

HockeyistType Hockeyist::getType() const {
    return type;
}

int Hockeyist::getStrength() const {
    return strength;
}

int Hockeyist::getEndurance() const {
    return endurance;
}

int Hockeyist::getDexterity() const {
    return dexterity;
}

int Hockeyist::getAgility() const {
    return agility;
}

double Hockeyist::getStamina() const {
    return stamina;
}

HockeyistState Hockeyist::getState() const {
    return state;
}

int Hockeyist::getOriginalPositionIndex() const {
    return originalPositionIndex;
}

int Hockeyist::getRemainingKnockdownTicks() const {
    return remainingKnockdownTicks;
}

int Hockeyist::getRemainingCooldownTicks() const {
    return remainingCooldownTicks;
}

int Hockeyist::getSwingTicks() const {
    return swingTicks;
}

ActionType Hockeyist::getLastAction() const {
    return lastAction;
}

int Hockeyist::getLastActionTick() const {
    return lastActionTick;
}
