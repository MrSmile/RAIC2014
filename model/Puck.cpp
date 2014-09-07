#include "Puck.h"

using namespace model;

Puck::Puck()
        : Unit(0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0), ownerHockeyistId(-1), ownerPlayerId(-1) { }

Puck::Puck(long long id, double mass, double radius, double x, double y, double speedX, double speedY,
        long long ownerHockeyistId, long long ownerPlayerId)
        : Unit(id, mass, radius, x, y, speedX, speedY, 0.0, 0.0), ownerHockeyistId(ownerHockeyistId),
        ownerPlayerId(ownerPlayerId) { }

long long Puck::getOwnerHockeyistId() const {
    return ownerHockeyistId;
}

long long Puck::getOwnerPlayerId() const {
    return ownerPlayerId;
}
