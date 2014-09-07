#pragma once

#ifndef _PUCK_H_
#define _PUCK_H_

#include "Unit.h"

namespace model {
    class Puck : public Unit {
    private:
        long long ownerHockeyistId;
        long long ownerPlayerId;
    public:
        Puck();
        Puck(long long id, double mass, double radius, double x, double y, double speedX, double speedY,
                long long ownerHockeyistId, long long ownerPlayerId);

        long long getOwnerHockeyistId() const;
        long long getOwnerPlayerId() const;
    };
}

#endif
