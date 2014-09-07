#pragma once

#ifndef _MOVE_H_
#define _MOVE_H_

#include "ActionType.h"

namespace model {
    class Move {
    private:
        double speedUp;
        double turn;
        ActionType action;
        double passPower;
        double passAngle;
        int teammateIndex;
    public:
        Move();

        double getSpeedUp() const;
        void setSpeedUp(const double speedUp);
        double getTurn() const;
        void setTurn(const double turn);
        ActionType getAction() const;
        void setAction(const ActionType action);
        double getPassPower() const;
        void setPassPower(const double passPower);
        double getPassAngle() const;
        void setPassAngle(const double passAngle);
        int getTeammateIndex() const;
        void setTeammateIndex(const int teammateIndex);
    };
}

#endif
