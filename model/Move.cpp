#include "Move.h"

using namespace model;

Move::Move()
        : speedUp(0.0), turn(0.0), action(NONE), passPower(1.0), passAngle(0.0), teammateIndex(-1) { }

double Move::getSpeedUp() const {
    return speedUp;
}

void Move::setSpeedUp(const double speedUp) {
    this->speedUp = speedUp;
}

double Move::getTurn() const {
    return turn;
}

void Move::setTurn(const double turn) {
    this->turn = turn;
}

ActionType Move::getAction() const {
    return action;
}

void Move::setAction(const ActionType action) {
    this->action = action;
}

double Move::getPassPower() const {
    return passPower;
}

void Move::setPassPower(const double passPower) {
    this->passPower = passPower;
}

double Move::getPassAngle() const {
    return passAngle;
}

void Move::setPassAngle(const double passAngle) {
    this->passAngle = passAngle;
}

int Move::getTeammateIndex() const {
    return teammateIndex;
}

void Move::setTeammateIndex(const int teammateIndex) {
    this->teammateIndex = teammateIndex;
}
