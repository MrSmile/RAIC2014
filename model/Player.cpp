#include "Player.h"

using namespace model;
using namespace std;

Player::Player()
        : id(-1), me(false), name(""), goalCount(-1), strategyCrashed(false), netTop(-1.0), netLeft(-1.0),
        netBottom(-1.0), netRight(-1.0), netFront(-1.0), netBack(-1.0), justScoredGoal(false), justMissedGoal(false) { }

Player::Player(long long id, bool me, const string& name, int goalCount, bool strategyCrashed, double netTop,
        double netLeft, double netBottom, double netRight, double netFront, double netBack, bool justScoredGoal,
        bool justMissedGoal)
        : id(id), me(me), name(name), goalCount(goalCount), strategyCrashed(strategyCrashed), netTop(netTop),
        netLeft(netLeft), netBottom(netBottom), netRight(netRight), netFront(netFront), netBack(netBack),
        justScoredGoal(justScoredGoal), justMissedGoal(justMissedGoal) { }

long long Player::getId() const {
    return id;
}

bool Player::isMe() const {
    return me;
}

const string& Player::getName() const {
    return name;
}

int Player::getGoalCount() const {
    return goalCount;
}

bool Player::isStrategyCrashed() const {
    return strategyCrashed;
}

double Player::getNetTop() const {
    return netTop;
}

double Player::getNetLeft() const {
    return netLeft;
}

double Player::getNetBottom() const {
    return netBottom;
}

double Player::getNetRight() const {
    return netRight;
}

double Player::getNetFront() const {
    return netFront;
}

double Player::getNetBack() const {
    return netBack;
}

bool Player::isJustScoredGoal() const {
    return justScoredGoal;
}

bool Player::isJustMissedGoal() const {
    return justMissedGoal;
}
