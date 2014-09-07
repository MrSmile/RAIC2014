#include "MyStrategy.h"

#define PI 3.14159265358979323846
#define _USE_MATH_DEFINES

#include <cmath>
#include <cstdlib>

using namespace model;
using namespace std;

void MyStrategy::move(const Hockeyist& self, const World& world, const Game& game, Move& move) {
    move.setSpeedUp(-1.0);
    move.setTurn(PI);
    move.setAction(STRIKE);
}

MyStrategy::MyStrategy() { }
