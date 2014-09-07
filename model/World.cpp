#include "World.h"

using namespace model;
using namespace std;

World::World()
        : tick(-1), tickCount(-1), width(-1.0), height(-1.0), players(vector<Player> ()),
        hockeyists(vector<Hockeyist> ()), puck(Puck ()) { }

World::World(int tick, int tickCount, double width, double height, const vector<Player>& players,
        const vector<Hockeyist>& hockeyists, const Puck& puck)
        : tick(tick), tickCount(tickCount), width(width), height(height), players(players), hockeyists(hockeyists),
        puck(puck) { }

int World::getTick() const {
    return tick;
}

int World::getTickCount() const {
    return tickCount;
}

double World::getWidth() const {
    return width;
}

double World::getHeight() const {
    return height;
}

const vector<Player>& World::getPlayers() const {
    return players;
}

const vector<Hockeyist>& World::getHockeyists() const {
    return hockeyists;
}

const Puck& World::getPuck() const {
    return puck;
}

Player World::getMyPlayer() const {
    for (int playerIndex = (int) players.size() - 1; playerIndex >= 0; --playerIndex) {
        Player player = players[playerIndex];
        if (player.isMe()) {
            return player;
        }
    }

    throw;
}

Player World::getOpponentPlayer() const {
    for (int playerIndex = (int) players.size() - 1; playerIndex >= 0; --playerIndex) {
        Player player = players[playerIndex];
        if (!player.isMe()) {
            return player;
        }
    }

    throw;
}
