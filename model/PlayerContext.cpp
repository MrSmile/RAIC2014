#include "PlayerContext.h"

using namespace model;
using namespace std;

PlayerContext::PlayerContext()
        : hockeyists(vector<Hockeyist> ()), world(World ()) { }

PlayerContext::PlayerContext(const vector<Hockeyist>& hockeyists, const World& world)
        : hockeyists(hockeyists), world(world) { }

const vector<Hockeyist>& PlayerContext::getHockeyists() const {
    return hockeyists;
}

const World& PlayerContext::getWorld() const {
    return world;
}
