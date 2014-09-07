#pragma once

#ifndef _WORLD_H_
#define _WORLD_H_

#include <vector>

#include "Hockeyist.h"
#include "Player.h"
#include "Puck.h"

namespace model {
    class World {
    private:
        int tick;
        int tickCount;
        double width;
        double height;
        std::vector<Player> players;
        std::vector<Hockeyist> hockeyists;
        Puck puck;
    public:
        World();
        World(int tick, int tickCount, double width, double height, const std::vector<Player>& players,
                const std::vector<Hockeyist>& hockeyists, const Puck& puck);

        int getTick() const;
        int getTickCount() const;
        double getWidth() const;
        double getHeight() const;
        const std::vector<Player>& getPlayers() const;
        const std::vector<Hockeyist>& getHockeyists() const;
        const Puck& getPuck() const;

        Player getMyPlayer() const;
        Player getOpponentPlayer() const;
    };
}

#endif
