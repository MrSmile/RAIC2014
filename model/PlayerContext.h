#pragma once

#ifndef _PLAYER_CONTEXT_H_
#define _PLAYER_CONTEXT_H_

#include <vector>

#include "Hockeyist.h"
#include "World.h"

namespace model {
    class PlayerContext {
    private:
        std::vector<Hockeyist> hockeyists;
        World world;
    public:
        PlayerContext();
        PlayerContext(const std::vector<Hockeyist>& hockeyists, const World& world);

        const std::vector<Hockeyist>& getHockeyists() const;
        const World& getWorld() const;
    };
}

#endif
