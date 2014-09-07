#pragma once

#ifndef _PLAYER_H_
#define _PLAYER_H_

#include <string>

namespace model {
    class Player {
    private:
        long long id;
        bool me;
        std::string name;
        int goalCount;
        bool strategyCrashed;
        double netTop;
        double netLeft;
        double netBottom;
        double netRight;
        double netFront;
        double netBack;
        bool justScoredGoal;
        bool justMissedGoal;
    public:
        Player();
        Player(long long id, bool me, const std::string& name, int goalCount, bool strategyCrashed, double netTop,
                double netLeft, double netBottom, double netRight, double netFront, double netBack, bool justScoredGoal,
                bool justMissedGoal);

        long long getId() const;
        bool isMe() const;
        const std::string& getName() const;
        int getGoalCount() const;
        bool isStrategyCrashed() const;
        double getNetTop() const;
        double getNetLeft() const;
        double getNetBottom() const;
        double getNetRight() const;
        double getNetFront() const;
        double getNetBack() const;
        bool isJustScoredGoal() const;
        bool isJustMissedGoal() const;
    };
}

#endif
