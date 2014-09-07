#pragma once

#ifndef _HOCKEYIST_STATE_H_
#define _HOCKEYIST_STATE_H_

namespace model {
    enum HockeyistState {
        _UNKNOWN_HOCKEYIST_STATE_ = -1,
        ACTIVE = 0,
        SWINGING = 1,
        KNOCKED_DOWN = 2,
        RESTING = 3,
        _HOCKEYIST_STATE_COUNT_ = 4
    };
}

#endif
