#pragma once

#ifndef _ACTION_TYPE_H_
#define _ACTION_TYPE_H_

namespace model {
    enum ActionType {
        _UNKNOWN_ACTION_TYPE_ = -1,
        NONE = 0,
        TAKE_PUCK = 1,
        SWING = 2,
        STRIKE = 3,
        CANCEL_STRIKE = 4,
        PASS = 5,
        SUBSTITUTE = 6,
        _ACTION_TYPE_COUNT_ = 7
    };
}

#endif
