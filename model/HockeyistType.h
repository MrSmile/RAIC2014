#pragma once

#ifndef _HOCKEYIST_TYPE_H_
#define _HOCKEYIST_TYPE_H_

namespace model {
    enum HockeyistType {
        _UNKNOWN_HOCKEYIST_TYPE_ = -1,
        GOALIE = 0,
        VERSATILE = 1,
        FORWARD = 2,
        DEFENCEMAN = 3,
        RANDOM = 4,
        _HOCKEYIST_TYPE_COUNT_ = 5
    };
}

#endif
