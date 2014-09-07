#pragma once

#ifndef _UNIT_H_
#define _UNIT_H_

namespace model {
    class Unit {
    private:
        long long id;
        double mass;
        double radius;
        double x;
        double y;
        double speedX;
        double speedY;
        double angle;
        double angularSpeed;
    protected:
        Unit(long long id, double mass, double radius, double x, double y, double speedX, double speedY, double angle,
                double angularSpeed);
    public:
        long long getId() const;
        double getMass() const;
        double getRadius() const;
        double getX() const;
        double getY() const;
        double getSpeedX() const;
        double getSpeedY() const;
        double getAngle() const;
        double getAngularSpeed() const;

        double getAngleTo(double x, double y) const;
        double getAngleTo(const Unit& unit) const;
        double getDistanceTo(double x, double y) const;
        double getDistanceTo(const Unit& unit) const;

        virtual ~Unit();
    };
}

#endif
