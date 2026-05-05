#include "vect3d.h"

#include <algorithm>
#include <cmath>

Vect3d::Vect3d(void) {
    Zero();
}

Vect3d::Vect3d(float x, float y, float z) {
    Set(x, y, z);
}

Vect3d::Vect3d(const float* newv) {
    Set(newv[0], newv[1], newv[2]);
}

Vect3d::Vect3d(const Vect3d& newv) {
    Set(newv.v[0], newv.v[1], newv.v[2]);
}

void Vect3d::Zero(void) {
    v[0] = v[1] = v[2] = 0.0f;
}

void Vect3d::One(void) {
    v[0] = v[1] = v[2] = 1.0f;
}

void Vect3d::Normalize() {
    float len = Length();
    if (len > 1e-12f) {
        v[0] /= len;
        v[1] /= len;
        v[2] /= len;
    }
}

Vect3d Vect3d::GetNormalized() const {
    Vect3d out(*this);
    out.Normalize();
    return out;
}

float Vect3d::Length() const {
    return std::sqrt(SquaredLength());
}

void Vect3d::RotateX(double angle) {
    *this = GetRotatedX(angle);
}

Vect3d Vect3d::GetRotatedX(double angle) const {
    float c = static_cast<float>(std::cos(angle));
    float s = static_cast<float>(std::sin(angle));
    return Vect3d(v[0], v[1] * c - v[2] * s, v[1] * s + v[2] * c);
}

void Vect3d::RotateY(double angle) {
    *this = GetRotatedY(angle);
}

Vect3d Vect3d::GetRotatedY(double angle) const {
    float c = static_cast<float>(std::cos(angle));
    float s = static_cast<float>(std::sin(angle));
    return Vect3d(v[0] * c + v[2] * s, v[1], -v[0] * s + v[2] * c);
}

void Vect3d::RotateZ(double angle) {
    *this = GetRotatedZ(angle);
}

Vect3d Vect3d::GetRotatedZ(double angle) const {
    float c = static_cast<float>(std::cos(angle));
    float s = static_cast<float>(std::sin(angle));
    return Vect3d(v[0] * c - v[1] * s, v[0] * s + v[1] * c, v[2]);
}

void Vect3d::RotateAxis(double angle, const Vect3d& axis) {
    *this = GetRotatedAxis(angle, axis);
}

Vect3d Vect3d::GetRotatedAxis(double angle, const Vect3d& axis) const {
    float norm = axis.Length();
    if (norm <= 1e-12f) {
        return *this;
    }
    Vect3d k = axis;
    k /= norm;
    return rotateAroundAxis(*this, k, static_cast<float>(angle));
}

void Vect3d::Saturate() {
    v[0] = std::clamp(v[0], 0.0f, 1.0f);
    v[1] = std::clamp(v[1], 0.0f, 1.0f);
    v[2] = std::clamp(v[2], 0.0f, 1.0f);
}

Vect3d Vect3d::GetSaturated() const {
    Vect3d out(*this);
    out.Saturate();
    return out;
}

Vect3d operator*(float scaleFactor, const Vect3d& rhs) {
    return Vect3d(rhs.v[0] * scaleFactor, rhs.v[1] * scaleFactor, rhs.v[2] * scaleFactor);
}

bool Vect3d::operator==(const Vect3d& rhs) const {
    return v[0] == rhs.v[0] && v[1] == rhs.v[1] && v[2] == rhs.v[2];
}

Vect3d Vect3d::rotateAroundAxisThroughPoint(const Vect3d& P,
    const Vect3d& axisPoint,
    const Vect3d& axisDirUnit,
    float angleRad) {
    Vect3d rel = P - axisPoint;
    Vect3d dir = axisDirUnit;
    float dirLen = dir.Length();
    if (dirLen > 1e-12f) {
        dir /= dirLen;
    }
    Vect3d rotated = rotateAroundAxis(rel, dir, angleRad);
    return axisPoint + rotated;
}

