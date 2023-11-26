#include "intersect.h"
// Suppress warnings in third-party code.
#include <framework/disable_all_warnings.h>
DISABLE_WARNINGS_PUSH()
#include <glm/geometric.hpp>
#include <glm/gtx/component_wise.hpp>
#include <glm/vector_relational.hpp>
DISABLE_WARNINGS_POP()
#include <cmath>
#include <limits>


bool pointInTriangle(const glm::vec3& v0, const glm::vec3& v1, const glm::vec3& v2, const glm::vec3& n, const glm::vec3& p) {
    float area = glm::length(glm::cross(v1 - v0, v2 - v0)) / 2;
    float a = (glm::length(glm::cross(p - v1, v2 - v1)) / 2) / area;
    float b = (glm::length(glm::cross(p - v0, v1 - v0)) / 2) / area;
    float c = (glm::length(glm::cross(p - v2, v0 - v2)) / 2) / area;

    if (a + b + c > 1) {
        return false;
    }

    return true;
}

bool intersectRayWithPlane(const Plane& plane, Ray& ray)
{
    float dot = glm::dot(plane.normal, ray.direction);
    if (dot == 0) {
        return false;
    }
    float t = (plane.D - glm::dot(plane.normal, ray.origin)) / dot;
    if (t <= 0) {
        return false;
    }
    if (t < ray.t) {
        ray.t = t;
        return true;
    }
    return false;
}

Plane trianglePlane(const glm::vec3& v0, const glm::vec3& v1, const glm::vec3& v2)
{
    Plane plane;
    glm::vec3 v = v0 - v2;
    glm::vec3 u = v1 - v2;
    glm::vec3 n = glm::normalize(glm::cross(v, u));
    float d = glm::dot(n, v0);

    plane.D = d;
    plane.normal = n;

    return plane;
}

/// Input: the three vertices of the triangle
/// Output: if intersects then modify the hit parameter ray.t and return true, otherwise return false
bool intersectRayWithTriangle(const glm::vec3& v0, const glm::vec3& v1, const glm::vec3& v2, Ray& ray, HitInfo& hitInfo)
{
    Plane plane = trianglePlane(v0, v1, v2);
    float t = ray.t;
    if (intersectRayWithPlane(plane, ray)) {
        if (pointInTriangle(v0, v1, v2, plane.normal, ray.origin + ray.t * ray.direction)) {
            return true;
        }
        ray.t = t;
    }
    return false;
}

/// Input: a sphere with the following attributes: sphere.radius, sphere.center
/// Output: if intersects then modify the hit parameter ray.t and return true, otherwise return false
bool intersectRayWithShape(const Sphere& sphere, Ray& ray, HitInfo& hitInfo)
{
    float f = glm::pow(2 * glm::dot(ray.direction, (ray.origin - sphere.center)), 2) - 4 * glm::pow(glm::length(ray.direction), 2) * (glm::pow(glm::length(ray.origin - sphere.center), 2) - glm::pow(sphere.radius, 2));
    if (f < 0) {
        return false;
    }
    float d = -2 * (glm::dot(ray.direction, ray.origin - sphere.center));

    float newt = (d - glm::sqrt(f)) / (2 * glm::pow(glm::length(ray.direction), 2));

    if (newt <= 0) {
        newt = (d + glm::sqrt(f)) / (2 * glm::pow(glm::length(ray.direction), 2));
    }

    if (newt <= 0 || newt > ray.t) {
        return false;
    }
    ray.t = newt;
    return true;
}

/// Input: an axis-aligned bounding box with the following parameters: minimum coordinates box.lower and maximum coordinates box.upper
/// Output: if intersects then modify the hit parameter ray.t and return true, otherwise return false
bool intersectRayWithShape(const AxisAlignedBox& box, Ray& ray)
{
    float txmin = (box.lower.x - ray.origin.x) / ray.direction.x;
    float tymin = (box.lower.y - ray.origin.y) / ray.direction.y;
    float tzmin = (box.lower.z - ray.origin.z) / ray.direction.z;
    float txmax = (box.upper.x - ray.origin.x) / ray.direction.x;
    float tymax = (box.upper.y - ray.origin.y) / ray.direction.y;
    float tzmax = (box.upper.z - ray.origin.z) / ray.direction.z;

    if (ray.direction.x == 0) {
        txmin = std::numeric_limits<float>::min();
        txmax = std::numeric_limits<float>::max();
    }
    if (ray.direction.y == 0) {
        tymin = std::numeric_limits<float>::min();
        tymax = std::numeric_limits<float>::max();
    }
    if (ray.direction.z == 0) {
        tzmin = std::numeric_limits<float>::min();
        tzmax = std::numeric_limits<float>::max();
    }

    float tinx = glm::min(txmin, txmax);
    float tiny = glm::min(tymin, tymax);
    float tinz = glm::min(tzmin, tzmax);
    float toutx = glm::max(txmin, txmax);
    float touty = glm::max(tymin, tymax);
    float toutz = glm::max(tzmin, tzmax);

    float tin = glm::max(glm::max(tinx, tiny), tinz);
    float tout = glm::min(glm::min(toutx, touty), toutz);

    if (tin > tout || tout < 0) {
        return false;
    }

    if (tin <= 0 && tout > 0 && tout < ray.t) {
        ray.t = tout;
        return true;
    }

    if (ray.t < tin) {
        return false;
    }
    ray.t = tin;
    return true;
}
