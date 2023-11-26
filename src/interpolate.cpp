#include "interpolate.h"
#include <glm/geometric.hpp>
#include "scene.h"

glm::vec3 computeBarycentricCoord (const glm::vec3& v0, const glm::vec3& v1, const glm::vec3& v2, const glm::vec3& p)
{
    float totalArea = glm::length(glm::cross(v2 - v0, v1 - v0)) / 2;

    float alpha = (glm::length(glm::cross(v2 - p, v1 - p)) / 2) / totalArea;
    float beta = (glm::length(glm::cross(v0 - p, v2 - p)) / 2) / totalArea;
    float gamma = (glm::length(glm::cross(v1 - p, v0 - p)) / 2) / totalArea;

    return glm::vec3 { alpha, beta, gamma };
    // TODO: implement this function.
}

glm::vec3 interpolateNormal (const glm::vec3& n0, const glm::vec3& n1, const glm::vec3& n2, const glm::vec3 barycentricCoord)
{
    float alpha = barycentricCoord.x;
    float beta = barycentricCoord.y;
    float gamma = barycentricCoord.z;

    glm::vec3 normal = alpha * n0 + beta * n1 + gamma * n2;
    return normal;
}


glm::vec2 interpolateTexCoord (const glm::vec2& t0, const glm::vec2& t1, const glm::vec2& t2, const glm::vec3 barycentricCoord)
{
// TODO: implement this function.
    glm::vec2 result { 0, 0 };
    result = t0 * barycentricCoord.x + t1 * barycentricCoord.y + t2 * barycentricCoord.z;

    return result;
}
