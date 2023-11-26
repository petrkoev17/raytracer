#include "light.h"
#include "config.h"
// Suppress warnings in third-party code.
#include <framework/disable_all_warnings.h>
DISABLE_WARNINGS_PUSH()
#include <glm/geometric.hpp>
DISABLE_WARNINGS_POP()
#include <cmath>
#include <iostream>
#include <random>

// samples a segment light source
// you should fill in the vectors position and color with the sampled position and color
void sampleSegmentLight(const SegmentLight& segmentLight, glm::vec3& position, glm::vec3& color)
{
    // Random sample gen
    std::random_device randomDev;
    std::mt19937 mtGen(randomDev());

    std::uniform_real_distribution<> generatorU(0.0, 1.0);
    float x = generatorU(mtGen);


    glm::vec3 direction = segmentLight.endpoint1 - segmentLight.endpoint0;

    glm::vec3 randomPos = direction * x;

    glm::vec3 randomSample = segmentLight.endpoint0 + randomPos;

    position = randomSample;

    glm::vec3 colorSample = x * segmentLight.color1 
        + (1 - x) * segmentLight.color0;

    color = colorSample;
}

// samples a parallelogram light source
// you should fill in the vectors position and color with the sampled position and color
void sampleParallelogramLight(const ParallelogramLight& parallelogramLight, glm::vec3& position, glm::vec3& color)
{
    // TODO: implement this function.
    std::random_device randomDevX;
    std::mt19937 mtGenX(randomDevX());

    std::uniform_real_distribution<> generatorU(0.0, 1.0);
    float x = generatorU(mtGenX);

    std::random_device randomDevY;
    std::mt19937 mtGenY(randomDevY());

    float y = generatorU(mtGenY);

    //Random positions from edges
    glm::vec3 randomPosX = parallelogramLight.edge01 * x;
    glm::vec3 randomPosY = parallelogramLight.edge02 * y;

    //Random Sample
    glm::vec3 randomSample = parallelogramLight.v0 + randomPosX + randomPosY;

    position = randomSample;


    //Random Color Sample
    glm::vec3 colorSample = (1 - y) * ((1 - x) * parallelogramLight.color0 + x * parallelogramLight.color1)
        + y * ((1 - x) * parallelogramLight.color2 + x * parallelogramLight.color3);

    color = colorSample;
}

// test the visibility at a given light sample
// returns 1.0 if sample is visible, 0.0 otherwise
float testVisibilityLightSample(const glm::vec3& samplePos, const glm::vec3& debugColor, const BvhInterface& bvh, const Features& features, Ray ray, HitInfo hitInfo)
{
    glm::vec3 point = ray.origin + ray.direction * ray.t;
    glm::vec3 shadowRayDirection = samplePos - point;
    Ray shadowRay = { point + 0.0001f*shadowRayDirection, shadowRayDirection }; 
    if (bvh.intersect(shadowRay, hitInfo, features)) {
        if (glm::distance(shadowRay.origin, shadowRay.origin + shadowRay.direction * shadowRay.t) < glm::distance(point, samplePos)) {
            if (shadowRay.direction.x != 0) {
                shadowRay.t = (samplePos.x - shadowRay.origin.x) / shadowRay.direction.x;
            }
            drawRay(shadowRay, { 1, 0, 0 });

            return 0;
        }
    }
    if (shadowRay.direction.x != 0) {
        shadowRay.t = (samplePos.x - shadowRay.origin.x) / shadowRay.direction.x;
    }
    drawRay(shadowRay, debugColor);
    return 1;
}

// given an intersection, computes the contribution from all light sources at the intersection point
// in this method you should cycle the light sources and for each one compute their contribution
// don't forget to check for visibility (shadows!)

// Lights are stored in a single array (scene.lights) where each item can be either a PointLight, SegmentLight or ParallelogramLight.
// You can check whether a light at index i is a PointLight using std::holds_alternative:
// std::holds_alternative<PointLight>(scene.lights[i])
//
// If it is indeed a point light, you can "convert" it to the correct type using std::get:
// PointLight pointLight = std::get<PointLight>(scene.lights[i]);
//
//
// The code to iterate over the lights thus looks like this:
// for (const auto& light : scene.lights) {
//     if (std::holds_alternative<PointLight>(light)) {
//         const PointLight pointLight = std::get<PointLight>(light);
//         // Perform your calculations for a point light.
//     } else if (std::holds_alternative<SegmentLight>(light)) {
//         const SegmentLight segmentLight = std::get<SegmentLight>(light);
//         // Perform your calculations for a segment light.
//     } else if (std::holds_alternative<ParallelogramLight>(light)) {
//         const ParallelogramLight parallelogramLight = std::get<ParallelogramLight>(light);
//         // Perform your calculations for a parallelogram light.
//     }
// }
//
// Regarding the soft shadows for **other** light sources **extra** feature:
// To add a new light source, define your new light struct in scene.h and modify the Scene struct (also in scene.h)
// by adding your new custom light type to the lights std::variant. For example:
// std::vector<std::variant<PointLight, SegmentLight, ParallelogramLight, MyCustomLightType>> lights;
//
// You can add the light sources programmatically by creating a custom scene (modify the Custom case in the
// loadScene function in scene.cpp). Custom lights will not be visible in rasterization view.
glm::vec3 computeLightContribution(const Scene& scene, const BvhInterface& bvh, const Features& features, Ray ray, HitInfo hitInfo)
{
    float hardContr = 1;
    int counter = 0;
    glm::vec3 result { 0 };
    glm::vec3 shadowContr { 1 };
    if (features.enableSoftShadow) {
        float f = 0;
        for (const auto& light : scene.lights) {
            if (std::holds_alternative<PointLight>(light)) {
                const PointLight pointLight = std::get<PointLight>(light);

            } else if (std::holds_alternative<SegmentLight>(light)) {
                const SegmentLight segmentLight = std::get<SegmentLight>(light);
                glm::vec3 position { 0 };
                glm::vec3 color { 0 };
                for (int i = 0; i < features.amountOfSamples; i++) {

                    sampleSegmentLight(segmentLight, position, color);
                    f += testVisibilityLightSample(position, color, bvh, features, ray, hitInfo);
                    counter++;
                }

            } else if (std::holds_alternative<ParallelogramLight>(light)) {
                const ParallelogramLight parallelogramLight = std::get<ParallelogramLight>(light);
                glm::vec3 position { 0 };
                glm::vec3 color { 0 };
                for (int i = 0; i < features.amountOfSamples; i++) {
                    sampleParallelogramLight(parallelogramLight, position, color);
                    f += testVisibilityLightSample(position, color, bvh, features, ray, hitInfo);
                    counter++;
                }
            }
        }

        shadowContr = glm::vec3 { f / counter };
    }

    if (features.enableHardShadow) {
        hardContr = 0;
        int counter = 0;
        for (const auto& light : scene.lights) {
            if (std::holds_alternative<PointLight>(light)) {
                const PointLight pointLight = std::get<PointLight>(light);
                // Perform your calculations for a point light.
                hardContr += testVisibilityLightSample(pointLight.position, pointLight.color, bvh, features, ray, hitInfo);
                counter++;
            }
        }
        hardContr /= counter;
    }
    if (features.enableShading) {
        // If shading is enabled, compute the contribution from all lights.
        int counter = 0;
        for (const auto& light : scene.lights) {
            if (std::holds_alternative<PointLight>(light)) {
                const PointLight pointLight = std::get<PointLight>(light);
                // Perform your calculations for a point light.
                counter++;
                result += computeShading(pointLight.position, pointLight.color, features, ray, hitInfo);

            } else if (std::holds_alternative<SegmentLight>(light)) {
                const SegmentLight segmentLight = std::get<SegmentLight>(light);
                // Perform your calculations for a segment light.
                glm::vec3 position { 0 };
                glm::vec3 color { 0 };
                glm::vec3 finalPos { 0 };
                glm::vec3 finalCol { 0 };
                for (int i = 0; i < features.amountOfSamples; i++) {
                    sampleSegmentLight(segmentLight, position, color);
                    finalPos += position;
                    finalCol += color;
                }

                finalCol /= features.amountOfSamples;
                finalPos /= features.amountOfSamples;
                counter++;
                result += computeShading(finalPos, finalCol, features, ray, hitInfo);

            } else if (std::holds_alternative<ParallelogramLight>(light)) {
                const ParallelogramLight parallelogramLight = std::get<ParallelogramLight>(light);
                // Perform your calculations for a parallelogram light.
                glm::vec3 position { 0 };
                glm::vec3 color { 0 };

                glm::vec3 finalPos { 0 };
                glm::vec3 finalCol { 0 };

                for (int i = 0; i < features.amountOfSamples; i++) {
                    sampleParallelogramLight(parallelogramLight, position, color);
                    finalPos += position;
                    finalCol += color;
                }

                finalPos /= features.amountOfSamples;
                finalCol /= features.amountOfSamples;
                counter++;
                result += computeShading(finalPos, finalCol, features, ray, hitInfo);
            }

        }
        result /= counter;
        // TODO: replace this by your own implementation of shading

        return result * shadowContr * hardContr;

    } else {
        // If shading is disabled, return the albedo of the material.
        return hitInfo.material.kd * hardContr;
    }
}