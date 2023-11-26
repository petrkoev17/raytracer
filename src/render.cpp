#include "render.h"
#include "intersect.h"
#include "light.h"
#include "screen.h"
#include "texture.h"
#include <framework/trackball.h>
#ifdef NDEBUG
#include <omp.h>
#endif
#include <iostream>

#include <cmath>
#include <iostream>
#include <random>

Square getSquare(glm::vec3& r)
{
    glm::vec3 u;
    if (r != glm::vec3 { 1, 0, 0 }) {
        u = glm::normalize(glm::cross(r, glm::vec3 { 1, 0, 0 }));
    } else {
        u = glm::normalize(glm::cross(r, glm::vec3 { 0, 1, 0 }));
    }
    glm::vec3 v = glm::normalize(glm::cross(u, r));

    return Square { r, u, v };
}

glm::vec3 getGlossyReflection(Ray ray, const BvhInterface& bvh, const Scene& scene, const Features& features)
{
    HitInfo hitinfo;



    if (bvh.intersect(ray, hitinfo, features) && hitinfo.material.ks != glm::vec3 {0.0f}) {
        Ray reflection = computeReflectionRay(ray, hitinfo);
        reflection.origin += 0.0001f * reflection.direction;
        HitInfo reflectionhit;
        glm::vec3 sumOfColors { 0.0f };
        bvh.intersect(reflection, reflectionhit, features);
        glm::vec3 colorReflection = getFinalColor(scene, bvh, reflection, features);
        sumOfColors += colorReflection;
        drawRay(reflection, colorReflection);
        // If shininess == 1, the reflection ray is sufficient
        if (hitinfo.material.shininess > 1) {
            glm::vec3 r = glm::normalize(reflection.direction);
            Square square = getSquare(r);
            std::random_device randomDev;
            std::mt19937 mtGen(randomDev());

            std::uniform_real_distribution<> generatorU(0.0, 1.0);
            float a = 0.20f * glm::min(1.0f, 1 / hitinfo.material.shininess);
            for (int i = 0; i < 15; i++) {
                float x = generatorU(mtGen);
                float y = generatorU(mtGen);

                float u1 = (-a / 2.0f) + x * a;
                float v1 = (-a / 2.0f) + y * a;

                glm::vec3 direction = square.r + u1 * square.u + v1 * square.v;
                Ray sampledRay { reflection.origin, direction };
                HitInfo sample {};
                bvh.intersect(sampledRay, sample, features);
                glm::vec3 colorSample = getFinalColor(scene, bvh, sampledRay, features);
                drawRay(sampledRay, colorSample);
                sumOfColors += colorSample;
            }

            return sumOfColors / 16.0f;
        } else {
            return sumOfColors;
        }
    }

    return glm::vec3 { 0.0f };
}

int max = 2;
int min = -2;
float RandomFloat(float min, float max)
{
    float random = ((float)rand()) / (float)RAND_MAX;
    float range = max - min;
    return (random * range) + min;
}

glm::vec3 getFinalColor(const Scene& scene, const BvhInterface& bvh, Ray ray, const Features& features, int rayDepth)
{
    HitInfo hitInfo;
    if (features.extra.enableMultipleRaysPerPixel && rayDepth == 0 && features.amountOfRays != 1) {
        glm::vec3 Lo = glm::vec3(0);
        for (int i = 0; i < features.amountOfRays; i++) {
            glm::vec3 randomVec = { RandomFloat(0.99, 1.01), RandomFloat(0.99, 1.01), RandomFloat(0.99, 1.01) };
            Ray randomRay = { ray.origin * randomVec, ray.direction };
            Lo += getFinalColor(scene, bvh, randomRay, features, 1);
        }
        Lo /= features.amountOfRays;
        return Lo;
    }

    if (bvh.intersect(ray, hitInfo, features)) {

        glm::vec3 Lo = computeLightContribution(scene, bvh, features, ray, hitInfo);

        if (features.extra.enableMotionBlur && enableDebugDraw) {
            glm::vec3 origin = ray.origin;
            Ray debugRay = ray;
            HitInfo old = hitInfo;
            for (int i = min; i <= max; i++) {
                debugRay.t = std::numeric_limits<float>::max();
                debugRay.origin = origin + (float) i * features.motion / 10.0f;
                bvh.intersect(debugRay, hitInfo, features);
                drawRay(debugRay, {1, 1, 1});
            }
            hitInfo = old;
        }

        if (features.enableRecursive) {
            Ray reflection = computeReflectionRay(ray, hitInfo);
            if ((hitInfo.material.ks.x != 0 || hitInfo.material.ks.y != 0 || hitInfo.material.ks.z != 0) && rayDepth < 3) {
                reflection.origin += 0.00001f * reflection.direction;

                glm::vec3 newColour = getFinalColor(scene, bvh, reflection, features, ++rayDepth);
                Lo += hitInfo.material.ks * newColour;
            }
        }

         if (features.extra.enableEnvironmentMapping && features.enableRecursive && hitInfo.material.kdTexture) {

            Ray reflected = computeReflectionRay(ray, hitInfo);
            glm::vec3 envMapRefl = cubeMap(normalize(reflected.direction));
            glm::vec2 texCoord = { envMapRefl.x, envMapRefl.y };
            float faceIndex = envMapRefl.z;

            glm::vec3 pixel = acquireTexel(scene.cubeMap[faceIndex], texCoord, features);
            return pixel;
        }




        if (features.enableTextureMapping) {
            if (hitInfo.material.kdTexture) {
                glm::vec3 result = (acquireTexel(*hitInfo.material.kdTexture.get(), hitInfo.texCoord, features)); 
                return result;
            }

        }

        // Draw a white debug ray if the ray hits.
        drawRay(ray, Lo);

        // Set the color of the pixel to white if the ray hits.
        return Lo;
    } else {
        if (features.extra.enableEnvironmentMapping) {
            glm::vec3 uvFaceIndex = cubeMap(ray.direction);

            glm::vec2 texCoord = { uvFaceIndex.x, uvFaceIndex.y };
            float faceIndex = uvFaceIndex.z;

            glm::vec3 pixel = acquireTexel(scene.cubeMap[faceIndex], texCoord, features);
            drawRay(ray, pixel);
            return pixel;
        }

       
        // Draw a red debug ray if the ray missed.
        drawRay(ray, glm::vec3(1.0f, 0.0f, 0.0f));
        // Set the color of the pixel to black if the ray misses.
        return glm::vec3(0.0f);
    }
}

void renderRayTracing(const Scene& scene, const Trackball& camera, const BvhInterface& bvh, Screen& screen, const Features& features)
{
    glm::ivec2 windowResolution = screen.resolution();
    // Enable multi threading in Release mode
#ifdef NDEBUG
#pragma omp parallel for schedule(guided)
#endif
    for (int y = 0; y < windowResolution.y; y++) {
        for (int x = 0; x != windowResolution.x; x++) {
            // NOTE: (-1, -1) at the bottom left of the screen, (+1, +1) at the top right of the screen.
            const glm::vec2 normalizedPixelPos {
                float(x) / float(windowResolution.x) * 2.0f - 1.0f,
                float(y) / float(windowResolution.y) * 2.0f - 1.0f
            };

            if (features.extra.enableMotionBlur) {
                int amountOfRays = 5;
                glm::vec3 r1 = glm::vec3(0);
                for (int i = 0; i < amountOfRays; i++) {
                    int randNum = rand() % (max - min + 1) + min;
                    glm::vec3 change = features.motion * (float)randNum / 30.0f;

                    Ray cameraRay = camera.generateRay(normalizedPixelPos, features);
                    cameraRay.origin += change;

                    r1 += getFinalColor(scene, bvh, cameraRay, features);
                }
                r1 /= amountOfRays;
                if (features.extra.enableGlossyReflection) {
                    r1 += getGlossyReflection(camera.generateRay(normalizedPixelPos, features), bvh, scene, features);
                }
                screen.setPixel(x, y, r1);
                continue;
            }
            const Ray cameraRay = camera.generateRay(normalizedPixelPos, features);
            if (features.extra.enableGlossyReflection) {
                screen.setPixel(x, y, getFinalColor(scene, bvh, cameraRay, features)+ getGlossyReflection(cameraRay, bvh, scene, features));
            } else {
                screen.setPixel(x, y, getFinalColor(scene, bvh, cameraRay, features));
            }
            
        }
    }
    if (features.extra.enableMotionBlur ) {
        std::vector<glm::vec3> pixels = screen.pixels();
        for (int y = 2; y < windowResolution.y - 2; y++) {
            for (int x = 2; x != windowResolution.x - 2; x++) {
                glm::vec3 newcol = glm::vec3(0);
                int counter = 0;
                for (int y1 = -2; y1 < 3; y1++) {
                    for (int x1 = -2; x1 < 3; x1++) {
                        int n = screen.indexAt(x + x1, y + y1);
                        counter++;
                        newcol += pixels[n];
                    }
                }
                newcol /= counter;
                screen.setPixel(x, y, newcol);
            }
        }
    }
}