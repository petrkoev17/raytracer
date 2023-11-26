#include "texture.h"
#include <cmath>
#include <glm/geometric.hpp>
#include <shading.h>

const glm::vec3 computeShading(const glm::vec3& lightPosition, const glm::vec3& lightColor, const Features& features, Ray ray, HitInfo hitInfo)
{
    // TODO: implement the Phong shading model.

    glm::vec3 result = glm::vec3 { 0, 0, 0 };

    

    //Diffuse
    glm::vec3 lightRay = normalize(lightPosition - (ray.origin + ray.direction * ray.t));

    float angle = glm::dot(lightRay, normalize(hitInfo.normal));
    angle = glm::clamp(angle, 0.0f, 1.0f);


    glm::vec3 diffuse = lightColor * angle * hitInfo.material.kd;


    //Specular
    glm::vec3 cameraRayS = normalize(ray.origin - (ray.origin + ray.direction * ray.t));

    glm::vec3 reflectedRay = glm::reflect(lightRay, hitInfo.normal);
    float specVal = glm::dot(reflectedRay, cameraRayS);
    specVal = glm::clamp(specVal, 0.0f, 1.0f);
    
    glm::vec3 specular = lightColor * hitInfo.material.ks * glm::pow(specVal, hitInfo.material.shininess);

    if (angle <= 0) {
        specular = glm::vec3 { 0 };
    }

    //Final Shading
    result = diffuse + specular;

    return result;
}


const Ray computeReflectionRay (Ray ray, HitInfo hitInfo)
{
    // Do NOT use glm::reflect!! write your own code.
    Ray reflectionRay {};

    glm::vec3 normal = glm::normalize(hitInfo.normal);
    reflectionRay.direction = ray.direction - 2.0f * (glm::dot(ray.direction, normal) * normal);
    reflectionRay.origin = ray.origin + ray.direction * ray.t;
    // TODO: implement the reflection ray computation.
    return reflectionRay;
}