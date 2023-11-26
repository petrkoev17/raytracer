#include "texture.h"
#include <framework/image.h>
#include <cmath>

glm::vec3 acquireTexel(const Image& image, const glm::vec2& texCoord, const Features& features)
{
    // TODO: implement this function.
    // Given texcoords, return the corresponding pixel of the image
    // The pixel are stored in a 1D array of row major order
    // you can convert from position (i,j) to an index using the method seen in the lecture
    // Note, the center of the first pixel is at image coordinates (0.5, 0.5)
    int result = 0;
    if (features.enableTextureMapping) {

        int i = texCoord.x * image.width;
        int j = std::clamp((1.0f - texCoord.y) * image.height, 0.0f, image.width - 1.0f); // 1 - texCoord.y to invert texture

        result = j * image.width + i;

        return image.pixels[result];

    }
    return image.pixels[result];
}


glm::vec3 cubeMap(glm::vec3 direction) {
    
    glm::vec3 result { 0 };

    float absX = std::abs(direction.x);
    float absY = std::abs(direction.y);
    float absZ = std::abs(direction.z);

    float uc = 0; 
    float vc = 0;

    float faceIndex = 0;
    
    float maxAxis = 0;
    
    //Top face
    if (direction.y > 0 && absY >= absZ && absY >= absX) {
        maxAxis = absY;
        uc = direction.x;
        vc = -direction.z;
        faceIndex = 2;
    }

    //Bottom face
    if (direction.y < 0 && absY >= absZ && absY >= absX) {
        maxAxis = absY;
        uc = direction.x;
        vc = direction.z;
        faceIndex = 3;
    }

    //Front face
    if (direction.x > 0 && absX >= absZ && absX >= absY) {
        maxAxis = absX;
        uc = -direction.z;
        vc = direction.y;
        faceIndex = 0;
    }

    //Back face
    if (direction.x < 0 && absX >= absZ && absX >= absY) {
        maxAxis = absX;
        uc = direction.z;
        vc = direction.y;
        faceIndex = 1;
    }
    //Left face
    if (direction.z > 0 && absZ >= absY && absZ >= absX) {
        maxAxis = absZ;
        uc = direction.x;
        vc = direction.y;
        faceIndex = 4;
    }

    //Right face
    if (direction.z < 0 && absZ >= absY && absZ >= absX) {
        maxAxis = absZ;
        uc = -direction.x;
        vc = direction.y;
        faceIndex = 5;
    }

    //Convert range from (-1 ; 1) to (0 ; 1)
    uc = (1 / 2.0f) * (uc / maxAxis + 1.0f);
    vc = (1 / 2.0f) * (vc / maxAxis + 1.0f);

    result = glm::vec3 { uc, vc, faceIndex };
    return result;
}
