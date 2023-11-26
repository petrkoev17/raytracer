#include "bounding_volume_hierarchy.h"
#include "draw.h"
#include "intersect.h"
#include "scene.h"
#include "texture.h"
#include "interpolate.h"
#include <glm/glm.hpp>
#include "iostream"
#define TIMEIT



// Calculates the total number of Triangles in the scene to reserve the space for all the indices
int calculateTotalNumberOfTriangles(Scene& scene) {
    int total = 0;
    for (int i = 0; i < scene.meshes.size(); i++) {
        total += scene.meshes[i].triangles.size();
    }

    return total;
}

// Calculates the mesh the given triangle index belongs to and its corresponding index for the triangle vector
glm::vec2 calculateMeshAndIndex(int index, Scene& scene)
{
    int i;
    for (i = 0; i <= scene.meshes.size() - 1; i++) {
        int size = scene.meshes.at(i).triangles.size();
        index -= size;
        if (index < 0) {
            index += size;
            break;
        }
    }
    return glm::vec2 { i, index };
}

glm::vec2 calculateMeshAndIndex(int index, Scene* scene)
{
    int i;
    for (i = 0; i <= scene->meshes.size() - 1; i++) {
        int size = scene->meshes.at(i).triangles.size();
        index -= size;
        if (index < 0) {
            index += size;
            break;
        }
    }
    return glm::vec2 { i, index };
}

// Calculates the centroid position of the triangle
glm::vec3 calculateCentroid(int index, Scene& scene) {
    glm::vec2 meshAndIndex = calculateMeshAndIndex(index, scene);
    int meshIdx = meshAndIndex.x;
    int triangleIdx = meshAndIndex.y;

    glm::uvec3 vertexPos = scene.meshes[meshIdx].triangles[triangleIdx];
    glm::vec3 v0 = scene.meshes[meshIdx].vertices[vertexPos.x].position;
    glm::vec3 v1 = scene.meshes[meshIdx].vertices[vertexPos.y].position;
    glm::vec3 v2 = scene.meshes[meshIdx].vertices[vertexPos.z].position;

    return (v0 + v1 + v2) * (1 / 3.0f);
}

// Struct to sort 
struct sortIt {
    int axis;
    Scene& scene;
    bool inline operator()(const int left, const int right) {
        glm::vec3 cA = calculateCentroid(left, scene);
        glm::vec3 cB = calculateCentroid(right, scene);
        switch (axis) {
        case 0:
            return cA.x < cB.x;
        case 1:
            return cA.y < cB.y;
        case 2:
            return cA.z < cB.z;
        }
    }
};

// Function that calculates the AABB for the given list of absolute triangle Indices
AxisAlignedBox computeAABB(std::vector<int>& indices, Scene& scene) {
    float xmin = 1;
    float ymin = 1;
    float zmin = 1;
    float xmax = -1;
    float ymax = -1;
    float zmax = -1;
    for (int i = 0; i < indices.size(); i++) {
        glm::vec2 meshAndIndex = calculateMeshAndIndex(indices[i], scene);
        int meshIdx = meshAndIndex.x;
        int triangleIdx = meshAndIndex.y;

        glm::uvec3 vertexPos = scene.meshes[meshIdx].triangles[triangleIdx];
        std::vector<int> pos = { (int)vertexPos.x, (int)vertexPos.y, (int)vertexPos.z };
        for (int i = 0; i < 3; i++) {
            glm::vec3 v = scene.meshes[meshIdx].vertices[pos[i]].position;

            if (v.x < xmin)
                xmin = v.x;
            if (v.y < ymin)
                ymin = v.y;
            if (v.z < zmin)
                zmin = v.z;
            if (v.x > xmax)
                xmax = v.x;
            if (v.y > ymax)
                ymax = v.y;
            if (v.z > zmax)
                zmax = v.z;
        }        
    }

    return AxisAlignedBox {
        glm::vec3 { xmin, ymin, zmin }, glm::vec3 { xmax, ymax, zmax }
    };
}

// Helper function to create the tree recursively
void helper(int objectNumber, int level, int& numOfLeaves, int& numOfLevels, int maxLevel, Scene& scene, std::vector<Node>& nodes) {
    Node current = nodes[objectNumber];
    if (current.indices.size() > 1 && level < maxLevel) {
        std::vector<int> triangleIndices = current.indices;
        sort(triangleIndices.begin(), triangleIndices.end(), sortIt(level % 3, scene));

        std::vector<int> leftTriangles(triangleIndices.begin(), triangleIndices.begin() + triangleIndices.size() / 2);
        std::vector<int> rightTriangles(triangleIndices.begin() + triangleIndices.size() / 2, triangleIndices.end());

        int idx = nodes.size();

        nodes[objectNumber].isLeaf = false;
        nodes[objectNumber].indices = std::vector<int> { idx, idx + 1 };

        Node left { true, computeAABB(leftTriangles, scene), leftTriangles };
        Node right { true, computeAABB(rightTriangles, scene), rightTriangles };

        nodes.push_back(left);
        nodes.push_back(right);

        helper(idx, level + 1, numOfLeaves, numOfLevels, maxLevel, scene, nodes);
        helper(idx+1, level + 1, numOfLeaves, numOfLevels, maxLevel, scene, nodes);
    } else {
        numOfLeaves += 1;
        if (level > numOfLevels) {
            numOfLevels = level;
        }
        return;
    }
}

float calculateCompareValue(Node& currentNode, Sha& sha) {
    switch (sha.axis) {
    case 0:
        return currentNode.bounds.lower.x + sha.splitRatio * (currentNode.bounds.upper.x - currentNode.bounds.lower.x);
    case 1:
        return currentNode.bounds.lower.y + sha.splitRatio * (currentNode.bounds.upper.y - currentNode.bounds.lower.y);
    case 2:
        return currentNode.bounds.lower.z + sha.splitRatio * (currentNode.bounds.upper.z - currentNode.bounds.lower.z);
    }
}

bool isSmaller(float& compareValue, glm::vec3& centroid, int& axis) {
    switch (axis) {
    case 0:
        return centroid.x <= compareValue;
    case 1:
        return centroid.y <= compareValue;
    case 2:
        return centroid.z <= compareValue;
    }
}

float computeVolume(AxisAlignedBox& aabb) {
    glm::vec3 lengths = aabb.upper - aabb.lower;
    return lengths.x * lengths.y * lengths.z;
}

Split calculateSplit(Scene& scene, Node& currentNode, Sha sha) 
{
    std::vector<int> leftIndices = {};
    std::vector<int> rightIndices = {};
    float compareValue = calculateCompareValue(currentNode, sha);
    for (int i = 0; i < currentNode.indices.size(); i++) {
        glm::vec3 centroid = calculateCentroid(currentNode.indices[i], scene);
        if (isSmaller(compareValue, centroid, sha.axis)) {
            leftIndices.push_back(currentNode.indices[i]);
        } else {
            rightIndices.push_back(currentNode.indices[i]);
        }
    }

    AxisAlignedBox leftAABB = computeAABB(leftIndices, scene);
    AxisAlignedBox rightAABB = computeAABB(rightIndices, scene);

    float cost = computeVolume(leftAABB) * leftIndices.size() + computeVolume(rightAABB) * rightIndices.size();
    return Split(sha, leftIndices, rightIndices, leftAABB, rightAABB, cost);
}

Split findBestSplit(Scene& scene, Node& currentNode, std::vector<Sha>& splits) 
{
    std::vector<Split> options = {};
    options.reserve(9);
    for (int i = 0; i < splits.size(); i++) {
        Split current = calculateSplit(scene, currentNode, splits[i]);
        if (current.leftIndices.size() > 0 && current.rightIndices.size() > 0) {
            options.push_back(current);
        }
    }
    if (options.size() == 0) {
        return Split {
            Sha { 0, 0.25f }, std::vector<int> {}, std::vector<int> {}, AxisAlignedBox { glm::vec3 { 0.0f }, glm::vec3 { 0.0f } }, AxisAlignedBox { glm::vec3 { 0.0f }, glm::vec3 { 0.0f } }, -1.0f
        };
    } else {
        float minCost = options[0].cost;
        float index = 0;
        for (int i = 1; i < options.size(); i++) {
            if (options[i].cost < minCost) {
                minCost = options[i].cost;
                index = i;
            }
        }
        return options[index];
    }
}

void shaAndBinning(std::vector<Node>& nodes, int& maxLevel, int& numOfLeaves, int& numOfLevels, Scene& scene, 
    int objectNumber, std::vector<int>& levels, std::vector<Sha>& splits) 
{
    if (objectNumber >= nodes.size()) {
        return;
    }
    int level = levels[objectNumber];
    if (level < maxLevel && nodes[objectNumber].indices.size() >= 1) {
        Split split = findBestSplit(scene, nodes[objectNumber], splits);
        if (split.cost < 0) {
            numOfLeaves++;
            if (level > numOfLevels) {
                numOfLevels = level;
            }
            shaAndBinning(nodes, maxLevel, numOfLeaves, numOfLevels, scene, objectNumber + 1, levels, splits);
        } else {
            int indexOfLeftChild = nodes.size();
            nodes[objectNumber].isLeaf = false;
            nodes[objectNumber].indices = std::vector { indexOfLeftChild, indexOfLeftChild + 1 };

            Node left { true, split.leftBox, split.leftIndices };
            Node right { true, split.rightBox, split.rightIndices };
            nodes.push_back(left);
            nodes.push_back(right);
            levels.push_back(level + 1);
            levels.push_back(level + 1);
            shaAndBinning(nodes, maxLevel, numOfLeaves, numOfLevels, scene, objectNumber + 1, levels, splits);
        }
    } else {
        numOfLeaves++;
        if (level > numOfLevels) {
            numOfLevels = level;
        }
        shaAndBinning(nodes, maxLevel, numOfLeaves, numOfLevels, scene, objectNumber + 1, levels, splits);
    }
}

// Constructor
BoundingVolumeHierarchy::BoundingVolumeHierarchy(Scene* pScene, const Features& features)
    : m_pScene(pScene)
{
#ifdef TIMEIT
    auto t3 = std::chrono::high_resolution_clock::now();
#endif
    m_numLeaves = 0;
    nodes = {};
    m_numLevels = 0;
    maxLevel = 9;
    Scene scene = *pScene;

    // Values to calculate the AABB around the entire scene
    float xmin = 1;
    float ymin = 1;
    float zmin = 1;
    float xmax = -1;
    float ymax = -1;
    float zmax = -1;
    // counter for the indices of the Triangles
    int counter = 0;
    std::vector<int> triangleIndices {};
    triangleIndices.reserve(calculateTotalNumberOfTriangles(scene));
    for (int meshIdx = 0; meshIdx < scene.meshes.size(); meshIdx++) {
        Mesh mesh = scene.meshes[meshIdx];
        for (int triangleIdx = 0; triangleIdx < mesh.triangles.size(); triangleIdx++) {
            triangleIndices.push_back(counter);
            counter += 1;

            glm::uvec3 vertexIndicesVec = mesh.triangles[triangleIdx];
            std::vector<int> vertexIndices = std::vector<int> { (int)vertexIndicesVec.x, (int)vertexIndicesVec.y, (int)vertexIndicesVec.z };
            for (int i = 0; i < 3; i++) {
                glm::vec3 v = mesh.vertices[vertexIndices[i]].position;
                if (v.x < xmin)
                    xmin = v.x;
                if (v.y < ymin)
                    ymin = v.y;
                if (v.z < zmin)
                    zmin = v.z;
                if (v.x > xmax)
                    xmax = v.x;
                if (v.y > ymax)
                    ymax = v.y;
                if (v.z > zmax)
                    zmax = v.z;
            }
        }
    }
    AxisAlignedBox entireScene {
        glm::vec3 { xmin, ymin, zmin }, glm::vec3 { xmax, ymax, zmax }
    };

    Node root = { true, entireScene, triangleIndices };
    nodes.push_back(root);
    if (features.extra.enableBvhSahBinning) {
        std::vector<Sha> splits = {};
        splits.reserve(9);
        for (int i = 0; i < 3; i++) {
            splits.push_back(Sha { i, 0.25f });
            splits.push_back(Sha { i, 0.50f });
            splits.push_back(Sha { i, 0.75f });
        }
        std::vector<int> levels {};
        levels.push_back(0);
        shaAndBinning(nodes, maxLevel, m_numLeaves, m_numLevels, scene, 0, levels, splits);
    } else {
        helper(0, 0, m_numLeaves, m_numLevels, maxLevel, scene, nodes);
    }
#ifdef TIMEIT
    auto t4
        = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> delta2 = t4 - t3;
    std::cout << "Creation time: " << delta2.count() << " ms" << std::endl;
#endif
}

// Return the depth of the tree that you constructed. This is used to tell the
// slider in the UI how many steps it should display for Visual Debug 1.
int BoundingVolumeHierarchy::numLevels() const
{
    return m_numLevels;
}

// Return the number of leaf nodes in the tree that you constructed. This is used to tell the
// slider in the UI how many steps it should display for Visual Debug 2.
int BoundingVolumeHierarchy::numLeaves() const
{
    return m_numLeaves;
}

// Return the number of internal nodes in the tree
int BoundingVolumeHierarchy::numOfInternalNodes() const 
{
    return nodes.size() - m_numLeaves;
}

// Gets the indices of the nodes of the level 

void calcLevelIndices(std::vector<Node>& levelNodes, int currentLvl, int index, int level, std::vector<Node>& nodes) {
    Node current = nodes[index];
    if (level == currentLvl) {
        levelNodes.push_back(current);
    } else if (current.isLeaf) {
        return;
    } else {
        calcLevelIndices(levelNodes, currentLvl + 1, current.indices[0], level, nodes);
        calcLevelIndices(levelNodes, currentLvl + 1, current.indices[1], level, nodes);
    }
}

// Use this function to visualize your BVH. This is useful for debugging. Use the functions in
// draw.h to draw the various shapes. We have extended the AABB draw functions to support wireframe
// mode, arbitrary colors and transparency.
void BoundingVolumeHierarchy::debugDrawLevel(int level)
{
    if (level < 0) {
        level = 0;
    }
    // Draw the AABB as a transparent green box.
    //AxisAlignedBox aabb{ glm::vec3(-0.05f), glm::vec3(0.05f, 1.05f, 1.05f) };
    //drawShape(aabb, DrawMode::Filled, glm::vec3(0.0f, 1.0f, 0.0f), 0.2f);

    // Draw the AABB as a (white) wireframe box.
    //AxisAlignedBox aabb { glm::vec3(0.0f), glm::vec3(0.0f, 1.05f, 1.05f) };
    //drawAABB(aabb, DrawMode::Wireframe);
    //drawAABB(aabb, DrawMode::Filled, glm::vec3(0.05f, 1.0f, 0.05f), 0.1f)
    if (level == 5) {
        int a = 4;
    }
    std::vector<Node> levelNodes = {};
    calcLevelIndices(levelNodes, 0, 0, level, nodes);
    for (int i = 0; i < levelNodes.size(); i++) {
        Node current = levelNodes[i];
        drawAABB(current.bounds, DrawMode::Wireframe);
    }
    
}


// Use this function to visualize your leaf nodes. This is useful for debugging. The function
// receives the leaf node to be draw (think of the ith leaf node). Draw the AABB of the leaf node and all contained triangles.
// You can draw the triangles with different colors. NoteL leafIdx is not the index in the node vector, it is the
// i-th leaf node in the vector.
void BoundingVolumeHierarchy::debugDrawLeaf(int leafIdx)
{
    // Draw the AABB as a transparent green box.
    //AxisAlignedBox aabb{ glm::vec3(-0.05f), glm::vec3(0.05f, 1.05f, 1.05f) };
    //drawShape(aabb, DrawMode::Filled, glm::vec3(0.0f, 1.0f, 0.0f), 0.2f);

    // Draw the AABB as a (white) wireframe box.
    //AxisAlignedBox aabb { glm::vec3(0.0f), glm::vec3(0.0f, 1.05f, 1.05f) };
    //drawAABB(aabb, DrawMode::Wireframe);
    //drawAABB(aabb, DrawMode::Filled, glm::vec3(0.05f, 1.0f, 0.05f), 0.1f);

    // once you find the leaf node, you can use the function drawTriangle (from draw.h) to draw the contained primitives
    Scene scene = *m_pScene;
    int counter = 0;
    Node leaf;
    for (int i = 0; i < nodes.size(); i++) {
        Node current = nodes[i];
        if (current.isLeaf && counter == leafIdx-1) {
            leaf = current;
            break;
        } else if (current.isLeaf) {
            counter++;
        }
    }

    drawAABB(leaf.bounds, DrawMode::Wireframe);
    for (int i = 0; i < leaf.indices.size(); i++) {
        glm::vec2 meshAndIndex = calculateMeshAndIndex(leaf.indices[i], scene);
        int meshIdx = meshAndIndex.x;
        int triangleIdx = meshAndIndex.y;
        glm::uvec3 vertexPos = scene.meshes[meshIdx].triangles[triangleIdx];
        Vertex v0 = scene.meshes[meshIdx].vertices[vertexPos.x];
        Vertex v1 = scene.meshes[meshIdx].vertices[vertexPos.y];
        Vertex v2 = scene.meshes[meshIdx].vertices[vertexPos.z];
        drawTriangle(v0, v1, v2);

    }
}


void getIndices(Ray& ray, int index, const std::vector<Node>& nodes, std::vector<std::pair<int, float>>& intersections)
{
    Node node = nodes[index];



    float oldt = ray.t;
    ray.t = std::numeric_limits<float>::max();

    if (node.isLeaf) {
        intersections.push_back(std::make_pair(index, oldt));
        return;
    }

    if (intersectRayWithShape(nodes[node.indices[0]].bounds, ray)) {

        getIndices(ray, node.indices[0], nodes, intersections);
    }
    ray.t = std::numeric_limits<float>::max();

    if (intersectRayWithShape(nodes[node.indices[1]].bounds, ray)) {

        getIndices(ray, node.indices[1], nodes, intersections);
    }
    
}

Node findInternalNode(std::vector<Node>& nodes, int& index) {
    int currentIndex = 0;
    int objectNumber = 0;
    for (int i = 0; i < nodes.size(); i++) {
        if (!nodes[i].isLeaf && currentIndex == index) {
            objectNumber = i;
            break;
        } else if (!nodes[i].isLeaf) {
            currentIndex += 1;
        }
    }
    return nodes[objectNumber];
}


// Draws the AABB of the internalNode and its corresponding children
void BoundingVolumeHierarchy::drawSplit(int internalNodeIdx) 
{
    std::vector<Node> list = nodes;
    Node internal = findInternalNode(list, internalNodeIdx);
    Node leftChild = nodes[internal.indices[0]];
    Node rightChild = nodes[internal.indices[1]];
    drawAABB(internal.bounds, DrawMode::Wireframe);
    drawAABB(leftChild.bounds, DrawMode::Wireframe, glm::vec3 { 1.0f, 0.0f, 0.0f });
    drawAABB(rightChild.bounds, DrawMode::Wireframe, glm::vec3 { 0.0f, 1.0f, 0.0f });
}

float gettout(const AxisAlignedBox& box, Ray& ray)
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


    float toutx = glm::max(txmin, txmax);
    float touty = glm::max(tymin, tymax);
    float toutz = glm::max(tzmin, tzmax);

    float tout = glm::min(glm::min(toutx, touty), toutz);
    return tout;

}

    // Return true if something is hit, returns false otherwise. Only find hits if they are closer than t stored
// in the ray and if the intersection is on the correct side of the origin (the new t >= 0). Replace the code
// by a bounding volume hierarchy acceleration structure as described in the assignment. You can change any
// file you like, including bounding_volume_hierarchy.h.
bool BoundingVolumeHierarchy::intersect(Ray& ray, HitInfo& hitInfo, const Features& features) const
{
    Vertex final0;
    Vertex final1;
    Vertex final2;
    // If BVH is not enabled, use the naive implementation.
    if (!features.enableAccelStructure) {
        bool hit = false;
        // Intersect with all triangles of all meshes.
        for (const auto& mesh : m_pScene->meshes) {
            for (const auto& tri : mesh.triangles) {
                const auto v0 = mesh.vertices[tri[0]];
                const auto v1 = mesh.vertices[tri[1]];
                const auto v2 = mesh.vertices[tri[2]];
                if (intersectRayWithTriangle(v0.position, v1.position, v2.position, ray, hitInfo)) {
                    hitInfo.material = mesh.material;
                    
                    hitInfo.barycentricCoord = computeBarycentricCoord(v0.position, v1.position, v2.position, ray.origin + ray.direction * ray.t);
                    hitInfo.normal = interpolateNormal(v0.normal, v1.normal, v2.normal, hitInfo.barycentricCoord);
                    hitInfo.texCoord = interpolateTexCoord(v0.texCoord, v1.texCoord, v2.texCoord, hitInfo.barycentricCoord);
                    
                    hit = true;
                    final0 = v0;
                    final1 = v1;
                    final2 = v2;
                }
            }
        }
        // Intersect with spheres.
        for (const auto& sphere : m_pScene->spheres)
            hit |= intersectRayWithShape(sphere, ray, hitInfo);
        if (features.enableNormalInterp) {
            hitInfo.barycentricCoord = computeBarycentricCoord(final0.position, final1.position, final2.position, ray.origin + (ray.direction * ray.t));
            hitInfo.normal = interpolateNormal(final0.normal, final1.normal, final2.normal, hitInfo.barycentricCoord);
        } 
        if (features.debugFlags.debugNormalInterpolation) {
            drawNormals(final0, final1, final2, ray.origin + (ray.direction * ray.t), hitInfo.normal);
        }
        return hit;
    } else {
        // TODO: implement here the bounding volume hierarchy traversal.
        // Please note that you should use `features.enableNormalInterp` and `features.enableTextureMapping`
        // to isolate the code that is only needed for the normal interpolation and texture mapping features.

        bool hit = false;
        bool flag = false;

        int counter = 0;

        std::vector<std::pair<int, float>> intersections;
        if (!intersectRayWithShape(nodes[0].bounds, ray)) {
            return false;
        }
        getIndices(ray, 0, nodes, intersections);

        std::sort(intersections.begin(), intersections.end(), [](const std::pair<int, float>& lhs, const std::pair<int, float>& rhs) { return lhs.second < rhs.second; });

        Vertex vert0;
        Vertex vert1;
        Vertex vert2;
        float tout = std::numeric_limits<float>::max();
         for (std::pair<int, float> nodeIndex : intersections) {
             
            if (hit && tout < nodeIndex.second) {
                if (features.showVisitedBoxes) {
                    drawAABB(nodes[nodeIndex.first].bounds, DrawMode::Wireframe, glm::vec3 { 0.0f, 1.0f, 0.0 });
                }
                continue;
            }
            drawAABB(nodes[nodeIndex.first].bounds, DrawMode::Wireframe);

            for (int index : nodes[nodeIndex.first].indices) {             

                 glm::vec2 meshAndIndex = calculateMeshAndIndex(index, m_pScene);
                 int meshIdx = meshAndIndex.x;
                 int triangleIdx = meshAndIndex.y;
                 
                 glm::uvec3& vertexPos = m_pScene->meshes[meshIdx].triangles[triangleIdx];
                 auto& v0 = m_pScene->meshes[meshIdx].vertices[vertexPos.x];
                 auto& v1 = m_pScene->meshes[meshIdx].vertices[vertexPos.y];
                 auto& v2 = m_pScene->meshes[meshIdx].vertices[vertexPos.z];


                 if (intersectRayWithTriangle(v0.position, v1.position, v2.position, ray, hitInfo)) {
                     hitInfo.material = m_pScene->meshes[meshIdx].material;
                     hitInfo.barycentricCoord = computeBarycentricCoord(v0.position, v1.position, v2.position, ray.origin + ray.direction * ray.t);
                     hitInfo.normal = interpolateNormal(v0.normal, v1.normal, v2.normal, hitInfo.barycentricCoord);
                     hitInfo.texCoord = interpolateTexCoord(v0.texCoord, v1.texCoord, v2.texCoord, hitInfo.barycentricCoord);
                     vert0 = m_pScene->meshes[meshIdx].vertices[vertexPos.x];
                     vert1 = m_pScene->meshes[meshIdx].vertices[vertexPos.y];
                     vert2 = m_pScene->meshes[meshIdx].vertices[vertexPos.z];
                     hit = true;
                     flag = true;
                 }
             }
            if (flag) {
                counter++;
            }
            if (counter == 1) {
                tout = gettout(nodes[nodeIndex.first].bounds, ray);
            }
           
         }
         if (hit) {
             drawTriangle(vert0, vert1, vert2);
         }
         if (features.enableNormalInterp) {
             hitInfo.barycentricCoord = computeBarycentricCoord(vert0.position, vert1.position, vert2.position, ray.origin + (ray.direction * ray.t));
             hitInfo.normal = interpolateNormal(vert0.normal, vert1.normal, vert2.normal, hitInfo.barycentricCoord);
         }
         if (features.debugFlags.debugNormalInterpolation) {
             drawNormals(vert0, vert1, vert2, ray.origin + (ray.direction * ray.t), hitInfo.normal);
         }

         return hit;
         
      
    }
}