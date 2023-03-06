//
//  Mesh.cpp
//  raytracer_ext
//
//  Created by Sofia Iannicelli on 3/4/23.
//

#include "Mesh.hpp"
#include "Triangle.hpp"
#include <vector>
#include <Eigen/Dense>

using Eigen::MatrixXd;

void Mesh::computeBBInfo() {
    bb_xmin = vertices.at(0).x();
    bb_xmax = vertices.at(0).x();
    bb_ymin = vertices.at(0).y();
    bb_ymax = vertices.at(0).y();
    bb_zmin = vertices.at(0).z();
    bb_zmax = vertices.at(0).z();
    for(vec3<double> vertex : vertices) {
        if(vertex.x() < bb_xmin) {
            bb_xmin = vertex.x();
        }
        else if(vertex.x() > bb_xmax) {
            bb_xmax = vertex.x();
        }
        
        if(vertex.y() < bb_ymin) {
            bb_ymin = vertex.y();
        }
        else if(vertex.y() > bb_ymax) {
            bb_ymax = vertex.y();
        }
        
        if(vertex.z() < bb_zmin) {
            bb_zmin = vertex.z();
        }
        else if(vertex.z() > bb_zmax) {
            bb_zmax = vertex.z();
        }
    }
}

bool Mesh::fallsInBoundingBox(vec3<double> point) {
    double x = point.x();
    double y = point.y();
    double z = point.z();
    
    if(x < bb_xmin || x > bb_xmax) {
        return false;
    }
    if(y < bb_ymin || y > bb_ymax) {
        return false;
    }
    if(z < bb_zmin || z > bb_zmax) {
        return false;
    }
    return true;
}

double Mesh::findRayRectangleIntersection(std::vector<vec3<double>> vertices, Ray ray) {
    // assuming the vertices are specified in CCW order
    vec3<double> vector1 = getUnitVector(vertices.at(0) - vertices.at(1));
    vec3<double> vector2 = getUnitVector(vertices.at(2) - vertices.at(1));
    
    vec3<double> normal = getUnitVector(cross(vector2, vector1));
    
    // check if the ray intersects the plane containing the rectangle
    // information about the plane equation
    double a = normal[0];
    double b = normal[1];
    double c = normal[2];
    
    vec3<double> vertex = vertices.at(0);
    double d = -a * vertex.x() - b * vertex.y() - c * vertex.z();
    
    // extract out the ray info
    vec3<double> ray_o = ray.getOrigin();
    vec3<double> ray_d = ray.getDirection();
    
    double denominator = (a * ray_d.x() + b * ray_d.y() + c * ray_d.z());
    if(denominator == 0) { return -1; }
    double t = -(a * ray_o.x() + b * ray_o.y() + c * ray_o.z() + d) / denominator;
    
    if(t < 0) {
        return t;
    }
    
    vec3<double> intersectionPt = ray.getPointOnRay(t);
    // check if the intersection point is inside the rectangle
    for(int i = 0; i < 4; ++i) {
        vec3<double> testVertex = vertices.at(i);
        vec3<double> nextVertex = vertices.at((i + 1) % 3);
        
        vec3<double> testNormal = cross((nextVertex - testVertex), (intersectionPt - testVertex));
        
        if(dot(normal, testNormal) < 0) {
            return -1;
        }
    }
    
    return t;
}

double Mesh::findRayBBIntersection(Ray ray) {
    // bounding box is 6 rectangles
    std::vector<vec3<double>> frontVertices;
    frontVertices.push_back(vec3<double>(bb_xmin, bb_ymin, bb_zmax));
    frontVertices.push_back(vec3<double>(bb_xmax, bb_ymin, bb_zmax));
    frontVertices.push_back(vec3<double>(bb_xmax, bb_ymax, bb_zmax));
    frontVertices.push_back(vec3<double>(bb_xmin, bb_ymin, bb_zmax));
    double frontResult = findRayRectangleIntersection(frontVertices, ray);
    if(frontResult > 0) {
        return frontResult;
    }
    
    std::vector<vec3<double>> backVertices;
    backVertices.push_back(vec3<double>(bb_xmin, bb_ymin, bb_zmin));
    backVertices.push_back(vec3<double>(bb_xmax, bb_ymin, bb_zmin));
    backVertices.push_back(vec3<double>(bb_xmax, bb_ymax, bb_zmin));
    backVertices.push_back(vec3<double>(bb_xmin, bb_ymin, bb_zmin));
    double backResult = findRayRectangleIntersection(backVertices, ray);
    if(backResult > 0) {
        return backResult;
    }
    
    std::vector<vec3<double>> leftVertices;
    leftVertices.push_back(vec3<double>(bb_xmin, bb_ymin, bb_zmin));
    leftVertices.push_back(vec3<double>(bb_xmin, bb_ymin, bb_zmax));
    leftVertices.push_back(vec3<double>(bb_xmin, bb_ymax, bb_zmax));
    leftVertices.push_back(vec3<double>(bb_xmin, bb_ymax, bb_zmin));
    double leftResult = findRayRectangleIntersection(leftVertices, ray);
    if(leftResult > 0) {
        return leftResult;
    }
    
    std::vector<vec3<double>> rightVertices;
    rightVertices.push_back(vec3<double>(bb_xmax, bb_ymin, bb_zmin));
    rightVertices.push_back(vec3<double>(bb_xmax, bb_ymin, bb_zmax));
    rightVertices.push_back(vec3<double>(bb_xmax, bb_ymax, bb_zmax));
    rightVertices.push_back(vec3<double>(bb_xmax, bb_ymax, bb_zmin));
    double rightResult = findRayRectangleIntersection(rightVertices, ray);
    if(rightResult > 0) {
        return rightResult;
    }
    
    std::vector<vec3<double>> topVertices;
    topVertices.push_back(vec3<double>(bb_xmin, bb_ymax, bb_zmax));
    topVertices.push_back(vec3<double>(bb_xmax, bb_ymax, bb_zmax));
    topVertices.push_back(vec3<double>(bb_xmax, bb_ymax, bb_zmin));
    topVertices.push_back(vec3<double>(bb_xmin, bb_ymax, bb_zmin));
    double topResult = findRayRectangleIntersection(topVertices, ray);
    if(topResult > 0) {
        return topResult;
    }
    
    std::vector<vec3<double>> bottomVertices;
    bottomVertices.push_back(vec3<double>(bb_xmin, bb_ymin, bb_zmax));
    bottomVertices.push_back(vec3<double>(bb_xmax, bb_ymin, bb_zmax));
    bottomVertices.push_back(vec3<double>(bb_xmax, bb_ymin, bb_zmin));
    bottomVertices.push_back(vec3<double>(bb_xmin, bb_ymin, bb_zmin));
    double bottomResult = findRayRectangleIntersection(bottomVertices, ray);
    if(bottomResult > 0) {
        return bottomResult;
    }
    
    return -1.0;
}

double Mesh::findRayObjectIntersection(Ray ray) {
    // check if intersects with bounding box
    if(findRayBBIntersection(ray) < 0) {
        return -1.0;
    }
    
    // else check all triangles
    for(int i = 0; i < vertices.size(); ++i) {
        std::vector<vec3<double>> triangle_vertices;
        triangle_vertices.push_back(vertices.at(i));
        triangle_vertices.push_back(vertices.at((i+1)%vertices.size()));
        triangle_vertices.push_back(vertices.at((i+2)%vertices.size()));
        
        Triangle triangle = Triangle(triangle_vertices);
        double testT = triangle.findRayObjectIntersection(ray);
        if(testT > 0) {
            // triangleIndex = i;
            mostRecentlyHitTriangleIndex = i;
            return testT;
        }
    }
    
    return -1.0;
}

int mod(const int x, const int y)
{
   int t = x - ((x / y) * y);
   if (t < 0) t += y;
   return t;
}

void Mesh::computeVertexNormals() {
    vec3<double> planeNormalSum = vec3<double>(0,0,0);
    // each vertex
    for(int i = 0; i < vertices.size(); ++i) {
        // each triangle connected to vertex
        for(int j = -2; j <= 0; ++j) {
            int idx0 = mod((i + j), (int)vertices.size());
            int idx1 = mod((i + j + 1), (int)vertices.size());
            int idx2 = mod((i + j + 2), (int)vertices.size());
            std::vector<vec3<double>> triangle_vertices;
            triangle_vertices.push_back(vertices.at(idx0));
            triangle_vertices.push_back(vertices.at(idx1));
            triangle_vertices.push_back(vertices.at(idx2));
            
            Triangle triangle = Triangle(triangle_vertices);
            planeNormalSum += triangle.getPlaneNormal();
        }
        vec3<double> vectorNormal = (1.0/3.0) * planeNormalSum;
        normals.push_back(getUnitVector(vectorNormal));
    }
}

vec3<double> Mesh::getIntersectionNormal(vec3<double> intersectionPoint) {
    double x = intersectionPoint.x();
    double y = intersectionPoint.y();
    double z = intersectionPoint.z();
    int firstVertexIdx = mostRecentlyHitTriangleIndex;
    int secondVertexIdx = (firstVertexIdx + 1) % vertices.size();
    int thirdVertexIdx = (firstVertexIdx + 2) % vertices.size();
    
    vec3<double> firstVertex = vertices.at(firstVertexIdx);
    vec3<double> secondVertex = vertices.at(secondVertexIdx);
    vec3<double> thirdVertex = vertices.at(thirdVertexIdx);
    
    vec3<double> firstVector = getUnitVector(firstVertex - thirdVertex);
    vec3<double> secondVector = getUnitVector(firstVertex - secondVertex);
    
    // stack exchange
    vec3<double> v2_1 = firstVertex - secondVertex;
    vec3<double> v2_3 = thirdVertex - secondVertex;
    vec3<double> v2_t = intersectionPoint - secondVertex;
    
    double d00 = dot(v2_1, v2_1);
    double d01 = dot(v2_1, v2_3);
    double d11 = dot(v2_3, v2_3);
    float denom = d00 * d11 - d01 * d01;
    
    float d20 = dot(v2_t, v2_1);
    float d21 = dot(v2_t, v2_3);
    
    double bary_0 = (d11 * d20 - d01 * d21) / denom; // weight related to p1
    double bary_1 = (d00 * d21 - d01 * d20) / denom; // weight related to p3
    double bary_2 = 1.0 - bary_0 - bary_1; // weight related to p2
    
    return bary_0 * normals.at(firstVertexIdx) + bary_2 * normals.at(secondVertexIdx) + bary_1 * normals.at(thirdVertexIdx);
    
}
