//
//  Mesh.hpp
//  raytracer_ext
//
//  Created by Sofia Iannicelli on 3/4/23.
//

#ifndef Mesh_hpp
#define Mesh_hpp

#include <stdio.h>
#include <vector>
#include "Object.hpp"

class Mesh : public Object {
private:
    std::vector<vec3<double>> vertices;
    std::vector<vec3<double>> normals;
    int mostRecentlyHitTriangleIndex;
    
    void computeVertexNormals();
    void computeBBInfo();
    
private:
    double bb_xmin;
    double bb_xmax;
    double bb_ymin;
    double bb_ymax;
    double bb_zmin;
    double bb_zmax;
    
    bool fallsInBoundingBox(vec3<double> point);
    double findRayBBIntersection(Ray ray);
    double findRayRectangleIntersection(std::vector<vec3<double>> vertices, Ray ray);
    
public:
    Mesh(std::vector<vec3<double>> vertices, double kd, double ks, double ka, double kgls, vec3<double> objectColor, vec3<double> objectSpecular, double refl, std::string description) : Object(kd, ks, ka, kgls, objectColor, objectSpecular, refl, description) {
        this->vertices = vertices;
        computeVertexNormals();
        computeBBInfo();
    }
    
    Mesh(std::vector<vec3<double>> vertices) : Object(0.7, 0.2, 0.1, 4.0, vec3<double>(1.0, 1.0, 1.0), vec3<double>(1.0, 1.0, 1.0), 0.2, "") {
        this->vertices = vertices;
        computeVertexNormals();
        computeBBInfo();
    }
    
    double findRayObjectIntersection(Ray ray);
    vec3<double> getIntersectionNormal(vec3<double> intersectionPoint);
    
    std::string toString() const {
        return "";
    }
    
    
};

#endif /* Mesh_hpp */
