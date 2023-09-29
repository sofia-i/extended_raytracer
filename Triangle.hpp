//
//  Triangle.h
//  raytracer_2
//
//  Created by Sofia Iannicelli on 3/4/23.
//

#ifndef Triangle_hpp
#define Triangle_hpp

#include <stdio.h>
#include <vector>
#include "vec3.hpp"
#include "Ray.hpp"

class Triangle {
private:
    double EPSILON = 2.0e-8;
    void initializeDistPlaneToOrigin() {
        double a = planeNormal[0];
        double b = planeNormal[1];
        double c = planeNormal[2];
        
        vec3<double>* vertex = vertices.at(0);
        
        distPlaneToOrigin = -a * vertex->x() - b * vertex->y() - c * vertex->z();
    }
    
protected:
    std::vector<vec3<double>*> vertices;
    vec3<double> planeNormal;
    double distPlaneToOrigin;
    
    vec3<double> calculatePlaneNormal();
    
public:
    Triangle(std::vector<vec3<double>*> vertices) {
        this->vertices = vertices;
        this->planeNormal = calculatePlaneNormal();
        initializeDistPlaneToOrigin();
    }
    
    vec3<double> getPlaneNormal() {
        return planeNormal;
    }
    std::vector<vec3<double>*> getVertices() {
        return vertices;
    }
    
    double findRayObjectIntersection(Ray ray);
    bool isInTriangle(vec3<double> intersectionPoint);
    bool isInPlane(vec3<double> point);
};


#endif /* Triangle_h */
