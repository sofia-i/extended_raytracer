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
protected:
    std::vector<vec3<double>> vertices;
    vec3<double> planeNormal;
    
    vec3<double> calculatePlaneNormal();
    
public:
    Triangle(std::vector<vec3<double>> vertices) {
        this->vertices = vertices;
        this->planeNormal = calculatePlaneNormal();
    }
    
    vec3<double> getPlaneNormal() {
        return planeNormal;
    }
    
    double findRayObjectIntersection(Ray ray);
};


#endif /* Triangle_h */
