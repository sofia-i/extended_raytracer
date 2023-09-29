//
//  TriangleMesh.hpp
//  raytracer_ext
//
//  Created by Sofia Iannicelli on 3/5/23.
//

#ifndef TriangleMesh_hpp
#define TriangleMesh_hpp

#include <stdio.h>
#include "Object.hpp"
#include "Triangle.hpp"

class TriangleMesh : public Object, public Triangle {
public:
    TriangleMesh(std::vector<vec3<double>*> vertices, double kd, double ks, double ka, vec3<double> objectColor, vec3<double> objectSpecular, double kgls, double refl): Object(kd, ks, ka, kgls, objectColor, objectSpecular, refl, ""), Triangle(vertices) {
    }
    
    TriangleMesh(std::vector<vec3<double>*> vertices, double kd, double ks, double ka, vec3<double> objectColor, vec3<double> objectSpecular, double kgls, double refl, std::string description): Object(kd, ks, ka, kgls, objectColor, objectSpecular, refl, description), Triangle(vertices){
    }
    
    vec3<double> getIntersectionNormal(vec3<double> intersectionPoint);
    double findRayObjectIntersection(Ray ray) {
        return Triangle::findRayObjectIntersection(ray);
    }
    
    std::string toString() const {
        std::string str = "";
        std::stringstream ss(str);
        
        ss << getDescription() << std::endl;
        // print out the vertices
        ss << Object::toString() << std::endl;
        
        return ss.str();
    }
};

#endif /* TriangleMesh_hpp */
