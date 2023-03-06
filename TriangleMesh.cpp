//
//  TriangleMesh.cpp
//  raytracer_ext
//
//  Created by Sofia Iannicelli on 3/5/23.
//

#include "TriangleMesh.hpp"
#include "Ray.hpp"

vec3<double> TriangleMesh::getIntersectionNormal(vec3<double> intersectionPoint) {
    return getUnitVector(planeNormal);
}
