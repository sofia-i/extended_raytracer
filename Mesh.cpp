//
//  Mesh.cpp
//  raytracer_ext
//
//  Created by Sofia Iannicelli on 3/4/23.
//

#include "Mesh.hpp"
#include "Triangle.hpp"
#include <vector>
#include <map>

void Mesh::computeBBInfo() {
    std::map<vec3<double>*, std::vector<Triangle*>>::iterator it;

    it = vertices.begin();
    vec3<double>* first_vertex = it->first;
    
    bb_xmin = first_vertex->x();
    bb_xmax = first_vertex->x();
    bb_ymin = first_vertex->y();
    bb_ymax = first_vertex->y();
    bb_zmin = first_vertex->z();
    bb_zmax = first_vertex->z();
    
    for (it = vertices.begin(); it != vertices.end(); it++)
    {
        vec3<double>* vertex = it->first;
        
        if(vertex->x() < bb_xmin) {
            bb_xmin = vertex->x();
        }
        else if(vertex->x() > bb_xmax) {
            bb_xmax = vertex->x();
        }
        
        if(vertex->y() < bb_ymin) {
            bb_ymin = vertex->y();
        }
        else if(vertex->y() > bb_ymax) {
            bb_ymax = vertex->y();
        }
        
        if(vertex->z() < bb_zmin) {
            bb_zmin = vertex->z();
        }
        else if(vertex->z() > bb_zmax) {
            bb_zmax = vertex->z();
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

/* REALLY BIG ASSUMPTION HERE IS THAT THE */
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

/* REALLY BIG ASSUMPTION HERE IS THAT THE */
double Mesh::findRayPlaneIntersection(std::vector<vec3<double>> vertices, Ray ray) {
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
    
    return t;
}

bool Mesh::intersectsBB(Ray ray) {
    // bounding box is 6 rectangles
    std::vector<vec3<double>> frontVertices{
        vec3<double>(bb_xmin, bb_ymin, bb_zmax),
        vec3<double>(bb_xmax, bb_ymin, bb_zmax),
        vec3<double>(bb_xmax, bb_ymax, bb_zmax),
        vec3<double>(bb_xmin, bb_ymin, bb_zmax)
    };
    double frontResult = findRayPlaneIntersection(frontVertices, ray);
    if(frontResult > 0) {
        vec3<double> intersectionPoint = ray.getPointOnRay(frontResult);
        // pass x and y bounds
        double intersection_x = intersectionPoint.x();
        double intersection_y = intersectionPoint.y();
        if(intersection_x >= bb_xmin && intersection_x <= bb_xmax &&
           intersection_y >= bb_ymin && intersection_y <= bb_ymax) {
            return true;
        }
    }
    
    std::vector<vec3<double>> backVertices {
        vec3<double>(bb_xmin, bb_ymin, bb_zmin),
        vec3<double>(bb_xmax, bb_ymin, bb_zmin),
        vec3<double>(bb_xmax, bb_ymax, bb_zmin),
        vec3<double>(bb_xmin, bb_ymin, bb_zmin)
    };
    double backResult = findRayPlaneIntersection(backVertices, ray);
    if(backResult > 0) {
        vec3<double> intersectionPoint = ray.getPointOnRay(backResult);
        // pass x and y bounds
        double intersection_x = intersectionPoint.x();
        double intersection_y = intersectionPoint.y();
        if(intersection_x >= bb_xmin && intersection_x <= bb_xmax &&
           intersection_y >= bb_ymin && intersection_y <= bb_ymax) {
            return true;
        }
    }
    
    std::vector<vec3<double>> leftVertices{
        vec3<double>(bb_xmin, bb_ymin, bb_zmin),
        vec3<double>(bb_xmin, bb_ymin, bb_zmax),
        vec3<double>(bb_xmin, bb_ymax, bb_zmax),
        vec3<double>(bb_xmin, bb_ymax, bb_zmin)
    };
    double leftResult = findRayPlaneIntersection(leftVertices, ray);
    if(leftResult > 0) {
        vec3<double> intersectionPoint = ray.getPointOnRay(leftResult);
        // pass y and z bounds
        double intersection_y = intersectionPoint.y();
        double intersection_z = intersectionPoint.z();
        if(intersection_y >= bb_ymin && intersection_y <= bb_ymax &&
           intersection_z >= bb_zmin && intersection_z <= bb_zmax) {
            return true;
        }
    }
    
    std::vector<vec3<double>> rightVertices{
        vec3<double>(bb_xmax, bb_ymin, bb_zmin),
        vec3<double>(bb_xmax, bb_ymin, bb_zmax),
        vec3<double>(bb_xmax, bb_ymax, bb_zmax),
        vec3<double>(bb_xmax, bb_ymax, bb_zmin)
    };
    double rightResult = findRayPlaneIntersection(rightVertices, ray);
    if(rightResult > 0) {
        vec3<double> intersectionPoint = ray.getPointOnRay(rightResult);
        // pass y and z bounds
        double intersection_y = intersectionPoint.y();
        double intersection_z = intersectionPoint.z();
        if(intersection_y >= bb_ymin && intersection_y <= bb_ymax &&
           intersection_z >= bb_zmin && intersection_z <= bb_zmax) {
            return true;
        }
    }
    
    std::vector<vec3<double>> topVertices{
        vec3<double>(bb_xmin, bb_ymax, bb_zmax),
        vec3<double>(bb_xmax, bb_ymax, bb_zmax),
        vec3<double>(bb_xmax, bb_ymax, bb_zmin),
        vec3<double>(bb_xmin, bb_ymax, bb_zmin)
    };
    double topResult = findRayPlaneIntersection(topVertices, ray);
    if(topResult > 0) {
        vec3<double> intersectionPoint = ray.getPointOnRay(topResult);
        // pass x and z bounds
        double intersection_x = intersectionPoint.x();
        double intersection_y = intersectionPoint.y();
        if(intersection_x >= bb_xmin && intersection_x <= bb_xmax &&
           intersection_y >= bb_ymin && intersection_y <= bb_ymax) {
            return true;
        }
    }
    
    std::vector<vec3<double>> bottomVertices{
        vec3<double>(bb_xmin, bb_ymin, bb_zmax),
        vec3<double>(bb_xmax, bb_ymin, bb_zmax),
        vec3<double>(bb_xmax, bb_ymin, bb_zmin),
        vec3<double>(bb_xmin, bb_ymin, bb_zmin)
    };
    double bottomResult = findRayPlaneIntersection(bottomVertices, ray);
    if(bottomResult > 0) {
        vec3<double> intersectionPoint = ray.getPointOnRay(bottomResult);
        // pass x and z bounds
        double intersection_x = intersectionPoint.x();
        double intersection_y = intersectionPoint.y();
        if(intersection_x >= bb_xmin && intersection_x <= bb_xmax &&
           intersection_y >= bb_ymin && intersection_y <= bb_ymax) {
            return true;
        }
    }
    
    return false;
}

double Mesh::findRayObjectIntersection(Ray ray) {
    // check if intersects with bounding box
    if(!intersectsBB(ray)) {
        return -1.0;
    }
    
    // else check all triangles
    for(int i = 0; i < triangles.size(); ++i) {
        Triangle* triangle = triangles.at(i);
        double testT = triangle->findRayObjectIntersection(ray);
        if(testT > 0) {
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
    
    // each vertex
    std::map<vec3<double>*, std::vector<Triangle*>>::iterator it;
    for (it = vertices.begin(); it != vertices.end(); it++) {
        vec3<double> planeNormalSum = vec3<double>(0,0,0);
        
        vec3<double>* vertex = it->first;
        std::vector<Triangle*> connected_triangles = it->second;
        
        // each triangle connected to vertex
        for(Triangle* triangle : connected_triangles) {
            planeNormalSum += triangle->getPlaneNormal();
        }
        
        vec3<double> vectorNormal = getUnitVector(planeNormalSum);
        vertex_normals.insert({vertex, vectorNormal});
    }
}

vec3<double> Mesh::getIntersectionNormal(vec3<double> intersectionPoint) {
    Triangle* hitTriangle = nullptr;
    // check all triangles
    for(int i = 0; i < triangles.size(); ++i) {
        Triangle* triangle = triangles.at(i);
        bool pointInPlane = triangle->isInPlane(intersectionPoint);
        bool pointInTriangle = triangle->isInTriangle(intersectionPoint);
        if(pointInPlane && pointInTriangle) {
            hitTriangle = triangle;
            break;
        }
    }
    
    std::vector<vec3<double>*> triangle_vertices = hitTriangle->getVertices();
    vec3<double>* firstVertex = triangle_vertices.at(0);
    vec3<double>* secondVertex = triangle_vertices.at(1);
    vec3<double>* thirdVertex = triangle_vertices.at(2);
    
    // stack exchange
    vec3<double> v2_1 = *firstVertex - *secondVertex;
    vec3<double> v2_3 = *thirdVertex - *secondVertex;
    vec3<double> v2_t = intersectionPoint - *secondVertex;
    
    double d00 = dot(v2_1, v2_1);
    double d01 = dot(v2_1, v2_3);
    double d11 = dot(v2_3, v2_3);
    float denom = d00 * d11 - d01 * d01;
    
    float d20 = dot(v2_t, v2_1);
    float d21 = dot(v2_t, v2_3);
    
    double bary_0 = (d11 * d20 - d01 * d21) / denom; // weight related to p1
    double bary_1 = (d00 * d21 - d01 * d20) / denom; // weight related to p3
    double bary_2 = 1.0 - bary_0 - bary_1; // weight related to p2
    
    vec3<double> firstVertexNormal = this->vertex_normals.at(firstVertex);
    vec3<double> secondVertexNormal = this->vertex_normals.at(secondVertex);
    vec3<double> thirdVertexNormal = this->vertex_normals.at(thirdVertex);
    
    return bary_0 * firstVertexNormal + bary_2 * secondVertexNormal + bary_1 * thirdVertexNormal;
    
}
