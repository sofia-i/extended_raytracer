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
#include <map>
#include "Object.hpp"
#include "Triangle.hpp"

class Mesh : public Object {
private:
    std::map<vec3<double>*, std::vector<Triangle*>> vertices;
    std::map<vec3<double>*, vec3<double>> vertex_normals;
    // std::vector<vec3<double>> vertex_normals;
    std::vector<Triangle*> triangles;
    
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
    bool intersectsBB(Ray ray);
    double findRayRectangleIntersection(std::vector<vec3<double>> vertices, Ray ray);
    double findRayPlaneIntersection(std::vector<vec3<double>> vertices, Ray ray);
    
public:
    Mesh(std::vector<vec3<double>*> vertices, double kd, double ks, double ka, double kgls, vec3<double> objectColor, vec3<double> objectSpecular, double refl, std::string description) : Object(kd, ks, ka, kgls, objectColor, objectSpecular, refl, description) {
        
        if(vertices.size() % 3 != 0) {
            std::cerr << "vertices size in Mesh " << description << " initialization not a multiple of 3." << std::endl;
        }
        
        for(int i = 0; i < vertices.size()-2; i += 3) {
            std::vector<vec3<double>*> triangle_vertices;
            triangle_vertices.push_back(vertices.at(i));
            triangle_vertices.push_back(vertices.at(i+1));
            triangle_vertices.push_back(vertices.at(i+2));
            
            // TODO: clear memory
            Triangle* triangle = new Triangle(triangle_vertices);
            this->triangles.push_back(triangle);
            
            for(int j = i; j < i+3; ++j) {
                vec3<double>* vertex = vertices.at(j);
                if(this->vertices.find(vertex) == this->vertices.end()) {
                    // not in map
                    std::vector<Triangle*> connected_triangles;
                    connected_triangles.push_back(triangle);
                    
                    this->vertices.insert({vertex, connected_triangles});
                }
                else {
                    // add the triangle to connected triangles for the vertex
                    (this->vertices.at(vertex)).push_back(triangle);
                }
            }
        }
        // this->vertices = vertices;
        computeVertexNormals();
        computeBBInfo();
    }
    
    Mesh(std::vector<vec3<double>*> vertices) : Object(0.7, 0.2, 0.1, 4.0, vec3<double>(1.0, 1.0, 1.0), vec3<double>(1.0, 1.0, 1.0), 0.2, "") {
        
        if(vertices.size() % 3 != 0) {
            std::cerr << "vertices size in Mesh " << description << " initialization not a multiple of 3." << std::endl;
        }
        
        for(int i = 0; i < vertices.size()-2; i += 3) {
            std::vector<vec3<double>*> triangle_vertices;
            triangle_vertices.push_back(vertices.at(i));
            triangle_vertices.push_back(vertices.at(i+1));
            triangle_vertices.push_back(vertices.at(i+2));
            
            // TODO: clear memory
            Triangle* triangle = new Triangle(triangle_vertices);
            this->triangles.push_back(triangle);
            
            for(int j = i; j < i+3; ++j) {
                vec3<double>* vertex = vertices.at(j);
                if(this->vertices.find(vertex) == this->vertices.end()) {
                    // not in map
                    std::vector<Triangle*> connected_triangles;
                    connected_triangles.push_back(triangle);
                    
                    this->vertices.insert({vertex, connected_triangles});
                }
                else {
                    // add the triangle to connected triangles for the vertex
                    (this->vertices.at(vertex)).push_back(triangle);
                }
            }
        }
        // this->vertices = vertices;
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
