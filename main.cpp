//
//  main.cpp
//  raytracer_ext
//
//  Created by Sofia Iannicelli on 2/15/23.
//

#include <iostream>
#include <fstream>
#include <string>

#include "Scene.hpp"
#include "SceneParser.hpp"
#include "Raytracer.hpp"
#include "PpmWriter.hpp"
#include "Mesh.hpp"

int main_helper() {
    // configure file info
    std::string input_suffix = "3";
    std::string outputFilePath = "/Users/sofiaiannicelli/Documents/BYU_WINTER_2023/graphics/raytracer_3/raytracer_ext/outputs/output_image_1_" + input_suffix + ".ppm";
    std::string inputFilePath = "/Users/sofiaiannicelli/Documents/BYU_WINTER_2023/graphics/raytracer_3/raytracer_ext/inputs/input_1_" + input_suffix + ".txt";

    // parse the scene information from the input file
    SceneParser sceneParser;
    Scene scene = sceneParser.parseFile(inputFilePath);
    std::cout << scene << std::endl;
    
    // send the results to a ppm file for output
    std::string magicNumber = "P3";
    int numColumns = 1000;
    int numRows = 1000;
    int maxColorVal = 255;
    
    Raytracer raytracer = Raytracer(scene, numColumns, numRows);
    
    // do the raytracing process to get the color results for each pixel
    int*** pixelColors = raytracer.raytrace(scene, numColumns, numRows);
    
    std::ofstream outputFile = PpmWriter::writePpm(outputFilePath, magicNumber, numColumns, numRows, maxColorVal, pixelColors);

    // close file
    outputFile.close();
    
    // deallocate memory
    for(int i = 0; i < numRows; ++i) {
        for(int j = 0; j < numColumns; ++j) {
            delete[] pixelColors[i][j];
        }
    }
    for(int i = 0; i < numRows; ++i) {
        delete[] pixelColors[i];
    }
    delete[] pixelColors;
    
    return 0;
}

std::vector<vec3<double>*> getHouseVertices() {
    vec3<double>* vertex1 = new vec3<double>(-0.5, 0.5, 0);
    vec3<double>* vertex2 = new vec3<double>(-0.5, -0.5, 0);
    vec3<double>* vertex3 = new vec3<double>(0.5, -0.5, 0);
    vec3<double>* vertex4 = new vec3<double>(0.5, 0.5, 0);
    vec3<double>* vertex5 = new vec3<double>(0, 0.9, 0);
    
    std::vector<vec3<double>*> vertices;
    vertices.push_back(vertex1);
    vertices.push_back(vertex2);
    vertices.push_back(vertex3);
    vertices.push_back(vertex1);
    vertices.push_back(vertex3);
    vertices.push_back(vertex4);
    vertices.push_back(vertex1);
    vertices.push_back(vertex4);
    vertices.push_back(vertex5);
    
    return vertices;
}

int main_helper2() {
    double pi = 3.1415926535;
    
    std::string input_suffix = "test";
    std::string outputFilePath = "/Users/sofiaiannicelli/Documents/BYU_WINTER_2023/graphics/raytracer_3/raytracer_ext/outputs/output_image_" + input_suffix + ".ppm";
    
    Scene scene = Scene();
    scene.setDirectionToLight(vec3<double>(0, 0.0, 1.0));
    scene.getCamera()->setCameraLookFrom(vec3<double>(0,0,5.0));
    
    std::vector<vec3<double>*> vertices;
    int divisions = 20;
    for(int i = 0; i < divisions; ++i) {
        double radius = 0.5;
        double x_coord = radius * cos(-pi + i * (2*pi/divisions));
        double z_coord = radius * sin(-(-pi + i * (2*pi/divisions)));
        vertices.push_back(new vec3<double>(x_coord, 0.3, z_coord));
        vertices.push_back(new vec3<double>(x_coord, -0.3, z_coord));
    }
    vec3<double>* top = new vec3<double>(0, 0.3, 0);
    vec3<double>* bottom = new vec3<double>(0, -0.3, 0);
    
    
    std::vector<vec3<double>*> mesh_vertices;
    
    
    // main body
     
    for(int i = 0; i < divisions/2; ++i) {
        vec3<double>* vertex0 = vertices.at(2*i);
        vec3<double>* vertex1 = vertices.at(2*i + 1);
        int start_of_next = 2*(i+1) % vertices.size();
        vec3<double>* vertex2 = vertices.at(start_of_next);
        vec3<double>* vertex3 = vertices.at(start_of_next + 1);
        
        mesh_vertices.push_back(vertex0);
        mesh_vertices.push_back(vertex1);
        mesh_vertices.push_back(vertex3);
        
        mesh_vertices.push_back(vertex0);
        mesh_vertices.push_back(vertex3);
        mesh_vertices.push_back(vertex2);
    }
    
    
    /* top
    for(int i = 0; i < divisions; ++i) {
        vec3<double>* vertex0 = vertices.at(2*i);
        vec3<double>* vertex2 = vertices.at((2*(i+1))%(vertices.size()));
        
        mesh_vertices.push_back(vertex0);
        mesh_vertices.push_back(vertex2);
        mesh_vertices.push_back(top);
    }
    */
    
    /* bottom
    for(int i = 0; i < divisions; ++i) {
        vec3<double>* vertex0 = vertices.at(2*i + 1);
        vec3<double>* vertex2 = vertices.at(((2*(i+1))%vertices.size()) + 1);
        
        mesh_vertices.push_back(vertex2);
        mesh_vertices.push_back(vertex0);
        mesh_vertices.push_back(bottom);
    }
    */
    
    // std::vector<vec3<double>*> mesh_vertices = getHouseVertices();
    
    Mesh* mesh = new Mesh(mesh_vertices);
    scene.addObject(mesh);
    std::cout << scene << std::endl;
    
    // send the results to a ppm file for output
    std::string magicNumber = "P3";
    int numColumns = 250;
    int numRows = 250;
    int maxColorVal = 255;
    
    Raytracer raytracer = Raytracer(scene, numColumns, numRows);
    
    // do the raytracing process to get the color results for each pixel
    int*** pixelColors = raytracer.raytrace(scene, numColumns, numRows);
    
    std::ofstream outputFile = PpmWriter::writePpm(outputFilePath, magicNumber, numColumns, numRows, maxColorVal, pixelColors);

    // close file
    outputFile.close();
    
    // deallocate memory
    for(int i = 0; i < numRows; ++i) {
        for(int j = 0; j < numColumns; ++j) {
            delete[] pixelColors[i][j];
        }
    }
    for(int i = 0; i < numRows; ++i) {
        delete[] pixelColors[i];
    }
    delete[] pixelColors;
    
    /*
    for(vec3<double>* vertex : vertices) {
        delete vertex;
    }
    /*
    delete vertex1;
    delete vertex2;
    delete vertex3;
    delete vertex4;
    delete vertex5;
     */
    delete mesh;
    
    return 0;
}

int main_helper3() {
    double pi = 3.1415926535;
    
    std::string input_suffix = "test";
    std::string outputFilePath = "/Users/sofiaiannicelli/Documents/BYU_WINTER_2023/graphics/raytracer_3/raytracer_ext/outputs/output_image_" + input_suffix + ".ppm";
    
    Scene scene = Scene();
    scene.setDirectionToLight(vec3<double>(0, 0.0, 1.0));
    
    std::vector<vec3<double>*> vertices;
    vertices.push_back(new vec3<double>(-0.5, 0.0, 0.0));
    vertices.push_back(new vec3<double>(0.0, 0.0, 0.5));
    vertices.push_back(new vec3<double>(0.0, 0.5, 0.0));
    
    vertices.push_back(new vec3<double>(0.5, 0.0, 0.0));
    
    vertices.push_back(new vec3<double>(0.0, 0.0, -0.5));
    
    std::vector<vec3<double>*> mesh_vertices;
    mesh_vertices.push_back(vertices.at(0));
    mesh_vertices.push_back(vertices.at(1));
    mesh_vertices.push_back(vertices.at(2));
    mesh_vertices.push_back(vertices.at(2));
    mesh_vertices.push_back(vertices.at(1));
    mesh_vertices.push_back(vertices.at(3));
    
    
    mesh_vertices.push_back(vertices.at(3));
    mesh_vertices.push_back(vertices.at(4));
    mesh_vertices.push_back(vertices.at(2));
    mesh_vertices.push_back(vertices.at(2));
    mesh_vertices.push_back(vertices.at(4));
    mesh_vertices.push_back(vertices.at(0));
    
    
    Mesh* mesh = new Mesh(mesh_vertices);
    scene.addObject(mesh);
    std::cout << scene << std::endl;
    
    // send the results to a ppm file for output
    std::string magicNumber = "P3";
    int numColumns = 250;
    int numRows = 250;
    int maxColorVal = 255;
    
    Raytracer raytracer = Raytracer(scene, numColumns, numRows);
    
    // do the raytracing process to get the color results for each pixel
    int*** pixelColors = raytracer.raytrace(scene, numColumns, numRows);
    
    std::ofstream outputFile = PpmWriter::writePpm(outputFilePath, magicNumber, numColumns, numRows, maxColorVal, pixelColors);

    // close file
    outputFile.close();
    
    // deallocate memory
    for(int i = 0; i < numRows; ++i) {
        for(int j = 0; j < numColumns; ++j) {
            delete[] pixelColors[i][j];
        }
    }
    for(int i = 0; i < numRows; ++i) {
        delete[] pixelColors[i];
    }
    delete[] pixelColors;
    
    
    for(vec3<double>* vertex : vertices) {
        delete vertex;
    }
    /*
    delete vertex1;
    delete vertex2;
    delete vertex3;
    delete vertex4;
    delete vertex5;
     */
    delete mesh;
    
    return 0;
}

int main() {
    return main_helper2();
}
