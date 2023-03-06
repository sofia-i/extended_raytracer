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

int main() {
    // main_helper();
    std::string input_suffix = "test";
    std::string outputFilePath = "/Users/sofiaiannicelli/Documents/BYU_WINTER_2023/graphics/raytracer_3/raytracer_ext/outputs/output_image_" + input_suffix + ".ppm";
    
    Scene scene = Scene();
    std::vector<vec3<double>> vertices;
    /*
    float vertices[] = {
                // positions            // colors
                -0.5f, 0.5f, 0.0f,      1.0f, 0.0f, 0.0f,
                -0.5f, -0.5f, 0.0f,     0.0f, 1.0f, 0.0f,
                0.5f, -0.5f, 0.0f,      0.0f, 0.0f, 1.0f,
                0.5f, 0.5f, 0.0f,       1.0f, 1.0f, 0.0f,
                0.0f, 0.9f, 0.0f,       1.0f, 0.0f, 0.0f
        };
     */
    vec3<double> vertex1 = vec3<double>(-0.5, 0.5, 0);
    vertices.push_back(vertex1);
    vec3<double> vertex2 = vec3<double>(-0.5, -0.5, 0);
    vertices.push_back(vertex2);
    vec3<double> vertex3 = vec3<double>(0.5, -0.5, 0);
    vertices.push_back(vertex3);
    vec3<double> vertex4 = vec3<double>(0.5, 0.5, 0);
    vertices.push_back(vertex4);
    vec3<double> vertex5 = vec3<double>(0, 0.9, 0);
    vertices.push_back(vertex5);
    Mesh* mesh = new Mesh(vertices);
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
    
    return 0;
}
