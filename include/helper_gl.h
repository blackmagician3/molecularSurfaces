#ifndef HELPER_GL_H
#define HELPER_GL_H
#include <glad/glad.h> // include glad to get all the required OpenGL headers
#include <glm/glm.hpp>
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include "shader.h"
#include "camera.hpp"

// TODO add to config
const unsigned int SCR_WIDTH = 1920;
const unsigned int SCR_HEIGHT = 1080;

float vertices[] =
    {
        // Coords     // texCoords
        1.0f, -1.0f, 1.0f, 0.0f,
        -1.0f, -1.0f, 0.0f, 0.0f,
        -1.0f, 1.0f, 0.0f, 1.0f,

        1.0f, 1.0f, 1.0f, 1.0f,
        1.0f, -1.0f, 1.0f, 0.0f,
        -1.0f, 1.0f, 0.0f, 1.0f};

unsigned int indices[] = {
    0, 1, 3, // first Triangle
    1, 2, 3  // second Triangle
};           // note that we start from 0!

// TODO: cleanup: this header is named poorly
int quadBuffers()
{
    // generate VAO, VBO
    unsigned int VAO, VBO;
    glGenVertexArrays(1, &VAO);
    glGenBuffers(1, &VBO);

    // binding
    glBindVertexArray(VAO);
    glBindBuffer(GL_ARRAY_BUFFER, VBO);
    glBufferData(GL_ARRAY_BUFFER, sizeof(vertices), &vertices, GL_STATIC_DRAW);

    // attribute array
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 4 * sizeof(float), (void *)0);
    glEnableVertexAttribArray(1);
    glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, 4 * sizeof(float), (void *)(2 * sizeof(float)));

    return VAO;
}

void draw(unsigned int width, unsigned int height, int VAO, GLuint *texture_ptr)
{
    // Render to the screen
    glClearColor(1.0f, 0.0f, 0.0f, 0.0f);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    glBindVertexArray(VAO);
    glDisable(GL_DEPTH_TEST);
    glBindTexture(GL_TEXTURE_2D, *texture_ptr);
    glDrawArrays(GL_TRIANGLES, 0, 6);
}

// test the buffer (write the first 20000 elements to a txt file)
void readTextureToFile(GLuint texture, int size)
{
    float4 *data = new float4[size];
    glGetTextureImage(texture, 0, GL_RGBA, GL_FLOAT, size * sizeof(float4), data);

    // access the data
    std::ofstream myfile("test.txt");
    if (myfile.is_open())
    {
        for (int i = 0; i < size; i++)
        {
            myfile << data[i].x << ", " << data[i].y << ", " << data[i].z << ", " << data[i].w << "\n";
        }
        myfile.close();
    }
    else
        std::cout << "Unable to open file" << std::endl;
}

/////////////////////////////////////////////////////////////////////////////
//////// Initialize GL //////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////

#endif