#ifndef INPUT_WINDOW_HPP
#define INPUT_WINDOW_HPP

#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include <glad/glad.h>
#include <GLFW/glfw3.h>
#include <application_settings.hpp>
#include "camera.h"

// Camera
Camera cam(0.0f, 0.0f, 8.0f);
void framebuffer_size_callback(GLFWwindow *window, int width, int height)
{
    glViewport(0, 0, width, height);
}

class InputWindow
{
public:
    //////////////////////////////////////////////////////////////////////////////////////////////
    // constructor
    InputWindow(unsigned int width, unsigned int height, const char *application_name)
    {
        iWindow = glfwCreateWindow(width, height, application_name, NULL, NULL);
        if (iWindow == NULL)
        {
            std::cout << "Failed to create GLFW window" << std::endl;
            glfwTerminate();
            window_error = true;
        }
        glfwMakeContextCurrent(iWindow);
        glfwSetFramebufferSizeCallback(iWindow, framebuffer_size_callback);

        if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress))
        {
            std::cout << "Failed to initialize GLAD" << std::endl;
            window_error = true;
        }

        glfwSetWindowUserPointer(iWindow, this);
        glfwSetKeyCallback(iWindow, onKey);
    }

    ~InputWindow()
    {
        glfwDestroyWindow(iWindow);
    }
    //////////////////////////////////////////////////////////////////////////////////////////////
    // functions
    void onKey(int key, int scancode, int action, int mods)
    {
        bool settings_updated = false;
        if (key == GLFW_KEY_P && action == GLFW_PRESS)
        {
            bool current_value = linkedSettings->getDebugMode();
            linkedSettings->changeDebugMode(!current_value);
            if (current_value)
            {
                std::cout << "debug mode set to OFF" << std::endl;
            }
            else
            {
                std::cout << "debug mode set to ON" << std::endl;
            }
            settings_updated = true;
        }

        if (key == GLFW_KEY_ESCAPE && action == GLFW_PRESS)
            glfwSetWindowShouldClose(iWindow, true);

        if (key == GLFW_KEY_M && action == GLFW_PRESS)
        {
            double xpos, ypos;
            glfwGetCursorPos(iWindow, &xpos, &ypos);
            linkedSettings->setMousePos(xpos, ypos);
            std::cout << "current mouse position is: " << xpos << ", " << ypos << std::endl;
            settings_updated = true;
        }

        if (key == GLFW_KEY_1 && action == GLFW_PRESS)
        {
            linkedSettings->setEpsilon(0.1f);
            std::cout << "epsilon changed to 0.1" << std::endl;
            settings_updated = true;
        }
        if (key == GLFW_KEY_2 && action == GLFW_PRESS)
        {
            linkedSettings->setEpsilon(0.01f);
            std::cout << "epsilon changed to 0.01" << std::endl;
            settings_updated = true;
        }
        if (key == GLFW_KEY_3 && action == GLFW_PRESS)
        {
            linkedSettings->setEpsilon(0.001f);
            std::cout << "epsilon changed to 0.001" << std::endl;
            settings_updated = true;
        }
        if (key == GLFW_KEY_4 && action == GLFW_PRESS)
        {
            linkedSettings->setEpsilon(0.0001f);
            std::cout << "epsilon changed to 0.0001" << std::endl;
            settings_updated = true;
        }
        if (key == GLFW_KEY_5 && action == GLFW_PRESS)
        {
            linkedSettings->setEpsilon(0.00001f);
            std::cout << "epsilon changed to 0.00001" << std::endl;
            settings_updated = true;
        }
        if (key == GLFW_KEY_X && action == GLFW_PRESS)
        {
            cam.resetCamera(glm::vec3(8.0f, 0.0f, 0.0f));
        }
        if (key == GLFW_KEY_Y && action == GLFW_PRESS)
        {
            cam.resetCamera(glm::vec3(0.0f, 8.0f, 0.0f));
        }
        if (key == GLFW_KEY_Z && action == GLFW_PRESS)
        {
            cam.resetCamera(glm::vec3(0.0f, 0.0f, 8.0f));
        }

        if (settings_updated)
            linkedSettings->update();
    }

    bool checkWindowCreation()
    {
        return window_error;
    }

    GLFWwindow *getWindowPointer()
    {
        return iWindow;
    }

    void linkSettings(AppSettings *settings)
    {
        linkedSettings = settings;
    }

private:
    GLFWwindow *iWindow;
    bool window_error = false;
    AppSettings *linkedSettings;

    static void onKey(GLFWwindow *window, int key, int scancode, int actions, int mods)
    {
        InputWindow *obj = (InputWindow *)glfwGetWindowUserPointer(window);
        obj->onKey(key, scancode, actions, mods);
    }
};

//////////////////////////////////////////////////////////////////////////////////////////////
// additional functions outside class

// camera is global in this header file to be easily accessible in callback functions
// this can be used to access the camera from other parts of the code
Camera *getCameraPointer()
{
    return &cam;
}

float lastX;
float lastY;
bool firstMouse = true;

glm::vec2 mouse_offset(GLFWwindow *window)
{
    double xpos, ypos;
    glfwGetCursorPos(window, &xpos, &ypos);
    if (firstMouse)
    {
        lastX = xpos;
        lastY = ypos;
        firstMouse = false;
    }

    float xoffset = xpos - lastX;
    float yoffset = lastY - ypos; // reversed since y-coordinates go from bottom to top

    lastX = xpos;
    lastY = ypos;
    // std::cout << xoffset << std::endl;

    return glm::vec2(xoffset, yoffset);
}

void mouse_callback(GLFWwindow *window, double xpos, double ypos)
{
    if (glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_RIGHT) == GLFW_PRESS)
    {
        // glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_DISABLED);
        if (firstMouse)
        {
            lastX = xpos;
            lastY = ypos;
            firstMouse = false;
        }

        float xOffset = xpos - lastX;
        float yOffset = ypos - lastY;
        lastX = xpos;
        lastY = ypos;

        cam.ProcessCameraRotation(xOffset, yOffset);
    }
    if (glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_RIGHT) == GLFW_RELEASE)
    {
        // glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_NORMAL);
    }
}

void scroll_callback(GLFWwindow *window, double xoffset, double yoffset)
{
    cam.ProcessMouseScroll(yoffset);
}

void processInput(GLFWwindow *window, AppSettings *settings, float deltaTime)
{
    bool settings_updated = false;

    // close on ESC
    if (glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS)
        glfwSetWindowShouldClose(window, GL_TRUE);

    // Up & Down Arrows to change SOlvent Radius
    if (glfwGetKey(window, GLFW_KEY_UP) == GLFW_PRESS)
    {
        settings->addToSolventRadius(0.01f);
        settings_updated = true;
    }
    if (glfwGetKey(window, GLFW_KEY_DOWN) == GLFW_PRESS)
    {
        settings->addToSolventRadius(-0.01f);
        settings_updated = true;
    }

    // camera movement with WASD + IK
    if (glfwGetKey(window, GLFW_KEY_W) == GLFW_PRESS)
        cam.ProcessCameraMovement(UPWARD, deltaTime);
    if (glfwGetKey(window, GLFW_KEY_S) == GLFW_PRESS)
        cam.ProcessCameraMovement(DOWNWARD, deltaTime);
    if (glfwGetKey(window, GLFW_KEY_A) == GLFW_PRESS)
        cam.ProcessCameraMovement(LEFT, deltaTime);
    if (glfwGetKey(window, GLFW_KEY_D) == GLFW_PRESS)
        cam.ProcessCameraMovement(RIGHT, deltaTime);
    if (glfwGetKey(window, GLFW_KEY_I) == GLFW_PRESS)
        cam.ProcessCameraMovement(FORWARD, deltaTime);
    if (glfwGetKey(window, GLFW_KEY_K) == GLFW_PRESS)
        cam.ProcessCameraMovement(BACKWARD, deltaTime);

    if (settings_updated)
        settings->update();
}

#endif