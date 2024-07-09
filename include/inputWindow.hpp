#ifndef INPUT_WINDOW_HPP
#define INPUT_WINDOW_HPP

#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <chrono>

#include <glad/glad.h>
#include <GLFW/glfw3.h>

#include "imgui.h"
#include "imgui_impl_glfw.h"
#include "imgui_impl_opengl3.h"
#include "portable_file_dialog.h"

#include <application_settings.hpp>
#include "camera.hpp"
#include <cuda_inits.cuh>

#ifndef MAXTHREADS
#define MAXTHREADS 512
#endif

struct resolution
{
    uint width;
    uint height;
};
// Camera
Camera cam(0.0f, 0.0f, 20.0f);
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
        if (key == GLFW_KEY_ESCAPE && action == GLFW_PRESS)
            glfwSetWindowShouldClose(iWindow, true);

        if (key == GLFW_KEY_M && action == GLFW_PRESS)
        {
            double xpos, ypos;
            int height, width;
            glfwGetCursorPos(iWindow, &xpos, &ypos);
            glfwGetWindowSize(iWindow, &width, &height);
            linkedSettings->setMousePos(xpos, height - ypos);
            std::cout << "current mouse position is: " << xpos << ", " << height - ypos << std::endl;
            settings_updated = true;
        }
        if (key == GLFW_KEY_T && action == GLFW_PRESS)
        {
            linkedSettings->setKnearest(15);
            std::cout << "k changed to 15" << std::endl;
            settings_updated = true;
        }

        if (key == GLFW_KEY_I && action == GLFW_PRESS)
        {
            linkedSettings->show_gui = !linkedSettings->show_gui;
        }
        if (key == GLFW_KEY_P && action == GLFW_PRESS)
        {
            linkedSettings->show_performance_tool = !linkedSettings->show_performance_tool;
        }
        if (key == GLFW_KEY_R && action == GLFW_PRESS)
        {
            linkedSettings->reset_view();
            cam.intializeCameraPosition(linkedSettings->getCameraFocus());
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

void mouse_callback(GLFWwindow *window, double xpos, double ypos)
{
    if (glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_RIGHT) == GLFW_PRESS)
    {
        // glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_DISABLED);
        if (firstMouse)
        {
            lastX = xpos;
            lastY = ypos;
            cam.updateMouseAxis();
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
        glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_NORMAL);
        firstMouse = true;
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

    // camera movement
    if (glfwGetKey(window, GLFW_KEY_UP) == GLFW_PRESS)
        cam.ProcessCameraMovement(UPWARD, deltaTime);
    if (glfwGetKey(window, GLFW_KEY_DOWN) == GLFW_PRESS)
        cam.ProcessCameraMovement(DOWNWARD, deltaTime);
    if (glfwGetKey(window, GLFW_KEY_LEFT) == GLFW_PRESS)
        cam.ProcessCameraMovement(LEFT, deltaTime);
    if (glfwGetKey(window, GLFW_KEY_RIGHT) == GLFW_PRESS)
        cam.ProcessCameraMovement(RIGHT, deltaTime);
    if (glfwGetKey(window, GLFW_KEY_PAGE_UP) == GLFW_PRESS)
        cam.ProcessCameraMovement(FORWARD, deltaTime);
    if (glfwGetKey(window, GLFW_KEY_PAGE_DOWN) == GLFW_PRESS)
        cam.ProcessCameraMovement(BACKWARD, deltaTime);

    if (settings_updated)
        settings->update();
}

// imgui
void setupImgui(AppSettings *settings)
{

    IM_ASSERT(ImGui::GetCurrentContext() != NULL && "ERROR: Missing Dear ImGui context.\n");

    // some tools
    // static bool show_tool_metrics = false;
    // static bool show_tool_debug_log = false;

    // some flags
    // static bool flag_1 = false;
    // static bool flag_2 = false;
    // static bool flag_3 = false;
    // static bool flag_4 = false;
    // static bool flag_5 = false;
    // static bool flag_6 = false;

    ImGuiWindowFlags window_flags = 0;
    window_flags |= ImGuiWindowFlags_NoTitleBar;
    // window_flags |= ImGuiWindowFlags_NoScrollbar;
    window_flags |= ImGuiWindowFlags_NoMove;
    window_flags |= ImGuiWindowFlags_NoResize;
    window_flags |= ImGuiWindowFlags_NoCollapse;
    window_flags |= ImGuiWindowFlags_NoNav;
    // window_flags |= ImGuiWindowFlags_NoBackground;
    // window_flags |= ImGuiWindowFlags_NoBringToFrontOnFocus;

    ImGuiTreeNodeFlags header_flags = ImGuiTreeNodeFlags_DefaultOpen;

    // size
    const ImGuiViewport *main_viewport = ImGui::GetMainViewport();

    ImGui::SetNextWindowPos(ImVec2(0, 0));
    ImGui::SetNextWindowSize(ImVec2(420, main_viewport->Size.y));

    if (!ImGui::Begin("Molecule Interface", &settings->show_gui, window_flags))
    {
        ImGui::End();
        return;
    }

    // ImGui::PushItemWidth(ImGui::GetFontSize() * -12);
    ImGui::PushItemWidth(ImGui::GetFontSize() * -16);

    // SimulationParams params = settings->getAllHostParams();
    SimulationParams *p_params = settings->getParamsPointer();
    // start of gui
    ImGui::Text("Imgui version (%s) (%d)", IMGUI_VERSION, IMGUI_VERSION_NUM);
    ImGui::Spacing();

    ImGuiTableFlags table_flags = ImGuiTableFlags_Borders;
    if (ImGui::CollapsingHeader("Molecule", header_flags))
    {
        ImGui::SeparatorText("INFORMATION:");
        if (ImGui::BeginTable("Molecule Information", 2, table_flags))
        {
            ImGui::TableNextColumn();
            ImGui::Text("Molecule designation");
            ImGui::TableNextColumn();
            ImGui::Text("tbd");
            ImGui::TableNextColumn();
            ImGui::Text("Number of atoms");
            ImGui::TableNextColumn();
            ImGui::Text("%i", settings->getAtomCount());
            ImGui::TableNextColumn();
            ImGui::Text("Molecule dimensions [A]");
            ImGui::TableNextColumn();
            ImGui::Text("(%.2f, %.2f, %.2f)", p_params->box_start.x, p_params->box_start.y, p_params->box_start.z);

            ImGui::EndTable();
        }

        // TODO: import molecule (button)
        if (ImGui::Button("Import"))
        {
            if (!pfd::settings::available())
            {
                std::cout << "ERROR: Portable File Dialogs are not available on this platform.\n";
            }
            else
            {
                auto selection = pfd::open_file("Select a file to load", ".",
                                                {"pdb", "*.pdb *.pdb1",
                                                 "All Files", "*"})
                                     .result();

                settings->loadMolecule(selection[0]);
            }
        }
        if (ImGui::BeginItemTooltip())
        {
            ImGui::Text("Import molecule from filepath. Currently the only supported fileformat is .pdb.");
            ImGui::EndTooltip();
        }
    }
    ImGui::Dummy(ImVec2(0.0f, 20.0f));
    if (ImGui::CollapsingHeader("Configuration", header_flags))
    {
        ImGui::SeparatorText("PARAMETERS:");
        // TODO: add fields

        // 1. use voxel (checkbox)
        if (ImGui::Checkbox("use grid", &p_params->use_voxel))
            settings->update();
        if (ImGui::BeginItemTooltip())
        {
            ImGui::Text("When using the grid option, individual atoms of the molecule are\npositioned in fields of a 3-dimensional grid based on their position.\nThis especially improves the performance of very large molecules.");
            ImGui::EndTooltip();
        }

        // 2. solver (combo box)
        bool no_arrow = false;
        bool no_preview = false;
        bool width_preview = false;
        bool h_small = false;
        bool h_regular = false;
        bool h_large = false;

        static ImGuiComboFlags combo_flags = 0;
        combo_flags |= ImGuiComboFlags_PopupAlignLeft;
        if (no_arrow)
            combo_flags |= ImGuiComboFlags_NoArrowButton;
        if (no_preview)
            combo_flags |= ImGuiComboFlags_NoPreview;
        if (width_preview)
            combo_flags |= ImGuiComboFlags_WidthFitPreview;
        if (h_small)
            combo_flags |= ImGuiComboFlags_HeightSmall;
        if (h_regular)
            combo_flags |= ImGuiComboFlags_HeightRegular;
        if (h_large)
            combo_flags |= ImGuiComboFlags_HeightLargest;

        const char *solver_items[] = {"analytic", "newton", "gradient"};

        const char *solver_preview_value = solver_items[p_params->solver];
        if (ImGui::BeginCombo("solver", solver_preview_value, combo_flags))
        {
            for (int n = 0; n < IM_ARRAYSIZE(solver_items); n++)
            {
                const bool is_selected = (p_params->solver == n);
                if (ImGui::Selectable(solver_items[n], is_selected))
                    p_params->solver = n;

                if (is_selected)
                    ImGui::SetItemDefaultFocus();
            }
            settings->update();
            ImGui::EndCombo();
        }
        if (ImGui::BeginItemTooltip())
        {
            ImGui::Text("This determines the way in which the surface points are calculated.\nUp to 3 atomic positions are taken into account.\n\nAnalytical: Uses algebraic relationships to determine the closest surface point. \nNumerical: Uses the Newton iteration to gradually approximate the position of the surface point.");
            ImGui::EndTooltip();
        }
        if (p_params->solver == 1 || p_params->solver == 2)
        {
            if (ImGui::InputInt("maximum solver iterations", &p_params->solver_iter_max, 1, 10))
            {
                if (p_params->solver_iter_max > 0)
                    settings->update();
            }
            if (ImGui::InputFloat("solver threshold", &p_params->solver_threshold, 0.005, 0.01))
            {
                if (p_params->solver_threshold > 0)
                    settings->update();
            }
        }

        // 3. k nearest (integer with step buttons)
        static int k_near = (int)p_params->k_nearest;
        if (ImGui::InputInt("k nearest", &k_near, 1, 10))
        {
            if (k_near > 0)
            {
                settings->setKnearest(k_near);
                settings->update();
            }
        }
        if (ImGui::BeginItemTooltip())
        {
            ImGui::Text("This parameter sets a maximum number of neighbouring atoms that are considered for the calculation.\nValues starting at 10 atoms are recommended.");
            ImGui::EndTooltip();
        }

        // 4. solvent radius (float with step buttons) -> value check
        if (ImGui::InputFloat("solvent radius", &p_params->solvent_radius, 0.1, 1))
        {
            if (p_params->solvent_radius > 0 && p_params->solvent_radius <= p_params->solvent_max)
            {
                settings->update();
            }
            else
            {
                if (p_params->solvent_radius <= 0)
                    p_params->solvent_radius = 0.01;
                if (p_params->solvent_radius > p_params->solvent_max)
                    p_params->solvent_radius = p_params->solvent_max;
            }
        }
        if (ImGui::BeginItemTooltip())
        {
            ImGui::Text("This is the radius of the probe molecule to yield the solvent accessible surface.\nIt's default radius is 1.4 Angstroms, equivalent to the radius of a hydrogen atom.\nThe probe radius can be changed in the range from 0 (exclusive) to 2 (inclusive) Angstrom.");
            ImGui::EndTooltip();
        }

        // 5. epsilon (float without step buttons)
        if (ImGui::InputFloat("epsilon", &p_params->epsilon))
        {
            if (p_params->epsilon > 0)
            {
                p_params->epsilon_squared = p_params->epsilon * p_params->epsilon;
                settings->update();
            }
        }
        if (ImGui::BeginItemTooltip())
        {
            ImGui::Text("This is the parameter to determine the precision of calculations. The smaller the epsilon, the more accurate the calculation.");
            ImGui::EndTooltip();
        }

        ImGui::SeparatorText("DISPLAY:");

        // 6. color scheme (combo box)
        const char *color_items[] = {"High Contrast", "Corey", "Koltun", "Jmol", "Rasmol (old)", "Rasmol (new)", "PubChem"};
        static int color_current_idx = (int)p_params->colorScheme;
        const char *color_preview_value = color_items[color_current_idx];

        if (ImGui::BeginCombo("color scheme", color_preview_value, combo_flags))
        {
            for (int n = 0; n < IM_ARRAYSIZE(color_items); n++)
            {
                const bool is_selected = (color_current_idx == n);
                if (ImGui::Selectable(color_items[n], is_selected))
                {
                    color_current_idx = n;
                    settings->setColorScheme(color_current_idx);
                    settings->update();
                    settings->loadMolecule(settings->getRecentFilePath());
                }

                if (is_selected)
                {
                    ImGui::SetItemDefaultFocus();
                }
            }
            ImGui::EndCombo();
        }
        if (ImGui::BeginItemTooltip())
        {
            ImGui::Text("This changes the color scheme of the program to fit popular color conventions for molecular models.\n\n IMPORTANT: This requires a reload of the molecule!");
            ImGui::EndTooltip();
        }

        if (ImGui::Checkbox("use color gradients", &p_params->use_interpolation))
            settings->update();

        // 7. window size for calculation (combo box with common resolutions)
        const char *resolution_items[] = {"1280 x 720 (16:9) [HD]", "1366 x 768", "1680 x 1050 (16:10)", "1920 x 1080 (16:9) [FHD]", "2560 x 1440 (16:9) [QHD]", "3440 x 1440 (21:9) [UWQHD]", "3840 x 2160 (16:9) [4K]"};
        const resolution res_pairs[] = {{1280, 720}, {1366, 768}, {1680, 1050}, {1920, 1080}, {2560, 1440}, {3440, 1440}, {3840, 2160}};
        static int resolution_current_idx = 3;
        const char *resolution_preview_value = resolution_items[resolution_current_idx];
        if (ImGui::BeginCombo("resolution", resolution_preview_value, combo_flags))
        {
            for (int n = 0; n < IM_ARRAYSIZE(resolution_items); n++)
            {
                const bool is_selected = (resolution_current_idx == n);
                if (ImGui::Selectable(resolution_items[n], is_selected))
                {
                    resolution_current_idx = n;
                    settings->setWindowWidth(res_pairs[resolution_current_idx].width);
                    settings->setWindowHeight(res_pairs[resolution_current_idx].height);
                    reset_textures();
                    createTextureBuffers(settings->getWindowWidth(), settings->getWindowHeight());
                    settings->update();
                }

                if (is_selected)
                    ImGui::SetItemDefaultFocus();
            }
            ImGui::EndCombo();
        }
        if (ImGui::BeginItemTooltip())
        {
            ImGui::Text("This sets the resolution that's used for the ray marching.");
            ImGui::Text("It defaults to 1920 x 1080 Pixels (Full HD resolution).");
            ImGui::EndTooltip();
        }
        // 8. limit framerate
        if (ImGui::BeginTable("framerate", 2))
        {
            ImGui::TableNextColumn();
            ImGui::Checkbox("limit framerate", &settings->limit_framerate);
            if (ImGui::BeginItemTooltip())
            {
                ImGui::Text("This allows limiting the maximum framerate.");
                ImGui::EndTooltip();
            }
            ImGui::TableNextColumn();
            if (settings->limit_framerate)
            {
                ImGui::SameLine();
                if (ImGui::InputInt(" ", &settings->frameLimit, 10, 20))
                {
                    if (settings->frameLimit < 1)
                        settings->frameLimit = 60;
                }
            }
            ImGui::EndTable();
        }

        ImGui::SeparatorText("ADVANCED:");

        // 9. TODO: overwrite block/thread dimensions (two fields, integer, step buttons)
        static int threads[2] = {(int)p_params->thread_x, (int)p_params->thread_y};
        if (ImGui::InputInt2("#threads", threads))
        {
            if (threads[0] > 0 && threads[1] > 0 && threads[0] <= MAXTHREADS && threads[1] <= MAXTHREADS)
            {
                p_params->thread_x = threads[0];
                p_params->thread_y = threads[1];
                settings->update();
            }
            else
            {
                if (threads[0] > 0 && threads[1] > 0)
                {
                    printf("WARNING: thread count must not exceed 512.");
                    threads[0] = threads[1] = 512;
                }
                else
                {
                    printf("WARNING: lowest thread number is 1.");
                    threads[0] = threads[1] = 1;
                }
            }
        }

        if (ImGui::BeginItemTooltip())
        {
            ImGui::Text("This overwrites the standard thread distribution in CUDA.");
            ImGui::EndTooltip();
        }
        // 10. debugmode (checkbox)
        if (ImGui::Checkbox("debugmode", &p_params->debug_mode))
        {
            if (p_params->debug_mode)
            {
                settings->setDebugFrame(0);
                printf("INFO: entering debug mode.\n");
            }
            else
            {
                settings->setDebugFrame(-1);
                printf("INFO: leaving debug mode.\n");
            }
            settings->update();
        }
        if (ImGui::BeginItemTooltip())
        {
            ImGui::Text("This settings is for debugging only. It enables additional debug information (e.g. counting frames, frames since debugmode was activated).");
            ImGui::EndTooltip();
        }
    }
    ImGui::Dummy(ImVec2(0.0f, 20.0f));
    if (ImGui::CollapsingHeader("Controls", header_flags))
    {
        ImGui::SeparatorText("KEYBINDINGS:");
        // TODO: table with keybindings
        if (ImGui::BeginTable("Molecule Information", 2, table_flags))
        {
            ImGui::TableNextColumn();
            ImGui::Text("Camera movement FORWARD");
            ImGui::TableNextColumn();
            ImGui::Text("Page Up");

            ImGui::TableNextColumn();
            ImGui::Text("Camera movement BACKWARDS");
            ImGui::TableNextColumn();
            ImGui::Text("Page Down");

            ImGui::TableNextColumn();
            ImGui::Text("Camera movement LEFT");
            ImGui::TableNextColumn();
            ImGui::Text("Arrow Left");

            ImGui::TableNextColumn();
            ImGui::Text("Camera movement RIGHT");
            ImGui::TableNextColumn();
            ImGui::Text("Arrow Right");

            ImGui::TableNextColumn();
            ImGui::Text("Camera movement UP");
            ImGui::TableNextColumn();
            ImGui::Text("Arrow Up");

            ImGui::TableNextColumn();
            ImGui::Text("Camera movement DOWN");
            ImGui::TableNextColumn();
            ImGui::Text("Arrow Down");

            ImGui::TableNextColumn();
            ImGui::Text("Reset camera");
            ImGui::TableNextColumn();
            ImGui::Text("R");

            ImGui::TableNextColumn();
            ImGui::Text("Camera Zoom In");
            ImGui::TableNextColumn();
            ImGui::Text("Scroll Wheel Up");

            ImGui::TableNextColumn();
            ImGui::Text("Camera Zoom Out");
            ImGui::TableNextColumn();
            ImGui::Text("Scroll Wheel Down");

            ImGui::TableNextColumn();
            ImGui::Text("Show/Hide User Interface");
            ImGui::TableNextColumn();
            ImGui::Text("I");

            ImGui::TableNextColumn();
            ImGui::Text("Show/Hide Performance");
            ImGui::TableNextColumn();
            ImGui::Text("P");

            ImGui::EndTable();
        }
    }
    ImGui::Dummy(ImVec2(0.0f, 20.0f));
    if (ImGui::CollapsingHeader("Help", header_flags))
    {
        // TODO: add github link
        ImGui::Text("GITHUB LINK:");
        ImGui::Text("............................");
        // ImGui::SeparatorText("IMGUI USER GUIDE:");
        // ImGui::ShowUserGuide();
    }

    // End of gui window
    ImGui::PopItemWidth();
    ImGui::End();
}
// imgui
void setupMetricsTool(AppSettings *settings, std::chrono::system_clock::time_point start_time)
{
    IM_ASSERT(ImGui::GetCurrentContext() != NULL && "ERROR: Missing Dear ImGui context.\n");

    ImGuiWindowFlags metric_window_flags = 0;
    metric_window_flags |= ImGuiWindowFlags_NoTitleBar;
    metric_window_flags |= ImGuiWindowFlags_NoScrollbar;
    // metric_window_flags |= ImGuiWindowFlags_NoMove;
    metric_window_flags |= ImGuiWindowFlags_NoResize;
    metric_window_flags |= ImGuiWindowFlags_NoCollapse;
    metric_window_flags |= ImGuiWindowFlags_NoNav;
    metric_window_flags |= ImGuiWindowFlags_AlwaysAutoResize;
    // metric_window_flags |= ImGuiWindowFlags_NoBackground;
    // metric_window_flags |= ImGuiWindowFlags_NoBringToFrontOnFocus;

    // size
    const ImGuiViewport *metric_viewport = ImGui::GetMainViewport();

    ImGui::SetNextWindowPos(ImVec2(metric_viewport->Size.x - 300, 0), ImGuiCond_FirstUseEver);
    ImGui::SetNextWindowSize(ImVec2(0, 0), ImGuiCond_FirstUseEver);
    if (!ImGui::Begin("metrics tools", &settings->show_performance_tool, metric_window_flags))
    {
        ImGui::End();
        return;
    }
    ImGui::PushItemWidth(ImGui::GetFontSize() * -16);

    if (ImGui::CollapsingHeader("Statistics", ImGuiTreeNodeFlags_DefaultOpen))
    {
        ImGui::SetWindowPos(ImVec2(metric_viewport->Size.x - ImGui::GetWindowSize().x, 0), true);
        float fps = ImGui::GetIO().Framerate;
        ImGui::Text("FPS: %.2f", fps);
        if (ImGui::BeginItemTooltip())
        {
            ImGui::Text("fps from rolling average over last 120 frames");
            ImGui::EndTooltip();
        }

        std::chrono::duration<double> time = std::chrono::duration_cast<std::chrono::seconds>(std::chrono::system_clock::now() - start_time);

        std::chrono::hours hrs = std::chrono::duration_cast<std::chrono::hours>(time);
        std::chrono::minutes mins = std::chrono::duration_cast<std::chrono::minutes>(time - hrs);
        std::chrono::seconds secs = std::chrono::duration_cast<std::chrono::seconds>(time - hrs - mins);

        ImGui::Text("Application time: %02i:%02i:%02i", hrs.count(), mins.count(), secs.count());
        ImGui::Text("Application frame: %i", settings->frame);
    }
    ImGui::SetWindowPos(ImVec2(metric_viewport->Size.x - ImGui::GetWindowSize().x, 0), true);

    if (ImGui::CollapsingHeader("FPS graph", ImGuiTreeNodeFlags_DefaultOpen))
    {
        ImGui::SetWindowPos(ImVec2(metric_viewport->Size.x - ImGui::GetWindowSize().x, 0), true);
        if (ImGui::BeginTable("  ", 2))
        {
            ImGui::TableNextColumn();
            ImGui::Dummy(ImVec2(0.0f, 3.0f));
            const int PLOTTED_FRAMES = 100;
            static float frames[2 * PLOTTED_FRAMES] = {};
            // static int plot_offset = 0;
            static int frame_offset1 = 0;
            static int frame_offset2 = 0;
            static int plot_offset = 0;

            float fps = 1.0f / ImGui::GetIO().DeltaTime;
            frames[frame_offset1] = fps;
            frames[frame_offset2] = fps;
            frame_offset1 = (frame_offset1 + 1) % IM_ARRAYSIZE(frames);
            frame_offset2 = (frame_offset1 + PLOTTED_FRAMES) % IM_ARRAYSIZE(frames);

            int baseline = ImGui::GetIO().Framerate;
            int len = log10((int)baseline);
            float div = pow(10, len);
            float rounded = ceil(baseline / div) * div;
            plot_offset = min(frame_offset1, frame_offset2);
            ImGui::PlotLines(" ", &frames[plot_offset], (IM_ARRAYSIZE(frames) / 2), 0, NULL, 0.5 * rounded, 1.5 * rounded, ImVec2(250.0f, 80.0f));

            ImGui::TableNextColumn();
            ImGui::Text("%.2f", 1.5 * rounded);
            ImGui::Dummy(ImVec2(0.0f, 16.5f));
            ImGui::Text("%.2f", fps);
            if (ImGui::BeginItemTooltip())
            {
                ImGui::Text("current fps");
                ImGui::EndTooltip();
            }
            ImGui::Dummy(ImVec2(0.0f, 16.5f));
            ImGui::Text("%.2f", 0.5 * rounded);

            ImGui::EndTable();
        }
    }

    // End of gui window
    ImGui::PopItemWidth();
    ImGui::End();
}

#endif