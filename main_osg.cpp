//#include <QCoreApplication>

#include <iostream>

#define GL_SILENCE_DEPRECATION

#include <osgViewer/Viewer>
#include <osgViewer/config/SingleWindow>

#include <imgui.h>
#include <imgui_impl_opengl3.h>

#include "OsgImGuiHandler.h"

#include "voxelGenerator/terraingenerator_roblox.h"

class ImGuiInitOperation : public osg::Operation
{
public:
    ImGuiInitOperation()
        : osg::Operation("ImGuiInitOperation", false)
    {
    }

    void operator()(osg::Object* object) override
    {
        osg::GraphicsContext* context = dynamic_cast<osg::GraphicsContext*>(object);
        if (!context)
            return;

        if (!ImGui_ImplOpenGL3_Init("#version 120"))
        {
            std::cout << "ImGui_ImplOpenGL3_Init() failed\n";
        }
    }
};

class ImGuiDemo : public OsgImGuiHandler
{
protected:
    void drawUi() override
    {
        // ImGui code goes here...
        //ImGui::ShowDemoWindow();
        ImGui::Begin("Hello, world!");
        if (ImGui::Button("Open File")) {
            ImGui::InputInt3("size", size_terrain);
            if (ImGui::Button("Create Terrain")) {
                TerrainGenerator_Roblox::getInstance()->setRange(Vector3(0), Vector3(size_terrain[0], size_terrain[1], size_terrain[2]));
            }
        }

        ImGui::End();
    }

private:
    int size_terrain[3] = {20, 20, 20};
};

int main() {
    osgViewer::Viewer viewer;
    viewer.apply(new osgViewer::SingleWindow(100, 100, 640, 480));
    viewer.setRealizeOperation(new ImGuiInitOperation);
    viewer.addEventHandler(new ImGuiDemo);
    return viewer.run();
}
