#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <imgui/imgui.h>
#include <iostream>
#include "scene.h"


Eigen::MatrixXd V;
Eigen::MatrixXi F;
igl::opengl::glfw::Viewer viewer;

double currTime = 0;

//initial values
double timeStep = 0.02;
double CRCoeff= 1.0;

Scene scene;

void createPlatform(Eigen::MatrixXd& platV, Eigen::MatrixXi& platF, Eigen::RowVector3d& platCOM, Eigen::RowVector4d& platOrientation)
{
    double platWidth=100.0;
    platCOM<<0.0,-5.0,-0.0;
    platV.resize(8,3);
    platF.resize(12,3);
    platV<<-platWidth,0.0,-platWidth,
            -platWidth,0.0,platWidth,
            platWidth,0.0,platWidth,
            platWidth,0.0, -platWidth,
            -platWidth,-platWidth/10.0,-platWidth,
            -platWidth,-platWidth/10.0,platWidth,
            platWidth,-platWidth/10.0,platWidth,
            platWidth,-platWidth/10.0, -platWidth;
    platF<<0,1,2,
            2,3,0,
            6,5,4,
            4,7,6,
            1,0,5,
            0,4,5,
            2,1,6,
            1,5,6,
            3,2,7,
            2,6,7,
            0,3,4,
            3,7,4;
    
    platOrientation<<1.0,0.0,0.0,0.0;
}

bool key_down(igl::opengl::glfw::Viewer &viewer, unsigned char key, int modifier)
{
    if (key == ' ')
    {
        viewer.core.is_animating = !viewer.core.is_animating;
        return true;
    }
    
    if (key == 'S')
    {
        if (!viewer.core.is_animating){
            scene.updateScene(timeStep, CRCoeff, V,F);
            currTime+=timeStep;
            std::cout <<"currTime: "<<currTime<<std::endl;
            return true;
        }
    }
    return false;
}


bool pre_draw(igl::opengl::glfw::Viewer &viewer)
{
    using namespace Eigen;
    using namespace std;
    
    if (viewer.core().is_animating){
        scene.updateScene(timeStep, CRCoeff, V,F);
        currTime+=timeStep;
        //cout <<"currTime: "<<currTime<<endl;
    }
    viewer.data().clear();
    viewer.data().set_mesh(V,F);

    return false;
}

class CustomMenu : public igl::opengl::glfw::imgui::ImGuiMenu
{
  float floatVariable = 0.1f; // Shared between two menus
  
  virtual void draw_viewer_menu() override
  {
    // Draw parent menu
    ImGuiMenu::draw_viewer_menu();
    
    // Add new group
    if (ImGui::CollapsingHeader("New Group", ImGuiTreeNodeFlags_DefaultOpen))
    {
      // Expose variable directly ...
      ImGui::InputFloat("float", &floatVariable, 0, 0, 3);
      
      // ... or using a custom callback
      static bool boolVariable = true;
      if (ImGui::Checkbox("bool", &boolVariable))
      {
        // do something
        std::cout << "boolVariable: " << std::boolalpha << boolVariable << std::endl;
      }
      
      // Expose an enumeration type
      enum Orientation { Up=0, Down, Left, Right };
      static Orientation dir = Up;
      ImGui::Combo("Direction", (int *)(&dir), "Up\0Down\0Left\0Right\0\0");
      
      // We can also use a std::vector<std::string> defined dynamically
      static int num_choices = 3;
      static std::vector<std::string> choices;
      static int idx_choice = 0;
      if (ImGui::InputInt("Num letters", &num_choices))
      {
        num_choices = std::max(1, std::min(26, num_choices));
      }
      if (num_choices != (int) choices.size())
      {
        choices.resize(num_choices);
        for (int i = 0; i < num_choices; ++i)
          choices[i] = std::string(1, 'A' + i);
        if (idx_choice >= num_choices)
          idx_choice = num_choices - 1;
      }
      ImGui::Combo("Letter", &idx_choice, choices);
      
      // Add a button
      if (ImGui::Button("Print Hello", ImVec2(-1,0)))
      {
        std::cout << "Hello\n";
      }
    }
  }
  
  virtual void draw_custom_window() override
  {
    // Define next window position + size
    ImGui::SetNextWindowPos(ImVec2(180.f * menu_scaling(), 10), ImGuiSetCond_FirstUseEver);
    ImGui::SetNextWindowSize(ImVec2(200, 160), ImGuiSetCond_FirstUseEver);
    ImGui::Begin(
                 "New Window", nullptr,
                 ImGuiWindowFlags_NoSavedSettings
                 );
    
    // Expose the same variable directly ...
    ImGui::PushItemWidth(-80);
    ImGui::DragFloat("float", &floatVariable, 0.0, 0.0, 3.0);
    ImGui::PopItemWidth();
    
    static std::string str = "bunny";
    ImGui::InputText("Name", str);
    
    ImGui::End();
  }
  
};



int main(int argc, char *argv[])
{
    using namespace Eigen;
    using namespace std;

    
    // Load scene
    if (argc<2){
        cout<<"Please provide name of scene file!"<<endl;
        return 0;
    }
    cout<<"scene file: "<<std::string("../data/")+std::string(argv[1])<<endl;
    scene.loadScene(std::string("../data/")+std::string(argv[1]));
    
    
    //create platform
    MatrixXd platV;
    MatrixXi platF;
    RowVector3d platCOM;
    RowVector4d platOrientation;
    createPlatform(platV, platF, platCOM, platOrientation);
    
    scene.addRigidObject(platV, platF, 10000.0, true, platCOM, platOrientation);
    scene.updateScene(0.0, CRCoeff, V,F);
    
    // Viewer Settings
    viewer.data().set_mesh(V,F);
    viewer.callback_pre_draw = &pre_draw;
    viewer.callback_key_down = &key_down;
    viewer.core.is_animating = false;
    viewer.core.animation_max_fps = 50.;
    
    
    //Adding options to GUI
    //To add new options, just modify the code below in a similar manner.
    viewer.callback_init = [&](igl::opengl::glfw::Viewer& viewer)
    {
        //algorithmic options
        viewer.ngui->addGroup("Simulation");
        viewer.ngui->addVariable("CR Coeff",CRCoeff);
        viewer.ngui->addVariable<double>("Time Step",[&](double val) {
            viewer.core.animation_max_fps = (((int)1.0/val));
            timeStep=val;
        },[&]() {
            return timeStep;
        });
        
        // call to generate menu
        viewer.screen->performLayout();
        return false;
    };
    
    
    cout<<"Press [space] to toggle continuous simulation" << endl;
    cout<<"Press 'S' to advance time step-by-step"<<endl;
    viewer.launch();
}
