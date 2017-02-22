#include <igl/viewer/Viewer.h>
#include "scene.h"


Eigen::MatrixXd V;
Eigen::MatrixXi F;

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

bool key_down(igl::viewer::Viewer &viewer, unsigned char key, int modifier)
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


bool pre_draw(igl::viewer::Viewer &viewer)
{
    using namespace Eigen;
    using namespace std;
    
    if (viewer.core.is_animating){
        scene.updateScene(timeStep, CRCoeff, V,F);
        currTime+=timeStep;
        //cout <<"currTime: "<<currTime<<endl;
    }
    viewer.data.clear();
    viewer.data.set_mesh(V,F);

    return false;
}



int main(int argc, char *argv[])
{
    using namespace Eigen;
    using namespace std;

    
    // Load scene
    if (argc<2){
        cout<<"Please provide name of scene file!"<<endl;
        return 0;
    }
    cout<<"scene file: "<<std::string(DATA_PATH)+std::string("/")+std::string(argv[1])<<endl;
    scene.loadScene(std::string(DATA_PATH)+std::string("/")+std::string(argv[1]));
    
    
    //create platform
    MatrixXd platV;
    MatrixXi platF;
    RowVector3d platCOM;
    RowVector4d platOrientation;
    createPlatform(platV, platF, platCOM, platOrientation);
    
    scene.addRigidObject(platV, platF, 10000.0, true, platCOM, platOrientation);
    scene.updateScene(0.0, CRCoeff, V,F);
    
    // Viewer Settings
    igl::viewer::Viewer viewer;
    viewer.data.set_mesh(V,F);
    viewer.callback_pre_draw = &pre_draw;
    viewer.callback_key_down = &key_down;
    viewer.core.is_animating = false;
    viewer.core.animation_max_fps = 50.;
    
    
    //Adding options to GUI
    //To add new options, just modify the code below in a similar manner.
    viewer.callback_init = [&](igl::viewer::Viewer& viewer)
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
