# Practical 1: Rigid-Body Simulation

The first exercise will be about simulating a moving and rotating rigid body. The objectives of the practical are:

1. Implement rigid velocity, position, and orientation integration in a discrete time iteration.
2. Implement impulse-based collision resolution.
3. Extend the framework with some chosen effects, such as friction, drag, or 

This is the repository for the skeleton on which you will build your first exercise. Using CMake allows you to work and submit your code in all platforms.

##Installation

The skeleton uses the following dependencies: [libigl](http://libigl.github.io/libigl/) and consequently [Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page), and [libccd](https://github.com/danfis/libccd). libigl viewer is using [nanogui](https://github.com/wjakob/nanogui) Everything is bundled as submodules (or just code) within the skeleton, and you do not have to take care of it. To get the library, use:

```bash
git clone --recursive https://github.com/avaxman/INFOMGP-Practical1.git
```

to compile the examples, go into the `practical1` folder and enter in a terminal:

```bash
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ../
make
```

In windows, you need to use [cmake-gui](https://cmake.org/runningcmake/). Pressing twice ``configure`` and then ``generate`` will generate a Visual Studio solution in which you can work. The active soution should be ``practical1_bin``.

##Using the dependencies

You do not need to utilize any dependency on your own or install anything other than the above. For the most part it is background or collision detection code, which is not a direct part of the exercise. The exception is [Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page) for the representation and manipulation of vectors and matrices, but it is a quite a shallow learning curve. Go through the getting started sectrion on the website (reading until and including "Dense matrix and array manipulation" should be enough).

##Working with the repository

All the code you need to update is in the ``practical1`` folder. Please do not attempt to commit any changes to here. <span style="color:red">You may ONLY fork the repository for your convenience and work on it if you can somehow make the forked repository PRIVATE afterwards</span>. Public open-source solutions to the exercise will disqualify the students! submission will be done in the "classical" department style of submission servers.

##The code

Most of the action happens in `scene.h`. The main function is:

```cpp
void updateScene(double timeStep, double CRCoeff, MatrixXd& fullV, MatrixXi& fullT){
        fullV.conservativeResize(numFullV,3);
        fullT.conservativeResize(numFullT,3);
        int currVIndex=0, currFIndex=0;
        
        //integrating velocity, position and orientation from forces and previous states
        for (int i=0;i<rigidObjects.size();i++)
            rigidObjects[i].integrate(timeStep);
            
        //detecting and handling collisions when found
        //This is done exhaustively: checking every two objects in the scene.
        double depth;
        RowVector3d contactNormal, penPosition;
        for (int i=0;i<rigidObjects.size();i++)
            for (int j=i+1;j<rigidObjects.size();j++)
                if (rigidObjects[i].isCollide(rigidObjects[j],depth, contactNormal, penPosition))
                    handleCollision(rigidObjects[i], rigidObjects[j],depth, contactNormal.normalized(), penPosition,CRCoeff);
        
        
        
        //Code for updating visualization meshes
        for (int i=0;i<rigidObjects.size();i++){
            fullT.block(currFIndex, 0, rigidObjects[i].T.rows(),3)=rigidObjects[i].T.array()+currVIndex;   //need to advance the indices, because every object is indexed independently
            fullV.block(currVIndex, 0, rigidObjects[i].currV.rows(),3)=rigidObjects[i].currV;
            currFIndex+=rigidObjects[i].T.rows();
            currVIndex+=rigidObjects[i].currV.rows();
        }
        currTime+=timeStep;
    } 
```
