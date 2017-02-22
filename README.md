# Practical 1: Rigid-Body Simulation

The first exercise will be about simulating a moving and rotating rigid body. The objectives of the practical are:

1. Implement rigid velocity, position, and orientation integration in a discrete time iteration.
2. Implement impulse-based collision resolution.
3. Extend the framework with some chosen effects, such as friction or drag.  

This is the repository for the skeleton on which you will build your first exercise. Using CMake allows you to work and submit your code in all platforms. The entire environment is in C++, but most of the "nastier" coding parts have been simplified; for the most part, you only need to work in the math.

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

##The coding environment

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

The two most important functions are ``integrate()`` and ``handleCollision()``. They all contain a mixture of written code, and code you have to complete. The code you have to complete is always marked as:

```cpp
/***************
TODO
***************/
```

For the most part, the description of the function will tell you what exactly you need to complete. See summary LINK for the detailed tasks 

###User interface

The program is loaded by giving a TXT file that describes the scene as an argument to the executable. The file should be in the `data` subfolder, which is automatically discovered by the CMake. The format of the file is:

```
#num_objects
object1.off     density1     is_fixed1    COM1     o1
object2.off     density2     is_fixed2    COM2     o2 
.....
```

Where:

1. objectX.off - an OFF file describing the geometry of the a triangulated mesh. [Meshlab](http://www.meshlab.net/) can view these files, but their format is pretty straightforward. The object is loaded, and translated to have a center of mass of $(0,0,0)$ in its object coordinates. 
2. density - the uniform density of the object. The program will automatically compute the mas by the volume.
3. is_fixed - if the object should be immobile (fixed in space) or not.
4. COM - the initial position in space where the object would be translated to. That means where the COM is in time $t=0$.
5. o - the initial orientation of the object, expressed as a quaternion that rotates the geometry to $o*object*inv(o)$ at time $t=0$.




