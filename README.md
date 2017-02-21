# INFOMGP - Skeleton for Practical 1 (Rigid-Body Simulation).

This is the repository for the skeleton on which you will build your first exercise. Using CMake allows you to work and submit your code in all platforms.

##Installation

The skeleton uses the following dependencies: [libigl](http://libigl.github.io/libigl/) and consequently [Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page), and [libccd](https://github.com/danfis/libccd). libigl viewer is using [nanogui](https://github.com/wjakob/nanogui) Everything is bundled as submodules (or just code) within the skeleton, and you do not have to take care of it. To get the library, use:

```bash
git clone --recursive https://github.com/avaxman/libhedra.git
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
