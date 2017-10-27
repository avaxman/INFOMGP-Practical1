#ifndef SCENE_HEADER_FILE
#define SCENE_HEADER_FILE

#include <vector>
#include <fstream>
#include <igl/bounding_box.h>
#include <igl/readOFF.h>
#include "ccd.h"
#include "volInt.h"
#include "auxfunctions.h"

using namespace Eigen;
using namespace std;


void support(const void *_obj, const ccd_vec3_t *_d, ccd_vec3_t *_p);
void stub_dir(const void *obj1, const void *obj2, ccd_vec3_t *dir);
void center(const void *_obj, ccd_vec3_t *dir);



//Impulse is defined as a pair <position, direction>
typedef std::pair<RowVector3d,RowVector3d> Impulse;


//the class the contains each individual rigid objects and their functionality
class RigidObject{
public:
    MatrixXd origV;   //original vertex positions, where COM=(0.0,0.0,0.0) - never change this!
    MatrixXd currV;   //current vertex position
    MatrixXi T;

    //position of object in space. We must always have that currV = QRot(origV, orientation)+ COM
    RowVector4d orientation; //current orientation
    RowVector3d COM;  //current center of mass
    Matrix3d invIT;  //Original *inverse* inertia tensor around the COM, defined in the rest state to the object (so to the canonical world system)

    //kinematics
    bool isFixed;  //is the object immobile
    double mass;
    double currDensity;
    RowVector3d comVelocity;  //the linear velocity of the center of mass
    RowVector3d angVelocity;  //the angular velocity of the object.
    RowVector3d frictionTan; // the tangent used for friction calculations
    RowVector3d totalVelocity;
    RowVector3d gravity;
    RowVector3d drag; //drag force

    //dynamics
    std::vector<Impulse> impulses;  //Gets updated by the collision process


    //checking collision between bounding boxes, and consequently the rigid objects if succeeds.
    //you do not need to update these functions (isBoxCollide and isCollide) unless you are doing a different collision

    bool isBoxCollide(const RigidObject& ro){
        RowVector3d VMin1=currV.colwise().minCoeff();
        RowVector3d VMax1=currV.colwise().maxCoeff();
        RowVector3d VMin2=ro.currV.colwise().minCoeff();
        RowVector3d VMax2=ro.currV.colwise().maxCoeff();

        //checking all axes for non-intersection of the dimensional interval
        for (int i=0;i<3;i++)
            if ((VMax1(i)<VMin2(i))||(VMax2(i)<VMin1(i)))
                return false;

        return true;  //all dimensional intervals are overlapping = intersection

    }

    bool isCollide(const RigidObject& ro, double& depth, RowVector3d& intNormal, RowVector3d& intPosition){


        //collision between bounding boxes
        if (!isBoxCollide(ro))
            return false;


        ccd_t ccd;
        ccd_vec3_t sep;

		CCD_INIT(&ccd);
        //sophisticated collision between convex triangle meshes
        ccd.support1       = support; // support function for first object
        ccd.support2       = support; // support function for second object
        ccd.center1         =center;
        ccd.center2         =center;

        ccd.first_dir       = stub_dir;
        ccd.max_iterations = 100;     // maximal number of iterations

        void* obj1=(void*)this;
        void* obj2=(void*)&ro;

        ccd_real_t _depth;
        ccd_vec3_t dir, pos;

       // int nonintersect = ccdGJKPenetration(obj1, obj2,&ccd, &_depth,&dir,&pos);
        int nonintersect = ccdMPRPenetration(obj1, obj2, &ccd, &_depth, &dir, &pos);
        if (nonintersect)
            return false;


        for (int i=0;i<3;i++){
            intNormal(i)=dir.v[i];
            intPosition(i)=pos.v[i];
        }

        depth=_depth;
        intPosition-=depth*intNormal/2.0;  //to bring it to (current) obj2
        return !nonintersect;
    }


    //return the current inverted inertia tensor around the current COM. Update it by applying the orientation
    Matrix3d getCurrInvInertiaTensor(){
        Matrix3d R;
        /***************
         TODO
         ***************/
         //Complete
         getCOMandInvIT(currV, T, currDensity, mass, COM, R);
         for(int i = 0; i < R.rows()-1; i++){
            R.row(i) = QRot(R.row(i), orientation);
         }
        return R;
    }

    //Update the current position and orientation by integrating the linear and angular velocities, and update currV accordingly
    //You need to modify this according to its purpose
    void updatePosition(double timeStep){
        /***************
         TODO
         ***************/
         //Complete
         //integrate velocity
         //COM = (comVelocity * timeStep) + COM;
         COM = (comVelocity * timeStep) + COM;
         //integrate orientation
         RowVector4d w;
         w <<0, (angVelocity.x() * timeStep), (angVelocity.y()*timeStep), (angVelocity.z()*timeStep);
         RowVector4d orientationStep = QExp(w);
         //update orientation
         orientation = QMult(orientation, orientationStep);
         //update currV
         for(int i = 0; i < origV.rows(); i++){
           currV.row(i) = QRot(origV.row(i), orientation) + COM;
         }
         invIT = getCurrInvInertiaTensor();
    }

    //Updating velocity *instantaneously*. i.e., not integration from acceleration, but as a result of a collision impulse from the "impulses" list
    //You need to modify this for that purpose.
    void updateImpulseVelocities(){

              if (isFixed){
                  comVelocity.setZero();
                  impulses.clear();
                  angVelocity.setZero();
                  return;
              }

              //update linear and angular velocity according to all impulses
              for (int i=0;i<impulses.size();i++){
              RowVector3d pos = std::get<0>(impulses[i]);
              RowVector3d dir = std::get<1>(impulses[i]);
              double dirMag = dir.norm();
              RowVector3d angImp = invIT*pos.cross(dir + frictionTan*dirMag).transpose();
              RowVector3d linImp = (dir + frictionTan*dirMag)/(mass/1000);

                  /***************
                   TODO
                   ***************/
              for(int j=0; j<3; j++){
              comVelocity[j] += linImp[j];

            //  angVelocity[j] += angImp[j];
              }

              }
              impulses.clear();
    }


    //Updating the linear and angular velocities of the object
    //You need to modify this to integrate from acceleration in the field (basically gravity)
    void updateVelocity(double timeStep, double airDrag){

        if (isFixed)
            return;

        /***************
         TODO
         ***************/

         comVelocity = comVelocity + (gravity*timeStep);
         comVelocity = comVelocity - (airDrag*comVelocity);
    }


    //the full integration for the time step (velocity + position)
    //You need to modify this if you are changing the integration
    void integrate(double timeStep, double airDrag){
        updateVelocity(timeStep, airDrag);
        updatePosition(timeStep);
    }


    RigidObject(const MatrixXd& _V, const MatrixXi& _T, const double density, const bool _isFixed, const RowVector3d& _COM, const RowVector4d& _orientation){
        origV=_V;
        T=_T;
        isFixed=_isFixed;
        COM=_COM;
        currDensity = density;
        orientation=_orientation;
		    comVelocity.setZero();
        angVelocity.setZero();

        RowVector3d naturalCOM;  //by the geometry of the object

        //initializes the original geometry (COM + IT) of the object
        getCOMandInvIT(origV, T, density, mass, naturalCOM, invIT);

        gravity << 0, -9.81, 0;

        origV.rowwise()-=naturalCOM;  //removing the natural COM of the OFF file (natural COM is never used again)

        currV.resize(origV.rows(), origV.cols());
        for (int i=0;i<currV.rows();i++)
            currV.row(i)<<QRot(origV.row(i), orientation)+COM;
    }

    ~RigidObject(){}
};






//This class contains the entire scene operations, and the engine time loop.
class Scene{
public:
    double currTime;
    int numFullV, numFullT;
    std::vector<RigidObject> rigidObjects;


    //adding an objects. You do not need to update this generally
    void addRigidObject(const MatrixXd& V, const MatrixXi& T, const double density, const bool isFixed, const RowVector3d& COM, const RowVector4d orientation){

        RigidObject ro(V,T, density, isFixed, COM, orientation);
        rigidObjects.push_back(ro);
        numFullV+=V.rows();
        numFullT+=T.rows();
    }

    /*********************************************************************
     This function handles a collision between objects ro1 and ro2 when found, by assigning impulses to both objects.
     Input: RigidObjects ro1, ro2
            depth: the depth of penetration
            contactNormal: the normal of the conact measured ro1->ro2
            penPosition: a point on ro2 such that if ro2 <= ro2 + depth*contactNormal, then penPosition+depth*contactNormal is the common contact point
            CRCoeff: the coefficient of restitution
     *********************************************************************/
    void handleCollision(RigidObject& ro1, RigidObject& ro2, const double depth, const RowVector3d& contactNormal, const RowVector3d& penPosition, const double CRCoeff, double friction){
      //Interpretation resolution: move each object by inverse mass weighting, unless either is fixed, and then move the other. Remember to respect the direction of contactNormal and update penPosition accordingly.
      RowVector3d contactPosition;
      /***************
       TODO
       ***************/
       RowVector3d depthVector = depth*contactNormal;
   contactPosition = penPosition + depthVector;
   if(ro1.isFixed){
       ro2.COM -= depthVector;
   }
   else if(ro2.isFixed){
       ro1.COM += depthVector;
   }
   else{
           ro1.COM -= contactNormal*depth*0.5;

           ro2.COM += contactNormal*depth*0.5;
   }
   //Create impulses and push them into ro1.impulses and ro2.impulses.
   /***************
    TODO
    ***************/
   RowVector3d ro1R = contactPosition - ro1.COM;
   RowVector3d ro2R = contactPosition - ro2.COM;
   //RowVector3d ro1RCrossNorm = ro1R.cross(contactNormal);
   //RowVector3d ro2RCrossNorm = ro2R.cross(contactNormal);
   double impulseMag = -(1+CRCoeff)*
                       ((ro1.comVelocity - ro2.comVelocity).dot(contactNormal)
                       /((1/(ro1.mass/1000) + 1/(ro2.mass/1000))));
                    //   +(ro1R.cross(contactNormal)*ro1.invIT*ro1R.cross(contactNormal).transpose())+(ro2R.cross(contactNormal)*ro2.invIT*ro2R.cross(contactNormal).transpose())));

   RowVector3d frictionTan = (contactNormal.cross(ro1.comVelocity - ro2.comVelocity)).cross(contactNormal)*friction;
   ro1.frictionTan = frictionTan;
   ro2.frictionTan = frictionTan;

   Impulse* ro1Impulse = new Impulse(ro1R, -impulseMag*contactNormal);
   Impulse* ro2Impulse = new Impulse(ro2R, impulseMag*contactNormal);
   ro1.impulses.push_back(*ro1Impulse);
   ro2.impulses.push_back(*ro2Impulse);
   delete ro1Impulse;
   delete ro2Impulse;

   //updating velocities according to impulses
   ro1.updateImpulseVelocities();
   ro2.updateImpulseVelocities();
    }



    /*********************************************************************
    This function handles a single time step by:
     1. Integrating velocities, positions, and orientations by the timeStep
     2. detecting and handling collisions with the coefficient of restitutation CRCoeff
     3. updating the visual scene in fullV and fullT
     *********************************************************************/
    void updateScene(double timeStep, double CRCoeff, MatrixXd& fullV, MatrixXi& fullT, double airDrag, double friction){
        fullV.conservativeResize(numFullV,3);
        fullT.conservativeResize(numFullT,3);
        int currVIndex=0, currFIndex=0;

        //integrating velocity, position and orientation from forces and previous states
        for (int i=0;i<rigidObjects.size();i++)
            rigidObjects[i].integrate(timeStep, airDrag);

        //detecting and handling collisions when found
        //This is done exhaustively: checking every two objects in the scene.
        double depth;
        RowVector3d contactNormal, penPosition;
        for (int i=0;i<rigidObjects.size();i++)
            for (int j=i+1;j<rigidObjects.size();j++)
                if (rigidObjects[i].isCollide(rigidObjects[j],depth, contactNormal, penPosition))
                    handleCollision(rigidObjects[i], rigidObjects[j],depth, contactNormal.normalized(), penPosition,CRCoeff,friction);



        //Code for updating visualization meshes
        for (int i=0;i<rigidObjects.size();i++){
            fullT.block(currFIndex, 0, rigidObjects[i].T.rows(),3)=rigidObjects[i].T.array()+currVIndex;   //need to advance the indices, because every object is indexed independently
            fullV.block(currVIndex, 0, rigidObjects[i].currV.rows(),3)=rigidObjects[i].currV;
            currFIndex+=rigidObjects[i].T.rows();
            currVIndex+=rigidObjects[i].currV.rows();
        }
        currTime+=timeStep;
    }

    //loading a scene from the scene .txt files
    //you do not need to update this function
    bool loadScene(const std::string sceneFileName){

        ifstream sceneFileHandle;
        sceneFileHandle.open(sceneFileName);
        if (!sceneFileHandle.is_open())
            return false;
        int numofObjects;

        currTime=0;
        numFullT=0;
        numFullV=0;
        sceneFileHandle>>numofObjects;
        for (int i=0;i<numofObjects;i++){
            MatrixXi objT;
            MatrixXd objV;
            std::string OFFFileName;
            bool isFixed;
            double density;
            RowVector3d COM;
            RowVector4d orientation;
            sceneFileHandle>>OFFFileName>>density>>isFixed>>COM(0)>>COM(1)>>COM(2)>>orientation(0)>>orientation(1)>>orientation(2)>>orientation(3);
            orientation.normalize();
            igl::readOFF(std::string(DATA_PATH)+std::string("/")+OFFFileName,objV,objT);
            addRigidObject(objV,objT,density, isFixed, COM, orientation);
			cout << "COM: " << COM <<endl;
            cout << "orientation: " << orientation <<endl;
        }
        return true;
    }


    Scene(){}
    ~Scene(){}
};


/*****************************Auxiliary functions for collision detection. Do not need updating********************************/

/** Support function for libccd*/
void support(const void *_obj, const ccd_vec3_t *_d, ccd_vec3_t *_p)
{
    // assume that obj_t is user-defined structure that holds info about
    // object (in this case box: x, y, z, pos, quat - dimensions of box,
    // position and rotation)
    //std::cout<<"calling support"<<std::endl;
    RigidObject *obj = (RigidObject *)_obj;
    RowVector3d p;
    RowVector3d d;
    for (int i=0;i<3;i++)
        d(i)=_d->v[i]; //p(i)=_p->v[i];


    d.normalize();
    //std::cout<<"d: "<<d<<std::endl;

    int maxVertex=-1;
    int maxDotProd=-32767.0;
    for (int i=0;i<obj->currV.rows();i++){
        double currDotProd=d.dot(obj->currV.row(i)-obj->COM);
        if (maxDotProd < currDotProd){
            maxDotProd=currDotProd;
            //std::cout<<"maxDotProd: "<<maxDotProd<<std::endl;
            maxVertex=i;
        }

    }
    //std::cout<<"maxVertex: "<<maxVertex<<std::endl;

    for (int i=0;i<3;i++)
        _p->v[i]=obj->currV(maxVertex,i);

    //std::cout<<"end support"<<std::endl;
}

void stub_dir(const void *obj1, const void *obj2, ccd_vec3_t *dir)
{
    dir->v[0]=1.0;
    dir->v[1]=0.0;
    dir->v[2]=0.0;
}

void center(const void *_obj,ccd_vec3_t *center)
{
    RigidObject *obj = (RigidObject *)_obj;
    for (int i=0;i<3;i++)
        center->v[i]=obj->COM(i);
}






#endif
