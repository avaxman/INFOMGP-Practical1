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
  RowVector3d comVelocity;  //the linear velocity of the center of mass
  RowVector3d angVelocity;  //the angular velocity of the object.
  
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
	R = Q2RotMatrix(orientation);
	R = R.transpose() * invIT * R;

    return R;
  }
  
  
  //Update the current position and orientation by integrating the linear and angular velocities, and update currV accordingly
  //You need to modify this according to its purpose
  void updatePosition(double timeStep){
    /***************
     TODO
     ***************/
	  RowVector4d temp(0, angVelocity[0], angVelocity[1], angVelocity[2]);
	  RowVector4d newOrient = orientation + (1/2 * timeStep * QMult(temp, orientation));
	  orientation = newOrient;

	  RowVector3d nextCOM = COM + comVelocity * timeStep;
	  COM = nextCOM;

	  for (int i = 0; i < currV.rows(); i++)
		currV.row(i) = QRot(origV.row(i), orientation) + COM;
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
      /***************
       TODO
       ***************/
		comVelocity += impulses[i].second * 1/mass;
		angVelocity += getCurrInvInertiaTensor() * (impulses[i].first - COM).cross(impulses[i].second).transpose();
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
	RowVector3d acceleration{0, -9.8, 0};

	RowVector3d nextAngVel = angVelocity + ((getCurrInvInertiaTensor() * RowVector3d(0,0,0).cross(acceleration * mass).transpose() * timeStep).transpose());
	angVelocity = nextAngVel;

	RowVector3d nextComVel = comVelocity + ((-airDrag * comVelocity + (acceleration * mass)) / mass) * timeStep;
	comVelocity = nextComVel;
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
    orientation=_orientation;
    comVelocity.setZero();
    angVelocity.setZero();
    
    RowVector3d naturalCOM;  //by the geometry of the object
    
    //initializes the original geometry (COM + IT) of the object
    getCOMandInvIT(origV, T, density, mass, naturalCOM, invIT);
    
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
  void handleCollision(RigidObject& ro1, RigidObject& ro2, const double depth, const RowVector3d& contactNormal, const RowVector3d& penPosition, const double CRCoeff, double frictionCoeff){
    
    //Interpretation resolution: move each object by inverse mass weighting, unless either is fixed, and then move the other. Remember to respect the direction of contactNormal and update penPosition accordingly.
    RowVector3d contactPosition = penPosition + depth * contactNormal;
    /***************
     TODO
     ***************/	
	RowVector3d Ff = contactNormal * frictionCoeff;

	if (ro1.isFixed)
		ro2.COM += depth * contactNormal - Ff;
	else if (ro2.isFixed)
		ro1.COM +=  -1 * depth * contactNormal - Ff;
	else
	{
		double invMassRatio = ro2.mass / ro1.mass;

		ro1.COM += -1 * (depth * invMassRatio) * contactNormal - Ff;
		ro2.COM += (depth * (1- invMassRatio)) * contactNormal - Ff;
	}

	for (int i = 0; i < ro1.currV.rows(); i++)
		ro1.currV.row(i) = QRot(ro1.origV.row(i), ro1.orientation) + ro1.COM;

	for (int i = 0; i < ro2.currV.rows(); i++)
		ro2.currV.row(i) = QRot(ro2.origV.row(i), ro2.orientation) + ro2.COM;

    //Create impulses and push them into ro1.impulses and ro2.impulses.
    
    /***************
     TODO
     ***************/
	//double linearImpMagn;
	//if (ro1.isFixed)
	//{
	//	linearImpMagn = -1 * (1 + CRCoeff) * (-ro2.comVelocity).dot(contactNormal) * ro2.mass;
	//}
	//else if (ro2.isFixed)
	//{
	//	linearImpMagn = -1 * (1 + CRCoeff) * (ro1.comVelocity).dot(contactNormal) * ro1.mass;
	//}
	//else
	//{
	//	linearImpMagn = -1 * (1 + CRCoeff)*(ro1.comVelocity - ro2.comVelocity).dot(contactNormal) * (1 / ro1.mass + 1 / ro2.mass);
	//}

	//ro1.impulses.push_back(Impulse(contactPosition, (linearImpMagn / ro1.mass) * contactNormal));
	//ro2.impulses.push_back(Impulse(contactPosition, -(linearImpMagn / ro2.mass) * contactNormal));
	double com1 = ro1.comVelocity.sum();
	double com2 = ro2.comVelocity.sum();

	if (com1 > 0 && com2 > 0)
	{
		double angularImpMagn;
		RowVector3d collisionArm1, collisionArm2;
		collisionArm1 = contactPosition - ro1.COM;
		collisionArm2 = contactPosition - ro2.COM;
		RowVector3d vel1, vel2, cross;
		vel1 = ro1.comVelocity + ro1.angVelocity.cross(collisionArm1);
		vel2 = ro2.comVelocity + ro2.angVelocity.cross(collisionArm2);
		cross = contactNormal.cross((vel1 - vel2).cross(contactNormal));
		cross.normalize();
		double numerator = (-1 * (1 + CRCoeff) * (vel1 - vel2).dot(contactNormal + cross));

		RowVector3d rA, rB;
		rA = penPosition - ro1.COM;
		rB = penPosition - ro2.COM;
		RowVector3d crossA, crossB;
		crossA = rA.cross(contactPosition + cross);
		crossB = rB.cross(contactPosition + cross);
		double denominator = (1 / ro1.mass + 1 / ro2.mass) + (rA * ro1.getCurrInvInertiaTensor() * crossA.transpose()) + (rB * ro2.getCurrInvInertiaTensor() * crossB.transpose());

		angularImpMagn = numerator / denominator;

		ro1.impulses.push_back(Impulse(contactPosition, angularImpMagn * (contactNormal + cross)));
		ro2.impulses.push_back(Impulse(contactPosition, -angularImpMagn * (contactNormal + cross)));
	}
    //updating velocities according to impulses	
    ro1.updateImpulseVelocities();
    ro2.updateImpulseVelocities();
  }
  
  
  
  /*********************************************************************
   This function handles a single time step by:
   1. Integrating velocities, positions, and orientations by the timeStep
   2. detecting and handling collisions with the coefficient of restitution CRCoeff
   3. updating the visual scene in fullV and fullT
   *********************************************************************/
  void updateScene(double timeStep, double CRCoeff, MatrixXd& fullV, MatrixXi& fullT, double airDrag, double frictionCoeff){
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
          handleCollision(rigidObjects[i], rigidObjects[j],depth, contactNormal.normalized(), penPosition,CRCoeff, frictionCoeff);
    
    
    
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
  bool loadScene(const std::string path, const std::string sceneFileName){
    
    ifstream sceneFileHandle;
    sceneFileHandle.open(path+std::string("/")+sceneFileName);
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
      igl::readOFF(path+std::string("/")+OFFFileName,objV,objT);
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