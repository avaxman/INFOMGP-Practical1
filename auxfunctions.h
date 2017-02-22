//
//  auxfunctions.h
//  practical1-test-project
//
//  Created by Amir Vaxman on 22/02/2017.
//  Copyright © 2017 Amir Vaxman. All rights reserved.
//

#ifndef auxfunctions_h
#define auxfunctions_h


/*****************************Quaternionic functions***********************************************************/


//returns the conjugate of a quaternion
inline Eigen::RowVector4d QConj(const Eigen::RowVector4d& q)
{
    Eigen::RowVector4d newq(q.rows(), q.cols());
    newq<<q(0), -q.tail(3);
    return newq;
}


//multiplies two quaternions
inline Eigen::RowVector4d QMult(const Eigen::RowVector4d& q1, const Eigen::RowVector4d& q2)
{
    Eigen::RowVector4d newq;
    double r1=q1(0);
    double r2=q2(0);
    Eigen::RowVector3d v1=q1.tail(3);
    Eigen::RowVector3d v2=q2.tail(3);
    newq<<r1*r2-v1.dot(v2), r1*v2+r2*v1+v1.cross(v2);
    return newq;
}


//The inverse of a quaternion
inline Eigen::RowVector4d QInv(const Eigen::RowVector4d& q)
{
    return(QConj(q)/q.squaredNorm());
}


//the exponent of a quaternion: from (0, theta*vec) (axis angle) to (cos(theta), sin(theta)*vec)) (quaternion representation of that rotation).
inline Eigen::RowVector4d QExp(const Eigen::RowVector4d& q)
{
    double nv=q.tail(3).norm();
    Eigen::RowVector4d expq;
    double expq0=exp(q(0));
    if (nv<10e-6){
        expq<<expq0,0.0,0.0,0.0;
        return expq;
    }
    
    expq<<cos(nv), q.tail(3).normalized()*sin(nv);
    expq*=expq0;
    
    return expq;
}


//Rotating an imaginary quaternion p (or a 3D vector) by q*p*q^{-1}
inline Eigen::RowVector3d QRot(const Eigen::RowVector3d p, const Eigen::RowVector4d q){
    Eigen::RowVector4d quatp; quatp<<0.0,p;
    Eigen::RowVector4d w=QMult(q, QMult(quatp, QInv(q)));
    return w.tail(3);
}


//returning a rotation matrix R for quaternion q so that R*vec.transpose() = QRot(vec, q).transpose() for a vec:1x3 (RowVector3d)
inline Eigen::Matrix3d Q2RotMatrix(const Eigen::RowVector4d q){
    Eigen::Matrix3d R;
    Eigen::Matrix3d corrv=q.tail(3).transpose()*q.tail(3);
    std::cout<<"q.norm(): "<<q.norm()<<std::endl;
    
    //diagonal elements
    R(0,0)=1.0-2.0*(corrv(1,1)+corrv(2,2));
    R(1,1)=1.0-2.0*(corrv(2,2)+corrv(0,0));
    R(2,2)=1.0-2.0*(corrv(0,0)+corrv(1,1));
    
    //off diagonal elements
    R(0,1)=2*(corrv(0,1)-q(0)*q(3));
    R(1,0)=2*(corrv(1,0)+q(0)*q(3));
    
    R(0,2)=2*(corrv(0,2)+q(0)*q(2));
    R(2,0)=2*(corrv(2,0)-q(0)*q(2));
    
    R(1,2)=2*(corrv(1,2)-q(0)*q(1));
    R(2,1)=2*(corrv(2,1)+q(0)*q(1));
    
    //testing procedure
    /*Eigen::RowVector3d rand=Eigen::RowVector3d::Random();
    if ((R*rand.transpose()-QRot(rand,q).transpose()).norm()>10e-6){
        std::cout<<"R*rand.transpose(): "<<R*rand.transpose()<<std::endl;
        std::cout<<"QRot(rand,q).transpose(): "<<QRot(rand,q).transpose()<<std::endl;
    }*/
    
    
    return R;
}


#endif /* auxfunctions_h */
