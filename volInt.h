

	/*******************************************************
        *                                                      *
	*  volInt.c                                            *
	*                                                      *
	*  This code computes volume integrals needed for      *
	*  determining mass properties of polyhedral bodies.   *
	*                                                      *
	*  For more information, see the accompanying README   *
	*  file, and the paper                                 *
	*                                                      *
	*  Brian Mirtich, "Fast and Accurate Computation of    *
	*  Polyhedral Mass Properties," journal of graphics    *
	*  tools, volume 1, number 1, 1996.                    *
	*                                                      *
	*  This source code is public domain, and may be used  *
	*  in any way, shape or form, free of charge.          *
	*                                                      *
	*  Copyright 1995 by Brian Mirtich                     *
	*                                                      *
	*  mirtich@cs.berkeley.edu                             *
	*  http://www.cs.berkeley.edu/~mirtich                 *
        *                                                      *
	*******************************************************/

/*
	Revision history

	26 Jan 1996	Program creation.

	 3 Aug 1996	Corrected bug arising when polyhedron density
			is not 1.0.  Changes confined to function main().
			Thanks to Zoran Popovic for catching this one.

	27 May 1997     Corrected sign error in translation of inertia
	                product terms to center of mass frame.  Changes 
			confined to function main().  Thanks to 
			Chris Hecker.
*/

#ifndef VOLINT_HEADER_FILE
#define VOLINT_HEADER_FILE

#include <stdio.h>
#include <iostream>
#include <math.h>
#include <igl/per_face_normals.h>

/*
   ============================================================================
   macros
   ============================================================================
*/

#define SQR(x) ((x)*(x))
#define CUBE(x) ((x)*(x)*(x))


/*
   ============================================================================
   globals
   ============================================================================
*/

static int A;   /* alpha */
static int B;   /* beta */
static int C;   /* gamma */

/* projection integrals */
static double P1, Pa, Pb, Paa, Pab, Pbb, Paaa, Paab, Pabb, Pbbb;

/* face integrals */
static double Fa, Fb, Fc, Faa, Fbb, Fcc, Faaa, Fbbb, Fccc, Faab, Fbbc, Fcca;

/* volume integrals */
static double T0;
static Eigen::RowVector3d T1, T2, TP;


/*
   ============================================================================
   read in a polyhedron
   ============================================================================
*/

/*void readPolyhedron(char *name, POLYHEDRON *p)
{
  FILE *fp;
  char line[200], *c;
  int i, j, n;
  double dx1, dy1, dz1, dx2, dy2, dz2, nx, ny, nz, len;
  FACE *f;

  
  if (!(fp = fopen(name, "r"))) {
    printf("i/o error\n");
    exit(1);
  }
  
  fscanf(fp, "%d", &p->numVerts);
  printf("Reading in %d vertices\n", p->numVerts);
  for (i = 0; i < p->numVerts; i++)
    fscanf(fp, "%lf %lf %lf", 
	   &p->verts[i][X], &p->verts[i][Y], &p->verts[i][Z]);

  fscanf(fp, "%d", &p->numFaces);
  printf("Reading in %d faces\n", p->numFaces);
  for (i = 0; i < p->numFaces; i++) {
    f = &p->faces[i];
    f->poly = p;
    fscanf(fp, "%d", &f->numVerts);
    for (j = 0; j < f->numVerts; j++) fscanf(fp, "%d", &f->verts[j]);

    /* compute face normal and offset w from first 3 vertices */
    /*dx1 = p->verts[f->verts[1]][X] - p->verts[f->verts[0]][X];
    dy1 = p->verts[f->verts[1]][Y] - p->verts[f->verts[0]][Y];
    dz1 = p->verts[f->verts[1]][Z] - p->verts[f->verts[0]][Z];
      
      
      
    dx2 = p->verts[f->verts[2]][X] - p->verts[f->verts[1]][X];
    dy2 = p->verts[f->verts[2]][Y] - p->verts[f->verts[1]][Y];
    dz2 = p->verts[f->verts[2]][Z] - p->verts[f->verts[1]][Z];
    nx = dy1 * dz2 - dy2 * dz1;
    ny = dz1 * dx2 - dz2 * dx1;
    nz = dx1 * dy2 - dx2 * dy1;
    len = sqrt(nx * nx + ny * ny + nz * nz);
    f->norm[X] = nx / len;
    f->norm[Y] = ny / len;
    f->norm[Z] = nz / len;
    f->w = - f->norm[X] * p->verts[f->verts[0]][X]
           - f->norm[Y] * p->verts[f->verts[0]][Y]
           - f->norm[Z] * p->verts[f->verts[0]][Z];

  }

  fclose(fp);

}*/

/*
   ============================================================================
   compute mass properties
   ============================================================================
*/


/* compute various integrations over projection of face */
void compProjectionIntegrals(const Eigen::MatrixXd& V, const Eigen::RowVectorXi& f)
{
  double a0, a1, da;
  double b0, b1, db;
  double a0_2, a0_3, a0_4, b0_2, b0_3, b0_4;
  double a1_2, a1_3, b1_2, b1_3;
  double C1, Ca, Caa, Caaa, Cb, Cbb, Cbbb;
  double Cab, Kab, Caab, Kaab, Cabb, Kabb;
  int i;

  P1 = Pa = Pb = Paa = Pab = Pbb = Paaa = Paab = Pabb = Pbbb = 0.0;

  for (i = 0; i < f.size(); i++) {
      a0 = V(f(i),A);
      b0 = V(f(i),B);
      a1 = V(f((i+1) % f.size()),A);
      b1 = V(f((i+1) % f.size()),B);
      da = a1 - a0;
      db = b1 - b0;
      a0_2 = a0 * a0; a0_3 = a0_2 * a0; a0_4 = a0_3 * a0;
      b0_2 = b0 * b0; b0_3 = b0_2 * b0; b0_4 = b0_3 * b0;
      a1_2 = a1 * a1; a1_3 = a1_2 * a1;
      b1_2 = b1 * b1; b1_3 = b1_2 * b1;

      C1 = a1 + a0;
      Ca = a1*C1 + a0_2; Caa = a1*Ca + a0_3; Caaa = a1*Caa + a0_4;
      Cb = b1*(b1 + b0) + b0_2; Cbb = b1*Cb + b0_3; Cbbb = b1*Cbb + b0_4;
      Cab = 3*a1_2 + 2*a1*a0 + a0_2; Kab = a1_2 + 2*a1*a0 + 3*a0_2;
      Caab = a0*Cab + 4*a1_3; Kaab = a1*Kab + 4*a0_3;
      Cabb = 4*b1_3 + 3*b1_2*b0 + 2*b1*b0_2 + b0_3;
      Kabb = b1_3 + 2*b1_2*b0 + 3*b1*b0_2 + 4*b0_3;

      P1 += db*C1;
      Pa += db*Ca;
      Paa += db*Caa;
      Paaa += db*Caaa;
      Pb += da*Cb;
      Pbb += da*Cbb;
      Pbbb += da*Cbbb;
      Pab += db*(b1*Cab + b0*Kab);
      Paab += db*(b1*Caab + b0*Kaab);
      Pabb += da*(a1*Cabb + a0*Kabb);
  }

  P1 /= 2.0;
  Pa /= 6.0;
  Paa /= 12.0;
  Paaa /= 20.0;
  Pb /= -6.0;
  Pbb /= -12.0;
  Pbbb /= -20.0;
  Pab /= 24.0;
  Paab /= 60.0;
  Pabb /= -60.0;
}

void compFaceIntegrals(const Eigen::MatrixXd& V, const Eigen::RowVectorXi& f, const Eigen::RowVector3d& n)
{
    double w;
    double k1, k2, k3, k4;

    compProjectionIntegrals(V, f);
    w = -n.dot(V.row(f(0)));
 
    k1 = 1 / n(C); k2 = k1 * k1; k3 = k2 * k1; k4 = k3 * k1;

    Fa = k1 * Pa;
    Fb = k1 * Pb;
    Fc = -k2 * (n(A)*Pa + n(B)*Pb + w*P1);

    Faa = k1 * Paa;
    Fbb = k1 * Pbb;
    Fcc = k3 * (SQR(n(A))*Paa + 2*n(A)*n(B)*Pab + SQR(n(B))*Pbb
	 + w*(2*(n(A)*Pa + n(B)*Pb) + w*P1));

    Faaa = k1 * Paaa;
    Fbbb = k1 * Pbbb;
    Fccc = -k4 * (CUBE(n(A))*Paaa + 3*SQR(n(A))*n(B)*Paab
	   + 3*n[A]*SQR(n(B))*Pabb + CUBE(n(B))*Pbbb
	   + 3*w*(SQR(n(A))*Paa + 2*n(A)*n(B)*Pab + SQR(n(B))*Pbb)
	   + w*w*(3*(n(A)*Pa + n(B)*Pb) + w*P1));

    Faab = k1 * Paab;
    Fbbc = -k2 * (n(A)*Pabb + n(B)*Pbbb + w*Pbb);
    Fcca = k3 * (SQR(n(A))*Paaa + 2*n(A)*n(B)*Paab + SQR(n(B))*Pabb
	 + w*(2*(n(A)*Paa + n(B)*Pab) + w*Pa));
}

void compVolumeIntegrals(const Eigen::MatrixXd& V, const Eigen::MatrixXi& T)
{
    Eigen::MatrixXd N;
    igl::per_face_normals(V,T, N);
 
    T0=0.0;
    T1.setZero();
    T2.setZero();
    TP.setZero();

    for (int i= 0; i < T.rows(); i++) {
        Eigen::RowVector3d absn=N.row(i).cwiseAbs();
        if (absn(0) > absn(1) && absn(0) > absn(2)) C = 0;
            else C = (absn(1) > absn(2)) ? 1 : 2;
    
        A = (C + 1) % 3;
        B = (A + 1) % 3;
        compFaceIntegrals(V, T.row(i), N.row(i));
        T0 += N(i,0) * ((A == 0) ? Fa : ((B == 0) ? Fb : Fc));

        T1[A] += N(i,A) * Faa;
        T1[B] += N(i,B) * Fbb;
        T1[C] += N(i,C) * Fcc;
        T2[A] += N(i,A) * Faaa;
        T2[B] += N(i,B) * Fbbb;
        T2[C] += N(i,C) * Fccc;
        TP[A] += N(i,A) * Faab;
        TP[B] += N(i,B) * Fbbc;
        TP[C] += N(i,C) * Fcca;
    }

    T1/= 2;
    T2 /= 3;
    TP/= 2;
}



//computes mass (by uniform density), center of mass, and inertia tensor (around the COM with the canonical axis system) for a polygon

void getCOMandInvIT(const Eigen::MatrixXd& V, const Eigen::MatrixXi& T, const double density, double& mass, Eigen::RowVector3d& COM, Eigen::Matrix3d& invIT){
    
    compVolumeIntegrals(V, T);

    std::cout<<"T0 = "<<T0<<std::endl;

    std::cout<<"Tx ="<<T1[0]<<std::endl;
    std::cout<<"Ty ="<<T1[1]<<std::endl;
    std::cout<<"Tz ="<<T1[2]<<std::endl;
  
    std::cout<<"Txx ="<<T2[0]<<std::endl;
    std::cout<<"Tyy = "<<T2[1]<<std::endl;
    std::cout<<"Tzz ="<<T2[2]<<std::endl;

    std::cout<<"Txy ="<<TP[0]<<std::endl;
    std::cout<<"Tyz ="<<TP[1]<<std::endl;
    std::cout<<"Tzx = "<<TP[2]<<std::endl;

    mass = density * T0;

    /* compute center of mass */
    COM = T1/T0;
    
    Eigen::Matrix3d IT;
 
    /* compute inertia tensor */
    IT(0,0) = density * (T2[1] + T2[2]);
    IT(1,1) = density * (T2[2] + T2[0]);
    IT(2,2) = density * (T2[0] + T2[1]);
    IT(0,1) = IT(1,0) = - density * TP[0];
    IT(1,2) = IT(2,1) = - density * TP[1];
    IT(2,0) = IT(0,2) = - density * TP[2];

    /* translate inertia tensor to center of mass */
    IT(0,0) -= mass * (COM[1]*COM[1] + COM[2]*COM[2]);
    IT(1,1) -= mass * (COM[2]*COM[2] + COM[0]*COM[0]);
    IT(2,2) -= mass * (COM[0]*COM[0] + COM[1]*COM[1]);
    IT(0,1) = IT(1,0) += mass * COM[0] * COM[1];
    IT(1,2) = IT(2,1) += mass * COM[1] * COM[2];
    IT(2,0) = IT(0,2) += mass * COM[2] * COM[0];

    std::cout<<"center of mass: "<<COM<<std::endl;
    std::cout<<"inertia tensor with origin at c.o.m. :"<<IT<<std::endl;
    
    invIT=IT.inverse();  //expensive operation! should only happen once
  
}

#endif
