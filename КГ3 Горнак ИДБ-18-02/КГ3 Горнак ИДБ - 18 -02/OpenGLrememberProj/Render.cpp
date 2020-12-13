

#include "Render.h"

#include <Windows.h>
#include <GL\GL.h>
#include <GL\GLU.h>
#include <math.h>
#include <vector>
#include "Vec3.h"

#define PI 3.14159265

double f(double p1, double p2, double p3, double t)
{
	return p1 * (1 - t)* (1 - t) * (1 - t) + 2 * p2 * t * (1 - t) + p3 * t * t; 
}

double f(double a, double b, double t)
{
	return a * (1 - t) + b * t;
}

double f(double p1, double p2, double p3, double p4, double t)
{
	return p1 * (1 - t) * (1 - t) * (1 - t) + 3 * p2 * t * (1 - t) * (1 - t) + 3 * p3 * t * t * (1 - t) + p4 * t * t * t; 
}

double f2(double p1, double p4, double r1, double r4, double t)
{
	return p1 * (2 * t * t * t - 3 * t * t + 1) + p4 * (-2 * t * t * t + 3 * t * t) + r1 * (t * t * t - 2 * t * t + t) + r4 * (t * t * t - t * t);
}

void bize(double P1[3], double P2[3], double P3[3], double P4[3])
{

	glLineWidth(3); 
	glBegin(GL_LINE_STRIP);
	for (double t = 0; t <= 1.0001; t += 0.01)
	{
		double P[3];
		P[0] = f(P1[0], P2[0], P3[0], P4[0], t);
		P[1] = f(P1[1], P2[1], P3[1], P4[1], t);
		P[2] = f(P1[2], P2[2], P3[2], P4[2], t);
		glVertex3dv(P); 
	}
	glEnd();
	glLineWidth(1); 

	glBegin(GL_LINES);
	glVertex3dv(P1);
	glVertex3dv(P2);
	glVertex3dv(P2);
	glVertex3dv(P3);
	glVertex3dv(P3);
	glVertex3dv(P4);
	glEnd();
}

void ermit(double P11[3], double P4[3], double R1[3], double R4[3])
{
	R1[0] = P11[0] + R1[0];
	R1[1] = P11[1] + R1[1];
	R1[2] = P11[2] + R1[2];

	R4[0] = P4[0] + R4[0];
	R4[1] = P4[1] + R4[1];
	R4[2] = P4[2] + R4[2];

	glBegin(GL_LINES);

	glVertex3dv(P11);
	glVertex3dv(R1);
	glEnd();

	glBegin(GL_LINES);
	glVertex3dv(P4);
	glVertex3dv(R4);
	glEnd();

	glLineWidth(3); 
	glBegin(GL_LINE_STRIP);
	for (double t = 0; t <= 1.0001; t += 0.01)
	{
		double PP[3];

		PP[0] = f2(P11[0], P4[0], R1[0], R4[0], t);
		PP[1] = f2(P11[1], P4[1], R1[1], R4[1], t);
		PP[2] = f2(P11[2], P4[2], R1[2], R4[2], t);

		glVertex3dv(PP);
	}
	glEnd();
	glLineWidth(1);

	
}

Vector3 bizeWithoutDraw(double P1[3], double P2[3], double P3[3], double P4[3], double t)
{
	Vector3 Vec;
	Vec.setCoords(f(P1[0], P2[0], P3[0], P4[0], t), f(P1[1], P2[1], P3[1], P4[1], t), f(P1[2], P2[2], P3[2], P4[2], t));
	return Vec;
}

Vector3 ermitWithoutDraw(double P1[3], double P2[3], double P3[3], double P4[3], double t)
{
	Vector3 Vec;
	Vec.setCoords(f2(P1[0], P2[0], P3[0], P4[0], t), f2(P1[1], P2[1], P3[1], P4[1], t), f2(P1[2], P2[2], P3[2], P4[2], t));
	return Vec;
}

void Draw_Cube(double P1[3])
{
	glBegin(GL_QUADS);


	glColor3d(0.4, 0.4, 0.4);
	glVertex3d(P1[0], P1[1], P1[2]);
	glVertex3d(P1[0] + 1, P1[1], P1[2]);
	glVertex3d(P1[0] + 1, P1[1] + 1, P1[2]);
	glVertex3d(P1[0], P1[1] + 1, P1[2]);

	
	glColor3d(0.1, 0.1, 0.1);
	glVertex3d(P1[0], P1[1], P1[2] + 1);
	glVertex3d(P1[0] + 1, P1[1], P1[2] + 1);
	glVertex3d(P1[0] + 1, P1[1] + 1, P1[2] + 1);
	glVertex3d(P1[0], P1[1] + 1, P1[2] + 1);


	glColor3d(0.2, 0.2, 0.2);
	glVertex3d(P1[0], P1[1], P1[2]);
	glVertex3d(P1[0], P1[1] + 1, P1[2]);
	glVertex3d(P1[0], P1[1] + 1, P1[2] + 1);
	glVertex3d(P1[0], P1[1], P1[2] + 1);

	glColor3d(0.3, 0.3, 0.3);
	glVertex3d(P1[0] + 1, P1[1], P1[2]);
	glVertex3d(P1[0] + 1, P1[1] + 1, P1[2]);
	glVertex3d(P1[0] + 1, P1[1] + 1, P1[2] + 1);
	glVertex3d(P1[0] + 1, P1[1], P1[2] + 1);

	glColor3d(0.4, 0.4, 0.4);
	glVertex3d(P1[0], P1[1], P1[2]);
	glVertex3d(P1[0] + 1, P1[1], P1[2]);
	glVertex3d(P1[0] + 1, P1[1], P1[2] + 1);
	glVertex3d(P1[0], P1[1], P1[2] + 1);

	glColor3d(0.5, 0.5, 0.5);
	glVertex3d(P1[0], P1[1] + 1, P1[2]);
	glVertex3d(P1[0] + 1, P1[1] + 1, P1[2]);
	glVertex3d(P1[0] + 1, P1[1] + 1, P1[2] + 1);
	glVertex3d(P1[0], P1[1] + 1, P1[2] + 1);
	glEnd();
	
	glBegin(GL_LINES);
	glColor3d(1, 0, 0);
	glVertex3d(0, 0, 0);
	glVertex3d(1, 0, 0);

	glColor3d(0, 1, 0);
	glVertex3d(0, 0, 0);
	glVertex3d(0, 1, 0);

	glColor3d(0, 0, 1);
	glVertex3d(0, 0, 0);
	glVertex3d(0, 0, 1);
	glEnd();
}

void LineA(double P1[3], double P2[3], double P3[3], double P4[3], double delta_time)
{
	static double t_max = 0;
	static bool flagReverse = false;
	
	if (!flagReverse)
	{
		t_max += delta_time / 5; 
		if (t_max > 1)
		{
			t_max = 1; 
			flagReverse = !flagReverse;
		}
	}
	else
	{
		t_max -= delta_time / 5; 
		if (t_max < 0)
		{
			t_max = 0; 
			flagReverse = !flagReverse;
		}
	}

	bize(P1, P2, P3, P4);

	
	Vector3 P_old = bizeWithoutDraw(P1, P2, P3, P4, !flagReverse ? t_max - delta_time : t_max + delta_time);
	Vector3 P = bizeWithoutDraw(P1, P2, P3, P4, t_max);
	Vector3 VecP_P_old = (P - P_old).normolize();

	Vector3 rotateX(VecP_P_old.X(), VecP_P_old.Y(), 0);
	rotateX = rotateX.normolize();

	Vector3 VecPrX = Vector3(1, 0, 0).vectProisvedenie(rotateX);
	double CosX = Vector3(1, 0, 0).ScalarProizv(rotateX);
	double SinAngleZ = VecPrX.Z() / abs(VecPrX.Z());
	double AngleOZ = acos(CosX) * 180 / PI * SinAngleZ;

	double AngleOY = acos(VecP_P_old.Z()) * 180 / PI - 90;

	double A[] = { -0.5,-0.5,-0.5 };
	glPushMatrix();
	glTranslated(P.X(), P.Y(), P.Z());
	glRotated(AngleOZ, 0, 0, 1);
	glRotated(AngleOY, 0, 1, 0);
	Draw_Cube(A);
	glPopMatrix();

	glColor3d(0, 0, 0);

}
void LineB(double P1[3], double P2[3], double P3[3], double P4[3], double delta_time)
{
	ermit(P1, P2, P3, P4);

}


void LineC(double P1[3], double P2[3], double P3[3], double P4[3], double delta_time)
{
	ermit(P1, P2, P3, P4);
}

void LineD(double P1[3], double P2[3], double P3[3], double P4[3], double delta_time)
{
	ermit(P1, P2, P3, P4);
}

void Render(double delta_time)
{    
	double P0[] = { 0,0,9 };
	double P1[] = { 4,6,9 };
	double P2[] = { 7,2,0 };
	double P3[] = { 10,-4,0 };

	double PP0[] = { 0,0,3 };
	double PP1[] = { 2,4,4 };
	double PP2[] = { 5,2,3 };
	double PP3[] = { 5,-2,5 };

	double P11[] = { 5,0,4 };
	double P4[] = { -5,0,4 };
	double R1[] = { 0,8,4};
	double R4[] = { 0,-8,4 };

	double P111[] = { -5,0,-1 };
	double P44[] = { 5,0,-1 };
	double R11[] = { 0,-8,-1 };
	double R44[] = { 0,8,-1 };

	LineA(P0, P1, P2, P3, delta_time);
	LineB(PP0, PP1, PP2, PP3, delta_time);
	LineC(P11, P4, R1, R4, delta_time);
	LineD(P111, P44, R11, R44, delta_time);
}   

