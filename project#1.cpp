// Project#1.cpp
// OS : Ubuntu 14.0
// Language: C++


/*
    Marufa Rahmi
    Course: COEN 290
*/

/* This program is the first project of COEN 290

	f - From point
	a - At point

	x - X coordinate
	y - Y coordinate
	z - Z coordinate

	i - increase the current parameter value
	d - decrease the current parameter value

	r - reset the At point to the default position
	c - reset the From point to the default position

	Additional keys:
	1 - to move left arm
	2 - to move right arm
	3 - to move left leg
	4 - to move right leg
	0 - All reset
*/


#include <stdlib.h>
#include <GL/glut.h>
#include <GL/glu.h>
#include <GL/gl.h>
#include <stdio.h>
#include <iostream>

using namespace std;

#define FROM_POINT 1
#define AT_POINT 2


#define X_COORD 1
#define Y_COORD 2
#define Z_COORD 3

struct Point3D{
    GLfloat x,y,z;
};

Point3D From, At, Up;
Point3D camera, object;

int ortho_proj = 1;

// window size and aspect ration
int Width = 600;
int Height = 600;
float aspct_rt = (float)Width / (float)Height;

// initial parameter is from point and coordinate is X
int para_type = FROM_POINT;
int coord = X_COORD;


// increament values for X, Y, znd Z coordinates
float x_inc = 0.2;
float y_inc = 0.2;
float z_inc = 0.2;

// angles for fours limbs
float left_shld_ang;
float right_shld_ang;
float left_leg_ang;
float right_leg_ang;

//function prototype
void init();



void head(){
/*
    Drawing Robot head. Head has three GLUT cube: head frame, eye, mouth piece and
    a GLU sphere for neck.
*/
    glPushMatrix();

        glColor3f(0,1,0);
        glTranslatef(0, 1.2, 0);

        // head frame
        glutWireCube(0.4);

        glPushMatrix();
            // drawing eye frame
            glColor3f(1,1,1);
            glTranslatef(0, 0.1, -0.2);
            glScalef(4, 1, 1);
            glutWireCube(0.1);
        glPopMatrix();

        glPushMatrix();
            // drawing neck
            glColor3f(0, 1, 0.5);
            glTranslatef(0, -0.28, 0);
            glRotatef(90, 1, 0, 0);
            GLUquadricObj* neck;
            neck=gluNewQuadric();
            gluQuadricNormals(neck, GLU_SMOOTH);
            gluQuadricTexture(neck, GL_TRUE);
            gluQuadricDrawStyle(neck, GLU_LINE);
            gluSphere(neck,0.12,10,8);
        glPopMatrix();

         glPushMatrix();
            // drawing mouth piece
            glColor3f(1,1,0);
            glTranslatef(0, -0.1, -0.2);
            glutWireCube(0.18);
        glPopMatrix();

    glPopMatrix();

}



void leftArm(){
/*
    This function draws left arm. arm has five component;
    a GLUT sphere for shoulder
    a GlUT scaled cube for arm
    a GLU cylinder as arm joint
    a GLUT cube  and
    a GLUT Cone
*/
    glPushMatrix();
        glColor3f(0, 1, 0.5);
        glTranslatef(-0.5, 0.65, 0);
        glRotatef(left_shld_ang, 1, 0, 0);
        glutWireSphere(0.1, 10, 10);

        glRotatef(45, 0,1,0);
        glColor3f(0, 1, 0);
        glTranslatef(-0.5,0,0);
        glPushMatrix();
            glScalef(8, 1, 1);
            glutWireCube(0.1);
        glPopMatrix();

        glColor3f(1, 1, 0);
        glTranslatef(-0.4,0,0);
        glRotatef(-90, 0,  1, 0);
        glPushMatrix();
            GLUquadricObj *arm_cyl;
            arm_cyl=gluNewQuadric();
            gluQuadricNormals(arm_cyl, GLU_SMOOTH);
            gluQuadricTexture(arm_cyl, GL_TRUE);
            gluQuadricDrawStyle(arm_cyl, GLU_LINE);
            gluCylinder(arm_cyl,0.05,0.05,0.1,5,5);
        glPopMatrix();

        glTranslatef(0.0,0,0.18);
        glPushMatrix();
            glColor3f(0, 1, 0);
            glutWireCube(0.16);
        glPopMatrix();

        glPushMatrix();
            glTranslatef(0.0, 0,0.18);
            glRotatef(180, 0,1,0);
            glScalef(2, 1,1);
            glutWireCone(0.1,0.1,5,5);
        glPopMatrix();

    glPopMatrix();
}



void rightArm(){
/*
    This function draws right arm. arm has five component;
    a GLUT sphere for shoulder
    a GlUT scaled cube for arm
    a GLU cylinder as arm joint
    a GLUT cube  and
    a GLUT Cone
*/
    glPushMatrix();
        glColor3f(0, 1, 0.5);
        glTranslatef(0.5, 0.65, 0);
        glRotatef(right_shld_ang, 1, 0, 0);
        glutWireSphere(0.1, 10, 10);

        glRotatef(-45, 0,1,0);
        glColor3f(0, 1, 0);
        glTranslatef(0.5,0,0);
        glPushMatrix();
            glScalef(8, 1, 1);
            glutWireCube(0.1);
        glPopMatrix();


        glColor3f(1, 1, 0);
        glTranslatef(0.4,0,0);
        glRotatef(90, 0,  1, 0);
        glPushMatrix();
            GLUquadricObj *arm_cyl;
            arm_cyl=gluNewQuadric();
            gluQuadricNormals(arm_cyl, GLU_SMOOTH);
            gluQuadricTexture(arm_cyl, GL_TRUE);
            gluQuadricDrawStyle(arm_cyl, GLU_LINE);
            gluCylinder(arm_cyl,0.05,0.05,0.1,5,5);
        glPopMatrix();

        glTranslatef(0.0,0,0.18);
        glPushMatrix();
            glColor3f(0, 1, 0);
            glutWireCube(0.16);
        glPopMatrix();

        glPushMatrix();
            glTranslatef(0.0, 0,0.18);
            glRotatef(180, 1,0,0);
            glScalef(2, 1,1);
            glutWireCone(0.1,0.1,5,5);
        glPopMatrix();

    glPopMatrix();
}



void leftLeg(){
/*
    This function draws left Leg. Leg has six component;
    a GLUT sphere to joint to main body
    three GLU cylinder for leg joints
    two GlU cylinder
*/
    glPushMatrix();
        glColor3f(0, 1, 0.5);
        glTranslatef(-0.2, 0.2, 0);
        glRotatef(left_leg_ang, 1, 0, 0);
        glutWireSphere(0.1, 10, 10);

        glColor3f(1, 1, 0);
        glTranslatef(0,0,0.1);
        glPushMatrix();
            GLUquadricObj *leg_joint1;
            leg_joint1=gluNewQuadric();
            gluQuadricNormals(leg_joint1, GLU_SMOOTH);
            gluQuadricTexture(leg_joint1, GL_TRUE);
            gluQuadricDrawStyle(leg_joint1, GLU_LINE);
            gluCylinder(leg_joint1,0.05,0.05,0.1,5,5);
        glPopMatrix();

        glColor3f(0, 1, 0);
        glTranslatef(0,0,0.1);
        glPushMatrix();
            GLUquadricObj *leg_cyl;
            leg_cyl=gluNewQuadric();
            gluQuadricNormals(leg_cyl, GLU_SMOOTH);
            gluQuadricTexture(leg_cyl, GL_TRUE);
            gluQuadricDrawStyle(leg_cyl, GLU_LINE);
            gluCylinder(leg_joint1,0.12,0.08,0.5,10,10);
        glPopMatrix();

        glColor3f(1, 1, 0);
        glTranslatef(0,0,0.5);
        glPushMatrix();
            GLUquadricObj *leg_joint2;
            leg_joint2=gluNewQuadric();
            gluQuadricNormals(leg_joint2, GLU_SMOOTH);
            gluQuadricTexture(leg_joint2, GL_TRUE);
            gluQuadricDrawStyle(leg_joint2, GLU_LINE);
            gluCylinder(leg_joint2,0.05,0.05,0.1,5,5);
        glPopMatrix();

        glColor3f(0, 1, 0);
        glTranslatef(0,0,0.2);
        glPushMatrix();
            glutWireCube(0.2);
        glPopMatrix();

        glColor3f(1, 1, 0);
        glTranslatef(0,0,0.1);
        glPushMatrix();
            GLUquadricObj *leg_joint3;
            leg_joint3=gluNewQuadric();
            gluQuadricNormals(leg_joint3, GLU_SMOOTH);
            gluQuadricTexture(leg_joint3, GL_TRUE);
            gluQuadricDrawStyle(leg_joint3, GLU_LINE);
            gluCylinder(leg_joint3,0.05,0.05,0.2,5,5);
        glPopMatrix();

        glColor3f(0, 1, 0);
        glTranslatef(0,0,0.2);
        glPushMatrix();
            GLUquadricObj *feet_cyl;
            feet_cyl=gluNewQuadric();
            gluQuadricNormals(feet_cyl, GLU_SMOOTH);
            gluQuadricTexture(feet_cyl, GL_TRUE);
            gluQuadricDrawStyle(feet_cyl, GLU_LINE);
            gluCylinder(feet_cyl,0.2,0.2,0.1,6,2);
        glPopMatrix();

    glPopMatrix();
}



void rightLeg(){
/*
    This function draws right Leg. Leg has six component;
    a GLUT sphere to joint to main body
    three GLU cylinder for leg joints
    two GlU cylinder
*/

    glPushMatrix();
        glColor3f(0, 1, 0.5);
        glTranslatef(0.2, 0.2, 0);
        glRotatef(right_leg_ang, 1, 0, 0);
        glutWireSphere(0.1, 10, 10);

        glColor3f(1, 1, 0);
        glTranslatef(0,0,0.1);
        glPushMatrix();
            GLUquadricObj *leg_joint1;
            leg_joint1=gluNewQuadric();
            gluQuadricNormals(leg_joint1, GLU_SMOOTH);
            gluQuadricTexture(leg_joint1, GL_TRUE);
            gluQuadricDrawStyle(leg_joint1, GLU_LINE);
            gluCylinder(leg_joint1,0.05,0.05,0.1,5,5);
        glPopMatrix();

        glColor3f(0, 1, 0);
        glTranslatef(0,0,0.1);
        glPushMatrix();
            GLUquadricObj *leg_cyl;
            leg_cyl=gluNewQuadric();
            gluQuadricNormals(leg_cyl, GLU_SMOOTH);
            gluQuadricTexture(leg_cyl, GL_TRUE);
            gluQuadricDrawStyle(leg_cyl, GLU_LINE);
            gluCylinder(leg_joint1,0.12,0.08,0.5,10,10);
        glPopMatrix();

        glColor3f(1, 1, 0);
        glTranslatef(0,0,0.5);
        glPushMatrix();
            GLUquadricObj *leg_joint2;
            leg_joint2=gluNewQuadric();
            gluQuadricNormals(leg_joint2, GLU_SMOOTH);
            gluQuadricTexture(leg_joint2, GL_TRUE);
            gluQuadricDrawStyle(leg_joint2, GLU_LINE);
            gluCylinder(leg_joint2,0.05,0.05,0.1,5,5);
        glPopMatrix();

        glColor3f(0, 1, 0);
        glTranslatef(0,0,0.2);
        glPushMatrix();
            glutWireCube(0.2);
        glPopMatrix();

        glColor3f(1, 1, 0);
        glTranslatef(0,0,0.1);
        glPushMatrix();
            GLUquadricObj *leg_joint3;
            leg_joint3=gluNewQuadric();
            gluQuadricNormals(leg_joint3, GLU_SMOOTH);
            gluQuadricTexture(leg_joint3, GL_TRUE);
            gluQuadricDrawStyle(leg_joint3, GLU_LINE);
            gluCylinder(leg_joint3,0.05,0.05,0.2,5,5);
         glPopMatrix();

        glColor3f(0, 1, 0);
        glTranslatef(0,0,0.2);
        glPushMatrix();
            GLUquadricObj *feet_cyl;
            feet_cyl=gluNewQuadric();
            gluQuadricNormals(feet_cyl, GLU_SMOOTH);
            gluQuadricTexture(feet_cyl, GL_TRUE);
            gluQuadricDrawStyle(feet_cyl, GLU_LINE);
            gluCylinder(feet_cyl,0.2,0.2,0.1,6,2);
         glPopMatrix();

    glPopMatrix();
}

void body(){
/*
    this function draws body.
    body has two component, a GlUT cone and a GLU disk.
*/
    glPushMatrix();
        glColor3f(0,1,0);
        glTranslatef(0, 0.85 , 0);
        glRotatef(90, 1, 0, 0);

        glutWireCone(0.5, 0.8, 15, 15);

        GLUquadricObj *disk;
        disk=gluNewQuadric();
        gluQuadricNormals(disk, GLU_SMOOTH);
        gluQuadricTexture(disk, GL_TRUE);
        gluQuadricDrawStyle(disk, GLU_LINE);
        gluDisk(disk,0.15, 0.5, 8, 4);

    glPopMatrix();

}



void myReshape(int w, int h)
{
/*
    Update the window if it is reshaped.
*/
    glViewport(0, 0, w, h);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();

	if (ortho_proj) {
		if (w <= h)
			glOrtho(-2.0, 2.0, -2.0 * (1.0/aspct_rt), 2.0 * (1.0/aspct_rt), -20.0, 20.0);
		else
			glOrtho(-2.0 * aspct_rt, 2.0 * aspct_rt, -2.0, 2.0, -20.0, 20.0);
	}
	else
        gluPerspective(50.0, 1.0, 3.0, 30.0);

    glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	gluLookAt(From.x, From.y, From.z, At.x, At.y, At.z, Up.x, Up.y, Up.z);
}



void key(unsigned char key, int x, int y)
{
/*
    Process keyboard inputs
*/

	int w, h;

	switch (key) {
		case 'a':
			para_type = AT_POINT;
			cout<<"At point selected"<<endl;
			break;

		case 'f':
			para_type = FROM_POINT;
			cout<<"From point selected"<<endl;
			break;

		case 'x':
			coord = X_COORD;
			cout<<"X coordinate selected"<<endl;
			break;

		case 'y':
			coord = Y_COORD;
			cout<<"Y coordinate selected"<<endl;
			break;

		case 'z':
			coord = Z_COORD;
			cout<<"Z coordinate selected"<<endl;
			break;

		case 'i':
			if (para_type == AT_POINT) {
				if (coord == X_COORD) At.x += x_inc;
				else if (coord == Y_COORD) At.y += y_inc;
				else At.z += z_inc;
				printf("At point: %.2f %.2f %.2f\n", At.x, At.y, At.z);
			}

			else if (para_type == FROM_POINT) {
				if (coord == X_COORD) From.x += x_inc;
				else if (coord == Y_COORD) From.y += y_inc;
				else From.z += z_inc;
				printf("From point: %.2f %.2f %.2f\n", From.x, From.y, From.z);
			}
			glutPostRedisplay();
			break;

		case 'd':
			if (para_type == AT_POINT) {
				if (coord == X_COORD) At.x -= x_inc;
				else if (coord == Y_COORD) At.y -= y_inc;
				else At.z -= z_inc;
				printf("At point: %.2f %.2f %.2f\n", At.x, At.y, At.z);
			}

			else if (para_type == FROM_POINT) {
				if (coord == X_COORD) From.x -= x_inc;
				else if (coord == Y_COORD) From.y -= y_inc;
				else From.z -= z_inc;
				printf("From point: %.2f %.2f %.2f\n", From.x, From.y, From.z);
			}
			glutPostRedisplay();
			break;

        case 'r':
			At = object;
			cout<<"At point changed to default"<<endl;
			glutPostRedisplay();
			break;

		case 'c':
			From = camera;
			cout<<"From point changed to default"<<endl;
			glutPostRedisplay();
			break;

        case '1':
            left_shld_ang += 2;
			cout<<"Moving Left hand"<<endl;
			glutPostRedisplay();
			break;

        case '2':
            right_shld_ang += 2;
			cout<<"Moving right hand"<<endl;
			glutPostRedisplay();
			break;

        case '3':
            left_leg_ang += 2;
			cout<<"Moving Left leg"<<endl;
			glutPostRedisplay();
			break;

        case '4':
            right_leg_ang += 2;
			cout<<"Moving right leg"<<endl;
			glutPostRedisplay();
			break;

        case '0':
            init();
            cout<<"All Reset"<<endl;
            glutPostRedisplay();
            break;

		case 27:
			exit(0);
			break;
	}
}

void display()
{
/* display callback, clear frame buffer and z buffer,
   rotate object and draw, swap buffers */

    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	//Set up viewing location
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	gluLookAt(From.x, From.y, From.z, At.x, At.y, At.z, Up.x, Up.y, Up.z);

/*
    glColor3f(0, 1.0, 1.0);
    glBegin(GL_LINES);
    for (GLfloat i = -2; i <= 2; i += 0.25) {
        glVertex3f(i, 0, 2);
        glVertex3f(i, 0, -2);
        glVertex3f(2, 0, i);
        glVertex3f(-2, 0, i);

    }
    glEnd();

*/

    head();
    body();
    leftArm();
    rightArm();
    leftLeg();
    rightLeg();

	glutSwapBuffers();
}


void init()
{
/*
    Initialize default parameters.
*/
    GLfloat w, h;

	camera.x = From.x = 1;
	camera.y = From.y = 1;
	camera.z = From.z = -5.0;
	object.x = At.x = 0.0;
	object.y = At.y = 0.0;
	object.z = At.z = 0.0;
	Up.x = 0.0;
	Up.y = 1.0;
	Up.z = 0.0;

	left_shld_ang = 90;
    right_shld_ang = 90;
    left_leg_ang = 90;
    right_leg_ang = 90;


	glMatrixMode(GL_PROJECTION);
    glLoadIdentity();


	if (ortho_proj) {
		if (w <= h)
			glOrtho(-2.0, 2.0, -2.0 * (1.0/aspct_rt), 2.0 * (1.0/aspct_rt), -20.0, 20.0);
		else
			glOrtho(-2.0 * aspct_rt, 2.0 * aspct_rt, -2.0, 2.0, -20.0, 20.0);
	}
	else
		gluPerspective(50.0, 1.0, 3.0, 30.0);

    glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	gluLookAt(From.x, From.y, From.z, At.x, At.y, At.z, Up.x, Up.y, Up.z);
}


main(int argc, char **argv)
{
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
    glutInitWindowSize(Width, Height);
    glutCreateWindow("Robot");
    glutReshapeFunc(myReshape);
    glutDisplayFunc(display);
    glutKeyboardFunc(key);
	init();

    glutMainLoop();
}






























