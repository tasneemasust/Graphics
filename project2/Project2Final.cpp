// Project#2.cpp
// OS : Ubuntu 14.0
// Language: C++


/*
    Marufa Rahmi
    Std ID: W1128039
    Course: COEN 290
*/

/* This program is the 2nd project of COEN 290
    Robot Animation

    features:
    - Hierarchical modeling
    - Shading
    - Mouse picking : user can select different eight parts of the character
    - Rotating selected part using mouse click
    - Video recording
    - playback

    To Play a previously recorded vedio:
    - Run the proram and click play button

    To record a new video:
    - Run the program
    - click record button (the red circle at the bottom)
    - click on the 3D character to select any body part
    - Move the selected part with left or right mouse button
    - click stop button to stop recording and saving the video
    - To playback click play button

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

#define SIZE 512
#define MAXEVENTS 500
#define RECORDSIZE 10

#define body 1
#define head 2
#define neck 21
#define left_arm 3
#define left_lower_arm 31
#define right_arm 4
#define right_lower_arm 41
#define left_leg 5
#define left_lower_leg 51
#define right_leg 6
#define right_lower_leg 61

#define play 1
#define stop 2
#define record 3


struct Point3D{
    GLfloat x,y,z;
};

Point3D From, At, Up;
Point3D camera, object;

int ortho_proj = 1;
bool wire_frame = false;
int obj_id, obj_angle;

// window size for both window
int Width = 600;
int Height = 600;
int lowWinWid = Width;
int lowWinHei = 520;


int cur_selection = body, selection_count = 0;
GLuint cur_selectBuf[512];

// Buffer contains the saved user inputs. Each event is of RECORDSIZE.
// E.g., if you only need to store the object_id and one transformation value (such as angle) per event, then RECORDSIZE is 2.

int event_buffer[MAXEVENTS*RECORDSIZE];

// event_ptr is for recording into the event_buffer array.
int event_ptr = 0;

// playback_ptr is for reading/playing back from the event_buffer array.
int playback_ptr=0;

//recordMode and playbackMode are flags used for recording and playback.
int recordMode = 0;
int playbackMode = 0;

//Recorded events are saved to <first nameLastInitial>
FILE *jFile = NULL;
char *fileName = "MarufaR.txt";

// angles for all body parts
float theta[100];

GLfloat mat_specular1[] ={0.628281, 0.555802, 0.366065, 1.0};
GLfloat mat_ambient1[] ={0.24725, 0.1995, 0.0745, 1.0};
GLfloat mat_diffuse1[] ={0.75164, 0.60648, 0.22648, 1.0};
GLfloat mat_shininess1[] ={128.0 * 0.4};

GLfloat mat_specular2[] ={0.508273, 0.508273, 0.508373};
GLfloat mat_ambient2[] ={0.19225, 0.19225, 0.19225};
GLfloat mat_diffuse2[] ={0.50754, 0.50754, 0.50754};
GLfloat mat_shininess2[] ={128.0 * 0.6};

GLfloat mat_specular3[] ={0.296648, 0.296648, 0.296648};
GLfloat mat_ambient3[] ={0.25, 0.20725, 0.20725};
GLfloat mat_diffuse3[] ={1, 0.829, 0.829};
GLfloat mat_shininess3[] ={128.0 * 0.088};

GLfloat mat_specular4[] ={0.633, 0.727811, 0.633};
GLfloat mat_ambient4[] ={0.0215, 0.1745, 0.0215};
GLfloat mat_diffuse4[] ={0.07568, 0.61424, 0.07568};
GLfloat mat_shininess4[] ={128 * 0.8};

GLfloat mat_specular5[] = {0.8F, 0.0F, 0.0F, 1.0F};
GLfloat mat_diffuse5[] = {0.6F, 0.0F, 0.0F, 1.0F};
GLfloat mat_ambient5[] = {0.8F, 0.5F, 0.5F, 1.0F};
GLfloat mat_shininess5[]= {100.0};



//function prototype
void init();
void display();
void drawObjects(GLenum mode);
void reset_angles();
void ortho(int w, int h);
void drawObjects2(GLenum mode);
void userEventAction(int key);
void setMenuEntries(bool init);
void SetMaterial(GLfloat spec[], GLfloat amb[], GLfloat diff[], GLfloat shin[]);


/**************************************************************************************************************/
// next six function draws the 3D character

void drawHead(GLenum mode){
/*
    Drawing Robot head. Head has three GLUT cube: head frame, eye, mouth piece and
    a GLU sphere for neck.
*/
    if(mode == GL_SELECT)glPushName(head);
    glPushMatrix();

        //glColor3f(0.4,0.3,0.3);
        glColor3f(0,1,0);
        glTranslatef(0, 1.2, 0);
        glRotatef(theta[head], 0, 1, 0);

        // head frame
        SetMaterial(mat_specular2, mat_diffuse2, mat_ambient2, mat_shininess4);
        if(wire_frame)glutWireCube(0.4);
        else glutSolidCube(0.4);

        glPushMatrix();
            // drawing eye frame
            SetMaterial(mat_specular1, mat_diffuse1, mat_ambient1, mat_shininess1);
            glColor3f(1,1,1);
            glTranslatef(0, 0.1, 0.2);
            glScalef(3, 1, 1);
            if(wire_frame)glutWireCube(0.1);
            else glutSolidCube(0.1);
        glPopMatrix();

        glPushMatrix();
            // drawing neck
            SetMaterial(mat_specular1, mat_diffuse1, mat_ambient1, mat_shininess1);
            glPushName(neck);
            glColor3f(0, 1, 0.5);
            glTranslatef(0, -0.28, 0);
            //glRotatef(90, 1, 0, 0);
            GLUquadricObj* neckObj;
            neckObj=gluNewQuadric();
            gluQuadricNormals(neckObj, GLU_SMOOTH);
            gluQuadricTexture(neckObj, GL_TRUE);
            if(wire_frame)gluQuadricDrawStyle(neckObj, GLU_LINE);
            else gluQuadricDrawStyle(neckObj, GLU_FILL);
            gluSphere(neckObj,0.12,10,8);
            glPopName();
        glPopMatrix();

         glPushMatrix();
            // drawing mouth piece
            SetMaterial(mat_specular3, mat_diffuse3, mat_ambient3, mat_shininess3);
            glColor3f(1,1,0);
            glTranslatef(0, -0.1, 0.2);
            if(wire_frame)glutWireCube(0.18);
            else glutSolidCube(0.18);
        glPopMatrix();

    glPopMatrix();
    if(mode == GL_SELECT)glPopName();

}


void drawBody(GLenum mode){
/*
    this function draws body.
    body has two component, a GlUT cone and a GLU disk.
*/
    if(mode == GL_SELECT)glPushName(body);
    SetMaterial(mat_specular2, mat_diffuse2, mat_ambient2, mat_shininess4);
    glPushMatrix();

        glColor3f(0,1,0);
        glTranslatef(0, 0.85 , 0);
        glRotatef(90, 1, 0, 0);

        if(wire_frame)glutWireCone(0.5, 0.8, 15, 15);
        else glutSolidCone(0.5, 0.8, 15, 15);

        glPushMatrix();
        GLUquadricObj *disk;
        disk=gluNewQuadric();
        gluQuadricNormals(disk, GLU_SMOOTH);
        gluQuadricTexture(disk, GL_TRUE);
        if(wire_frame)gluQuadricDrawStyle(disk, GLU_LINE);
        else gluQuadricDrawStyle(disk, GLU_FILL);
        gluDisk(disk,0.15, 0.5, 8, 4);
        glPopMatrix();
    glPopMatrix();
   if(mode == GL_SELECT) glPopName();
}


void leftArm(GLenum mode){
/*
    This function draws left arm. arm has five component;
    a GLUT sphere for shoulder
    a GlUT scaled cube for arm
    a GLU cylinder as arm joint
    a GLUT cube  and
    a GLUT Cone
*/
    if(mode == GL_SELECT)glPushName(left_arm);
    glPushMatrix();
        SetMaterial(mat_specular1, mat_diffuse1, mat_ambient1, mat_shininess1);
        glColor3f(0, 1, 0.5);
        glTranslatef(-0.5, 0.65, 0);
        glRotatef(theta[left_arm], 1, 0, 0);
        if(wire_frame)glutWireSphere(0.1, 10, 10);
        else glutSolidSphere(0.1, 10, 10);

        glRotatef(60, 0,1,0);
        glColor3f(0, 1, 0);
        glTranslatef(-0.5,0,0);
        glPushMatrix();
            SetMaterial(mat_specular2, mat_diffuse2, mat_ambient2, mat_shininess4);
            glScalef(8, 1, 1);
            if(wire_frame)glutWireCube(0.1);
            else glutSolidCube(0.1);
        glPopMatrix();
        if(mode == GL_SELECT)glPopName();

        if(mode == GL_SELECT)glPushName(left_lower_arm);
        glColor3f(1, 1, 0);
        glTranslatef(-0.4,0,0);
        glRotatef(-90, 0,  1, 0);
        glRotatef(theta[left_lower_arm], 1, 0 , 0);
        glPushMatrix();
            SetMaterial(mat_specular1, mat_diffuse1, mat_ambient1, mat_shininess1);
            GLUquadricObj *arm_cyl;
            arm_cyl=gluNewQuadric();
            gluQuadricNormals(arm_cyl, GLU_SMOOTH);
            gluQuadricTexture(arm_cyl, GL_TRUE);
            if(wire_frame)gluQuadricDrawStyle(arm_cyl, GLU_LINE);
            else gluQuadricDrawStyle(arm_cyl, GLU_FILL);
            gluCylinder(arm_cyl,0.05,0.05,0.1,5,5);
        glPopMatrix();

        glTranslatef(0.0,0,0.18);
        glPushMatrix();
            SetMaterial(mat_specular2, mat_diffuse2, mat_ambient2, mat_shininess4);
            glColor3f(0, 1, 0);
            if(wire_frame)glutWireCube(0.16);
            else glutSolidCube(0.16);
        glPopMatrix();

        glPushMatrix();
            SetMaterial(mat_specular1, mat_diffuse1, mat_ambient1, mat_shininess1);
            glTranslatef(0.0, 0,0.18);
            glRotatef(180, 0,1,0);
            glScalef(2, 1,1);
            glutWireCone(0.1,0.1,5,5);
        glPopMatrix();
        if(mode == GL_SELECT)glPopName();

    glPopMatrix();
}



void rightArm(GLenum mode){
/*
    This function draws right arm. arm has five component;
    a GLUT sphere for shoulder
    a GlUT scaled cube for arm
    a GLU cylinder as arm joint
    a GLUT cube  and
    a GLUT Cone
*/
    if(mode == GL_SELECT)glPushName(right_arm);
    glPushMatrix();
        SetMaterial(mat_specular1, mat_diffuse1, mat_ambient1, mat_shininess1);
        glColor3f(0, 1, 0.5);
        glTranslatef(0.5, 0.65, 0);
        glRotatef(theta[right_arm], 1, 0, 0);
        if(wire_frame)glutWireSphere(0.1, 10, 10);
        else glutSolidSphere(0.1, 10, 10);

        glRotatef(-60, 0,1,0);
        glColor3f(0, 1, 0);
        glTranslatef(0.5,0,0);
        glPushMatrix();
            SetMaterial(mat_specular2, mat_diffuse2, mat_ambient2, mat_shininess4);
            glScalef(8, 1, 1);
            if(wire_frame)glutWireCube(0.1);
            else glutSolidCube(0.1);
        glPopMatrix();
        if(mode == GL_SELECT)glPopName();

        glPushName(right_lower_arm);
        glColor3f(1, 1, 0);
        glTranslatef(0.4,0,0);
        glRotatef(90, 0, 1, 0);
        glRotatef(theta[right_lower_arm], 1,  0, 0);
        glPushMatrix();
            SetMaterial(mat_specular1, mat_diffuse1, mat_ambient1, mat_shininess1);
            GLUquadricObj *arm_cyl;
            arm_cyl=gluNewQuadric();
            gluQuadricNormals(arm_cyl, GLU_SMOOTH);
            gluQuadricTexture(arm_cyl, GL_TRUE);
            if(wire_frame)gluQuadricDrawStyle(arm_cyl, GLU_LINE);
            else gluQuadricDrawStyle(arm_cyl, GLU_FILL);
            gluCylinder(arm_cyl,0.05,0.05,0.1,5,5);
        glPopMatrix();

        glTranslatef(0.0,0,0.18);
        glPushMatrix();
            SetMaterial(mat_specular2, mat_diffuse2, mat_ambient2, mat_shininess4);
            glColor3f(0, 1, 0);
            if(wire_frame)glutWireCube(0.16);
            else glutSolidCube(0.16);
        glPopMatrix();

        glPushMatrix();
            SetMaterial(mat_specular1, mat_diffuse1, mat_ambient1, mat_shininess1);
            glTranslatef(0.0, 0,0.18);
            glRotatef(180, 1,0,0);
            glScalef(2, 1,1);
            glutWireCone(0.1,0.1,5,5);
        glPopMatrix();
        if(mode == GL_SELECT)glPopName();

    glPopMatrix();
}



void leftLeg(GLenum mode){
/*
    This function draws left Leg. Leg has six component;
    a GLUT sphere to joint to main body
    three GLU cylinder for leg joints
    two GlU cylinder
*/
    if(mode == GL_SELECT)glPushName(left_leg);
    SetMaterial(mat_specular1, mat_diffuse1, mat_ambient1, mat_shininess1);
    glPushMatrix();
        glColor3f(0, 1, 0.5);
        glTranslatef(-0.2, 0.2, 0);
        glRotatef(theta[left_leg], 1, 0, 0);
        if(wire_frame)glutWireSphere(0.1, 10, 10);
        else glutSolidSphere(0.1, 10, 10);

        glColor3f(1, 1, 0);
        glTranslatef(0,0,0.1);
        glPushMatrix();
            SetMaterial(mat_specular1, mat_diffuse1, mat_ambient1, mat_shininess1);
            GLUquadricObj *leg_joint1;
            leg_joint1=gluNewQuadric();
            gluQuadricNormals(leg_joint1, GLU_SMOOTH);
            gluQuadricTexture(leg_joint1, GL_TRUE);
            if(wire_frame)gluQuadricDrawStyle(leg_joint1, GLU_LINE);
            else gluQuadricDrawStyle(leg_joint1, GLU_FILL);
            gluCylinder(leg_joint1,0.05,0.05,0.1,5,5);
        glPopMatrix();

        glColor3f(0, 1, 0);
        glTranslatef(0,0,0.1);
        glPushMatrix();
            SetMaterial(mat_specular2, mat_diffuse2, mat_ambient2, mat_shininess2);
            GLUquadricObj *leg_cyl;
            leg_cyl=gluNewQuadric();
            gluQuadricNormals(leg_cyl, GLU_SMOOTH);
            gluQuadricTexture(leg_cyl, GL_TRUE);
            if(wire_frame)gluQuadricDrawStyle(leg_cyl, GLU_LINE);
            else gluQuadricDrawStyle(leg_cyl, GLU_FILL);
            gluCylinder(leg_joint1,0.12,0.08,0.5,10,10);
        glPopMatrix();
        if(mode == GL_SELECT)glPopName();

        glPushName(left_lower_leg);
        glColor3f(1, 1, 0);
        glTranslatef(0,0,0.5);
        glRotatef(theta[left_lower_leg], 1, 0, 0);
        glPushMatrix();
            SetMaterial(mat_specular1, mat_diffuse1, mat_ambient1, mat_shininess1);
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
            SetMaterial(mat_specular2, mat_diffuse2, mat_ambient2, mat_shininess4);
            if(wire_frame)glutWireCube(0.2);
            else glutSolidCube(0.2);
        glPopMatrix();

        glColor3f(1, 1, 0);
        glTranslatef(0,0,0.1);
        glPushMatrix();
            SetMaterial(mat_specular1, mat_diffuse1, mat_ambient1, mat_shininess1);
            GLUquadricObj *leg_joint3;
            leg_joint3=gluNewQuadric();
            gluQuadricNormals(leg_joint3, GLU_SMOOTH);
            gluQuadricTexture(leg_joint3, GL_TRUE);
            if(wire_frame)gluQuadricDrawStyle(leg_joint3, GLU_LINE);
            else gluQuadricDrawStyle(leg_joint3, GLU_FILL);
            gluCylinder(leg_joint3,0.05,0.05,0.2,5,5);
        glPopMatrix();

        glColor3f(0, 1, 0);
        glTranslatef(0,0,0.2);
        glPushMatrix();
            SetMaterial(mat_specular2, mat_diffuse2, mat_ambient2, mat_shininess4);
            GLUquadricObj *feet_cyl;
            feet_cyl=gluNewQuadric();
            gluQuadricNormals(feet_cyl, GLU_SMOOTH);
            gluQuadricTexture(feet_cyl, GL_TRUE);
            if(wire_frame)gluQuadricDrawStyle(feet_cyl, GLU_LINE);
            else gluQuadricDrawStyle(feet_cyl, GLU_FILL);
            gluCylinder(feet_cyl,0.2,0.2,0.1,6,2);
        glPopMatrix();

        glPushMatrix();
            SetMaterial(mat_specular2, mat_diffuse2, mat_ambient2, mat_shininess4);
            GLUquadricObj *disk;
            disk=gluNewQuadric();
            gluQuadricNormals(disk, GLU_SMOOTH);
            gluQuadricTexture(disk, GL_TRUE);
            if(wire_frame)gluQuadricDrawStyle(disk, GLU_LINE);
            else gluQuadricDrawStyle(disk, GLU_FILL);
            gluDisk(disk,0.05, 0.2, 6, 4);
        glPopMatrix();
        if(mode == GL_SELECT)glPopName();

    glPopMatrix();
}



void rightLeg(GLenum mode){
/*
    This function draws right Leg. Leg has six component;
    a GLUT sphere to joint to main body
    three GLU cylinder for leg joints
    two GlU cylinder
*/
    if(mode == GL_SELECT)glPushName(right_leg);
    glPushMatrix();
        SetMaterial(mat_specular1, mat_diffuse1, mat_ambient1, mat_shininess1);
        glColor3f(0, 1, 0.5);
        glTranslatef(0.2, 0.2, 0);
        glRotatef(theta[right_leg], 1, 0, 0);
        if(wire_frame)glutWireSphere(0.1, 10, 10);
        else glutSolidSphere(0.1, 10, 10);

        glColor3f(1, 1, 0);
        glTranslatef(0,0,0.1);
        glPushMatrix();
            SetMaterial(mat_specular1, mat_diffuse1, mat_ambient1, mat_shininess1);
            GLUquadricObj *leg_joint1;
            leg_joint1=gluNewQuadric();
            gluQuadricNormals(leg_joint1, GLU_SMOOTH);
            gluQuadricTexture(leg_joint1, GL_TRUE);
            if(wire_frame)gluQuadricDrawStyle(leg_joint1, GLU_LINE);
            else gluQuadricDrawStyle(leg_joint1, GLU_FILL);
            gluCylinder(leg_joint1,0.05,0.05,0.1,5,5);
        glPopMatrix();

        glColor3f(0, 1, 0);
        glTranslatef(0,0,0.1);
        glPushMatrix();
            SetMaterial(mat_specular2, mat_diffuse2, mat_ambient2, mat_shininess4);
            GLUquadricObj *leg_cyl;
            leg_cyl=gluNewQuadric();
            gluQuadricNormals(leg_cyl, GLU_SMOOTH);
            gluQuadricTexture(leg_cyl, GL_TRUE);
            if(wire_frame)gluQuadricDrawStyle(leg_cyl, GLU_LINE);
            else gluQuadricDrawStyle(leg_cyl, GLU_FILL);
            gluCylinder(leg_joint1,0.12,0.08,0.5,10,10);
        glPopMatrix();
        if(mode == GL_SELECT)glPopName();

        glPushName(right_lower_leg);
        glColor3f(1, 1, 0);
        glTranslatef(0,0,0.5);
        glRotatef(theta[right_lower_leg], 1, 0,0);
        glPushMatrix();
            SetMaterial(mat_specular1, mat_diffuse1, mat_ambient1, mat_shininess1);
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
            SetMaterial(mat_specular2, mat_diffuse2, mat_ambient2, mat_shininess4);
            if(wire_frame)glutWireCube(0.2);
            else glutSolidCube(0.2);
        glPopMatrix();

        glColor3f(1, 1, 0);
        glTranslatef(0,0,0.1);
        glPushMatrix();
            SetMaterial(mat_specular1, mat_diffuse1, mat_ambient1, mat_shininess1);
            GLUquadricObj *leg_joint3;
            leg_joint3=gluNewQuadric();
            gluQuadricNormals(leg_joint3, GLU_SMOOTH);
            gluQuadricTexture(leg_joint3, GL_TRUE);
            if(wire_frame)gluQuadricDrawStyle(leg_joint3, GLU_LINE);
            else gluQuadricDrawStyle(leg_joint3, GLU_FILL);
            gluCylinder(leg_joint3,0.05,0.05,0.2,5,5);
         glPopMatrix();

        glColor3f(0, 1, 0);
        glTranslatef(0,0,0.2);
        glPushMatrix();
            SetMaterial(mat_specular2, mat_diffuse2, mat_ambient2, mat_shininess4);
            GLUquadricObj *feet_cyl;
            feet_cyl=gluNewQuadric();
            gluQuadricNormals(feet_cyl, GLU_SMOOTH);
            gluQuadricTexture(feet_cyl, GL_TRUE);
            if(wire_frame)gluQuadricDrawStyle(feet_cyl, GLU_LINE);
            else gluQuadricDrawStyle(feet_cyl, GLU_FILL);
            gluCylinder(feet_cyl,0.2,0.2,0.1,6,2);
         glPopMatrix();

         glPushMatrix();
            SetMaterial(mat_specular2, mat_diffuse2, mat_ambient2, mat_shininess4);
            GLUquadricObj *disk;
            disk=gluNewQuadric();
            gluQuadricNormals(disk, GLU_SMOOTH);
            gluQuadricTexture(disk, GL_TRUE);
            if(wire_frame)gluQuadricDrawStyle(disk, GLU_LINE);
            else gluQuadricDrawStyle(disk, GLU_FILL);
            gluDisk(disk,0.05, 0.2, 6, 4);
        glPopMatrix();
         if(mode == GL_SELECT)glPopName();

    glPopMatrix();
}

/* This function is unused*/
void drawBall(){
    glPushMatrix();
    SetMaterial(mat_specular3, mat_diffuse3, mat_ambient3, mat_shininess1);
    glTranslatef(-0.6,-1,0);
    glRotatef(45,1,0,0);
    glutSolidSphere(0.15, 20, 20);
    glPopMatrix();
}

/*********************************************************************************************************************/
// Next four functions manage the mouse picking and storing state of the object

// This fuction save all angles of all body parts
void saveState(){
    event_buffer[event_ptr++] = head;
    event_buffer[event_ptr++] = theta[head];
    event_buffer[event_ptr++] = body;
    event_buffer[event_ptr++] = theta[body];
    event_buffer[event_ptr++] = left_arm;
    event_buffer[event_ptr++] = theta[left_arm];
    event_buffer[event_ptr++] = left_lower_arm;
    event_buffer[event_ptr++] = theta[left_lower_arm];
    event_buffer[event_ptr++] = right_leg;
    event_buffer[event_ptr++] = theta[right_leg];
    event_buffer[event_ptr++] = right_lower_leg;
    event_buffer[event_ptr++] = theta[right_lower_leg];
    event_buffer[event_ptr++] = right_arm;
    event_buffer[event_ptr++] = theta[right_arm];
    event_buffer[event_ptr++] = right_lower_arm;
    event_buffer[event_ptr++] = theta[right_lower_arm];
    event_buffer[event_ptr++] = left_leg;
    event_buffer[event_ptr++] = theta[left_leg];
    event_buffer[event_ptr++] = left_lower_leg;
    event_buffer[event_ptr++] = theta[left_lower_leg];
}

//this fuction rotate the selected object
void processHitsMove (GLint hits,  int click)
{
    float angle;
    unsigned int i, j;
    int obj,n;
    n = selection_count;

    for (j = 0; j < n; j++) { /*  for each object */
        obj = cur_selectBuf[j];
        //cout<<"obj "<< obj<<" n "<<n<<endl;

        if(click == 0 ){
            theta[obj] +=10.0;
            if( theta[obj] > 360.0 ) theta[obj] -= 360.0;
        }
        else{
            theta[obj] -= 10.0;
            if( theta[obj] < 360.0 ) theta[obj] += 360.0;
        }
        if(recordMode == 1){
            saveState(); //saving the state for every movement
        }
    }
}


void processHits (GLint hits, GLuint buffer[])
{
    unsigned int i, j;
    GLuint ii, jj, names, *ptr;

    ptr = (GLuint *) buffer;
    if(hits > 0)selection_count = 0;
    for (i = 0; i < hits; i++) {	/*  for each hit  */
        names = *ptr;
        ptr+=3;
        for (j = 0; j < names; j++) { /*  for each name */

            if(*ptr == head) cout<<"Head selected"<<endl;
            if(*ptr == neck) cout<<"Neck selected"<<endl;
            if(*ptr == body) cout<<"Body selected"<<endl;
            if(*ptr == left_arm) cout<<"Left arm selected"<<endl;
            if(*ptr == left_lower_arm) cout<<"Left lower arm selected"<<endl;
            if(*ptr == right_arm) cout<<"Right arm selected"<<endl;
            if(*ptr == right_lower_arm) cout<<"Right lower arm selected"<<endl;
            if(*ptr == left_leg) cout<<"Left leg selected"<<endl;
            if(*ptr == left_lower_leg) cout<<"Left lower leg selected"<<endl;
            if(*ptr == right_leg) cout<<"Right leg selected"<<endl;
            if(*ptr == right_lower_leg) cout<<"Right lower leg selected"<<endl;
            cur_selectBuf[selection_count++] = *ptr;
            //cout<<*ptr<<endl;
            ptr++;
        }
        printf ("\n");
    }
}


void mouse(int button, int state, int x, int y)
{
    GLuint selectBuf[SIZE];
    GLint hits;
    GLint viewport[4];

    if (button == GLUT_LEFT_BUTTON && state == GLUT_DOWN){
        // set picking window as upperwindow
        glutSetWindow(2);
        glGetIntegerv (GL_VIEWPORT, viewport);

        glSelectBuffer (SIZE, selectBuf);
        glRenderMode(GL_SELECT);

        glInitNames();
        glPushName(0);

        glPushMatrix ();
        glLoadIdentity ();
    // create 5x5 pixel picking region near cursor location
        gluPickMatrix ((GLdouble) x, (GLdouble) (viewport[3] - y), 5.0, 5.0, viewport);
        ortho(Width, Height);
        drawObjects(GL_SELECT);


        glMatrixMode (GL_PROJECTION);
        glPopMatrix ();
        glFlush ();

        hits = glRenderMode (GL_RENDER);
        processHits (hits, selectBuf);
        if(hits == 0)processHitsMove(hits,  0);

        glutPostRedisplay();

    }
    if(button == GLUT_RIGHT_BUTTON && state == GLUT_DOWN){
        processHitsMove(hits, 1);
        glutPostRedisplay();
    }
}


/****************************************************************************************************************/
// next seven functions manage video operations e.g. record, save, play etc

void timerFunc(int val)
{
    // Check if playback_ptr has reached the last event
    if(playback_ptr<event_ptr)
    {
        int event_counter = 10;

        // redisplay waits until 10 objects are fetched from file
        for(int i = 0; i<event_counter; i++){
            obj_id = event_buffer[playback_ptr++];
            obj_angle = event_buffer[playback_ptr++];

		// Update the object's transformation value

            switch(obj_id)
            {
                case body:
                    theta[body] = obj_angle;
                    break;
                case head:
                    theta[head] = obj_angle;
                    break;
                case left_arm:
                    theta[left_arm] = obj_angle;
                    break;
                case left_lower_arm:
                    theta[left_lower_arm] = obj_angle;
                    break;
                case right_arm:
                    theta[right_arm] = obj_angle;
                    break;
                case right_lower_arm:
                    theta[right_lower_arm] = obj_angle;
                    break;
                case left_leg:
                    theta[left_leg] = obj_angle;
                    break;
                case left_lower_leg:
                    theta[left_lower_leg] = obj_angle;
                    break;
                case right_leg:
                    theta[right_leg] = obj_angle;
                    break;
                case right_lower_leg:
                    theta[right_lower_leg] = obj_angle;
                    break;
            }
        }
		// Update the screen with the current transformation retrieved from the event_buffer.
		glutSetWindow(2);
		display();

        // Call TimerFunc to read another event
		glutTimerFunc(200, timerFunc, 1);
    }
	else{
		playback_ptr=0;
    }
}


void recordBegin(){
    ortho_proj = 1;
    recordMode = 1;
    playbackMode = 0;
    printf("Recording...  Press Sotp button or 'e' to end\n");
    event_ptr=0;
}


void recordEnd(){
    if(recordMode == 1){
        recordMode = 0;
        printf("Recording stopped.\n");
    }
}

//load event_buffer from file
void loadFile(){
    recordMode = 0;
    playbackMode = 0;

    event_ptr=0;
    playback_ptr=0;
    reset_angles();
    printf("Loading file %s\n", fileName);

    jFile = fopen(fileName, "r");
    if ( jFile == NULL ){
        printf("Warning: Could not open %s\n", fileName);
        playbackMode = 0;
    }
    else {
					// Store the events to event_buffer
        while((fscanf(jFile, "%d ", &event_buffer[event_ptr])) != EOF){
            event_ptr++;
        }
        fclose(jFile);
        playbackMode = 1;
    }
}


void playButton(){
    if(playbackMode==1){
        glutTimerFunc(4,timerFunc,1);
    }
}


//saving event buffer to file
void saveFile(){
    recordMode = 0;
    playbackMode = 0;

	jFile = fopen(fileName, "w");
	if (jFile == NULL){
		printf("Warning: Could not open %s\n", fileName);
	}
	else {
		for(int j=0;j<event_ptr;j++){
			fprintf(jFile, "%d ", event_buffer[j]);
		}
        fclose(jFile);
		printf("\nEvents saved in %s\n", fileName);
    }
	playback_ptr=0;
}


void key(unsigned char key, int x, int y){

   // Process keyboard inputs
   //video recording has two options by keys or buttons on the screen


	switch(key)
	{
		case 'b': //begin recording. Set recordMode
			recordBegin();
			break;

		case 'e': //stop recording. Reset recordMode.
			recordEnd();
			//Save file.
            saveFile();
			break;

		case 'r': //Reset everything.
			recordMode = 0;
			playbackMode = 0;

			event_ptr=0;
			playback_ptr=0;
			reset_angles();
			break;

		case 'p': //Playback
                //Load file. Reset everything and load the contents of the file into the buffer.
            loadFile();
            playButton();
			break;

        case '0':
            init();
            cout<<"All Reset"<<endl;
            glutPostRedisplay();
            break;

        //some saved movement of the 3D character with keys x,z,a,c
        case 'x':
            reset_angles();
            theta[body] += 5;
            theta[left_arm] += 45;
            theta[left_lower_arm] -=90;
            theta[right_leg] += 45;
            theta[right_lower_leg] += 45;

            theta[right_arm] -= 45;
            theta[right_lower_arm] -=90;
            theta[left_leg] -= 45;
            theta[left_lower_leg] += 45;
            if(recordMode == 1){
                saveState();
            }
            glutPostRedisplay();
            break;

        case 'z':
            reset_angles();
            theta[body] -= 5;
            theta[left_arm] -= 45;
            theta[left_lower_arm] -=90;
            theta[right_leg] -= 45;
            theta[right_lower_leg] += 45;

            theta[right_arm] += 45;
            theta[right_lower_arm] -=90;
            theta[left_leg] += 45;
            theta[left_lower_leg] += 45;
            if(recordMode == 1){
                saveState();
            }
            glutPostRedisplay();
            break;
        case 'a':
            theta[left_arm] -= 90;

            if(recordMode == 1){
                saveState();
            }
            glutPostRedisplay();
            break;
        case 'c':
            theta[right_arm] -= 90;
            if(recordMode == 1){
                saveState();
            }
            glutPostRedisplay();
            break;

		case 27:
			exit(0);
			break;
	}
}


/**************************************************************************************************************/
// Next Five functions are callback functions of the main window or display window


void init(){

    //Initialize default parameters of upper window/ view window.

	camera.x = From.x = 0.5;
	camera.y = From.y = 1;
	camera.z = From.z = 5.0;
	object.x = At.x = 0.0;
	object.y = At.y = 0.0;
	object.z = At.z = 0.0;
	Up.x = 0.0;
	Up.y = 1.0;
	Up.z = 0.0;

    reset_angles();
    if(!ortho_proj)setMenuEntries(true);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();

    ortho(Width, Height);

    glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	gluLookAt(From.x, From.y, From.z, At.x, At.y, At.z, Up.x, Up.y, Up.z);
}


// Floor of viewing window
void Enviro(char solid){

    int i, j;
    float x = 4;
    glColor3f(1.0, 1.0, 0.0);
    glTranslatef(0, -1.4, 0);
    SetMaterial(mat_specular4, mat_diffuse4, mat_ambient4, mat_shininess4);
    glBegin(GL_POLYGON);
    glVertex3f(-4, -1, -1);
    glVertex3f(-1, .3, -4);
    glVertex3f(1, .3, -4 );
    glVertex3f(3, -1, -1);
  glEnd();

}

void drawObjects(GLenum mode){

    glPushMatrix();
    glRotatef(theta[body], 0, 1, 0);
    drawBody(mode);
	drawHead(mode);
	leftArm(mode);
	rightArm(mode);
	leftLeg(mode);
	rightLeg(mode);

	glPopMatrix();
	//drawBall(); // addition object
}


void display()
{
    //display callback, clear frame buffer and z buffer,rotate object and draw, swap buffers

    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	//Set up viewing location
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	gluLookAt(From.x, From.y, From.z, At.x, At.y, At.z, Up.x, Up.y, Up.z);

    drawObjects(GL_RENDER);
    if(!ortho_proj);
    Enviro(1);
   // glutSwapBuffers();
	glFlush();
}


void myReshape(int w, int h){

    //Update the window if it is reshaped.

    glViewport(0, 0, w, h);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();

    ortho(Width, Height);

    glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	gluLookAt(From.x, From.y, From.z, At.x, At.y, At.z, Up.x, Up.y, Up.z);
}

void SetMaterial(GLfloat spec[], GLfloat amb[], GLfloat diff[], GLfloat shin[]){
     glShadeModel(GL_SMOOTH);
     glMaterialfv(GL_FRONT, GL_SPECULAR, spec);
     glMaterialfv(GL_FRONT, GL_SHININESS, shin);
     glMaterialfv(GL_FRONT, GL_AMBIENT, amb);
     glMaterialfv(GL_FRONT, GL_DIFFUSE, diff);
}

void lighting(){
    glShadeModel(GL_SMOOTH);
	//GLfloat lightpos[] = {-10, 10, 10, 1.0};
    //glLightfv(GL_LIGHT0, GL_POSITION, lightpos);

    //glColorMaterial ( GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE ) ;
    //glEnable ( GL_COLOR_MATERIAL ) ;


	glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT0);
	glEnable(GL_DEPTH_TEST);
}

/******************************************************************************************************************/
// Next two functions manage the mouse picking of the buttons (lower window)

void processHits2(GLint hits, GLuint buffer[]){
    unsigned int i, j;
    GLuint ii, jj, names, *ptr;

    ptr = (GLuint *) buffer;
    if(hits > 0)selection_count = 0;
    for (i = 0; i < hits; i++) {	// for each hit
        names = *ptr;
        ptr+=3;
        for (j = 0; j < names; j++) { //  for each name

            if(*ptr == play){
                cout<<"Play Button pressed"<<endl;
                loadFile();
                playButton();
            }
            if(*ptr == stop){
                cout<<"Stop Button pressed"<<endl;
                recordEnd();
                saveFile();
            }
            if(*ptr == record){
                cout<<"Record Button pressed"<<endl;
                recordBegin();
            }
            ptr++;
        }
    }
}


void mouse2(int button, int state, int x, int y){
    GLuint selectBuf[SIZE];
    GLint hits;
    GLint viewport[4];

    if (button == GLUT_LEFT_BUTTON && state == GLUT_DOWN){
        glutSetWindow(1);

        glGetIntegerv (GL_VIEWPORT, viewport);

        glSelectBuffer (SIZE, selectBuf);
        glRenderMode(GL_SELECT);

        glInitNames();
        glPushName(0);

        glPushMatrix ();
        glLoadIdentity ();
// create 5x5 pixel picking region near cursor location
        gluPickMatrix ((GLdouble) x, (GLdouble) (viewport[3] - y), 1.0, 1.0, viewport);
        ortho(Width, Height);

        drawObjects2(GL_SELECT);

        glMatrixMode (GL_PROJECTION);
        glPopMatrix ();
        glFlush ();

        hits = glRenderMode (GL_RENDER);
        processHits2 (hits, selectBuf);

        glutPostRedisplay();
    }
}


/*******************************************************************************************************************/
// Next three functions are to set initial parameters and display of the lower window

void initSecondWindow()
{
/*
    Initialize default parameters of the lower window.
*/
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    int w = glutGet(GLUT_WINDOW_WIDTH);
    int h = glutGet(GLUT_WINDOW_HEIGHT);
    //cout<<w<<h;

    if (1) {
		if (w <= h)
			gluOrtho2D(-2.0, 2.0, -2.0 * (GLfloat) h / (GLfloat) w,
				2.0 * (GLfloat) h / (GLfloat) w);
		else
			gluOrtho2D(-2.0 * (GLfloat) w / (GLfloat) h,
				2.0 * (GLfloat) w / (GLfloat) h, -2.0, 2.0);
	}

}

// drawing buttons for the lower window
void drawObjects2(GLenum mode){

    glColor3f(1,1,1);
    glPushMatrix();
    glTranslatef(0, -1.6, 0);
    glScalef(lowWinWid/(0.1), 1,1);
    glutSolidCube(0.03);
    glPopMatrix();

    // Drawing stop button
    glTranslatef(0, -1.8, 0);
    if(mode == GL_SELECT)glPushName(stop);
    glRectf(-0.05, -0.07, 0.12, 0.1);
    if(mode == GL_SELECT)glPopName();

    // Drawing record button
    if(mode == GL_SELECT)glPushName(record);
    glColor3f(1,0,0);
    glTranslatef(0.5, 0, 0);
    glutSolidSphere(0.08,15,15);
    if(mode == GL_SELECT)glPopName();

    // Drawing play button
    glTranslatef(0, 0, 0);
    glColor3f(1,1,1);
    if(mode == GL_SELECT)glPushName(play);
    glBegin(GL_TRIANGLES);
        glVertex3f(-1, -0.1, 0);
        glVertex3f(-1, 0.11, 0);
        glVertex3f(-0.8, 0.02, 0);
    glEnd();
    if(mode == GL_SELECT)glPopName();



}


void display2()
{
/* display callback, clear frame buffer and z buffer, swap buffers */

    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	//Set up viewing location
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

    drawObjects2(GL_RENDER);

    glutSwapBuffers();
	glFlush();
}


/**************************************************************************************************************/
//next function calculate the projection window

void ortho(int w, int h){

    w = glutGet(GLUT_WINDOW_WIDTH);
    h = glutGet(GLUT_WINDOW_HEIGHT);

    if (ortho_proj) {
		if (w <= h)
			glOrtho(-2.0, 2.0, -2.0 * (GLfloat) h / (GLfloat) w,
				2.0 * (GLfloat) h / (GLfloat) w, -20.0, 20.0);
		else
			glOrtho(-2.0 * (GLfloat) w / (GLfloat) h,
				2.0 * (GLfloat) w / (GLfloat) h, -2.0, 2.0, -20.0, 20.0);
	}
	else
        gluPerspective(50.0, 1.0, 3.0, 30.0);
}


/*****************************************************************************************************************/
// this function reset all the angles of the 3D character

void reset_angles(){

    theta[head] = 0.0;
    //theta[body] = 0.0;
    theta[left_arm] = 90;
    theta[left_lower_arm] = 0;
    theta[right_arm] = 90;
    theta[right_lower_arm] = 0;
    theta[left_leg] = 90;
    theta[left_lower_leg] = 0.0;
    theta[right_leg] = 90;
    theta[right_lower_leg] = 0.0;

    display();
}


/*****************************************************************************************************************/
// Menu setting

typedef struct menuEntryStruct {
    char *label;
    int key;
} menuEntryStruct;

static menuEntryStruct mainMenu[] = {
    "Body", 		    '0',
    "Head", 		    '1',
    "Left arm", 	    '2',
    "Left lower arm",   '3',
    "Right arm", 		'4',
    "Right lower arm", 	'5',
    "Left leg", 		'6',
    "Left lower leg", 	'7',
    "Right leg",    	'8',
    "Right lower leg", 	'9',

    "quit", 			27,
};
int mainMenuEntries = sizeof(mainMenu)/sizeof(menuEntryStruct);

void selectMain(int choice)
{
    userEventAction(mainMenu[choice].key);
}


static menuEntryStruct videoMenu[] = {
    "Record", 		    'r',
    "Stop record",		's',
    "Play",              'p',
};

int videoMenuEntries = sizeof(videoMenu)/sizeof(menuEntryStruct);

void selectVideo(int choice)
{
    userEventAction(videoMenu[choice].key);
}

void setMenuEntries(bool init)
{
    int i, sub;

    if (init) {
	sub = glutCreateMenu(selectVideo);
	for (i=0; i < videoMenuEntries; i++) {
	    glutAddMenuEntry(videoMenu[i].label, i);
	}
	glutCreateMenu(selectMain);
	for (i=0; i < mainMenuEntries; i++) {
	    glutAddMenuEntry(mainMenu[i].label, i);
	}
	glutAddSubMenu("Video", sub);
	glutAttachMenu(GLUT_RIGHT_BUTTON);
    } else {
    }
}

void userEventAction(int key) {
    switch(key) {
    case '0':
      cur_selectBuf[selection_count++] = body;
      break;
    case '1':
      cur_selectBuf[selection_count++] = head;
      break;
    case '2':
      cur_selectBuf[selection_count++] = left_arm;
      break;
    case '3':
      cur_selectBuf[selection_count++] = left_lower_arm;
      break;
    case '4':
        cur_selectBuf[selection_count++] = right_arm;
        break;
    case '5':
      cur_selectBuf[selection_count++] = right_lower_arm;
      break;
    case '6':
        cur_selectBuf[selection_count++] = left_leg;
      break;
    case '7':
      cur_selectBuf[selection_count++] = left_lower_leg;
      break;
    case '8':
        cur_selectBuf[selection_count++] = right_leg;
      break;
    case '9':
      cur_selectBuf[selection_count++] = right_lower_leg;

    case 'r':
      recordBegin();
      break;
    case 's':
        saveFile();
        recordEnd();
      break;
    case 'p':
        loadFile();
    	playButton();
      break;

    case 27:
      exit(0);
    default:
      break;
    }
    glutPostRedisplay();
}


int main(int argc, char **argv)
{
    glutInit(&argc, argv);
    glutInitDisplayMode (GLUT_SINGLE | GLUT_RGB|GLUT_DEPTH);
    glutInitWindowSize(Width, Height);
    glutCreateWindow("Robot");
    initSecondWindow();
	glutDisplayFunc(display2);
    glutMouseFunc (mouse2);

	/*Create a sub window at the top. */
	glutCreateSubWindow(1, 0, 0, glutGet(GLUT_WINDOW_WIDTH), lowWinHei);
    glutReshapeFunc(myReshape);
    glutDisplayFunc(display);
    glutKeyboardFunc(key);
    glutMouseFunc (mouse);
	init();
	lighting();


    glutMainLoop();
    return 0;
}






























