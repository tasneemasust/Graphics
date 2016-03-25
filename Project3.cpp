// Project#2.cpp
// OS : Ubuntu 14.0
// Language: C++


/*
    Marufa Rahmi
    Std ID: W1128039
    Course: COEN 290
*/

/* This program is the 3rd project of COEN 290
    B-spline surface

    1. generate B-Spline curve
    2. generate a B-Spline surface
    3. Texture mapping

    Keyboard inputs:
    1 - select X-axis
    2 - select Y-axis
    s - Spin around selected axis

    note: to see the shading properly its better to draw control points
    at the right side of Y-axis

*/

#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <math.h>

#include <GL/glut.h>
#include <GL/glu.h>
#include <GL/gl.h>

#define PI 2*acos(0)
#define MAX_CPTS  25
#define MAX_KNOT_VALUES MAX_CPTS+5	/* Max number of knot values. */
#define MAX_B_POINTS MAX_KNOT_VALUES*100	/* Assume that the t increment in drawBsplineCurve() will not be less than 0.01.*/

using namespace std;

int width = 500, height = 500;

struct point{
    GLfloat x,y,z;
};

struct triangle{
    point v1, v2, v3;
    point normal;
};

int gridOn = 1;			        /* Draw the background grid */
int cpolygonOn = 1;		        /* Draw the control polygon */
int BsplineOn = 1;		        /* Draw the B-spline curve. */

/* State variables */
int ctrlPointOn = 1;              /* Allowed to add control points*/
int ctrlPolygonOn = 1;            /* Displaying control polygon*/
int B_SplineCurveOn = 1;          /* Displaying B-spline curve*/
int selectPointOn = -1;           /* Not allowed to select points*/
int deletePointOn = -1;           /* Not allowed to delete any point*/
int movePointOn = -1;             /* Not allowed to move any point*/
int shadeOn = -1;
int wireframeOn = -1;
int textureMapOn = -1;


GLfloat cpts[MAX_CPTS][3];      /*array to store control points*/
int ncpts = 0;                  /*number of control points*/
int numBpts;                    /* Points in the B-spline curve. */
float B_points[MAX_B_POINTS][2];/* Array to store Points in the B-spline curve. */
float	knot[MAX_KNOT_VALUES];  /* Knot vector */

int numCurves = 0;
point triangleMesh[MAX_B_POINTS][20];

int selectedPoint = -1;
float curPos[3];

static GLfloat theta[] = {0.0,0.0,0.0};
static GLint axis = 1;

GLfloat mat_specular1[] ={0.628281, 0.555802, 0.366065, 1.0};
GLfloat mat_ambient1[] ={0.24725, 0.1995, 0.0745, 1.0};
GLfloat mat_diffuse1[] ={0.75164, 0.60648, 0.22648, 1.0};
GLfloat mat_shininess1[] ={128.0 * 0.4};

GLfloat mat_specular[] = {1.0, 1.0, 1.0, 1.0};
GLfloat mat_diffuse[] = {0.8, 0.0, 0.0, 1.0};
GLfloat mat_ambient[] = {0.2, 0.0, 0.0, 1.0};
GLfloat mat_shininess[]= {100.0};


/* Function prototypes */
void userEventAction(int key);
void SetMaterial(GLfloat spec[], GLfloat amb[], GLfloat diff[], GLfloat shin[]);

void drawGrid()
{
	int i;
	float x,  y;

		glBegin(GL_LINES);

			glColor3f(0.2, 0.2, 0.2);
			for (i=0; i<9; i++) {
				y = -10.0 + i*2.5;
				glVertex3f(-9.9, y, 0.0);
				glVertex3f(9.9, y, 0.0);
			}

			for (i=0; i<9; i++) {
				x = -10.0 + i*2.5;
				glVertex3f(x, -9.9,  0.0);
				glVertex3f(x, 9.9,  0.0);
			}
            //draw axis
			glColor3f(1.0, 0.0, 0.0);
			glVertex3f(-9.9,  0.0, 0.0);
			glVertex3f(9.9,  0.0, 0.0);

			glVertex3f(0.0, -9.9,   0.0);
			glVertex3f(0.0, 9.9,  0.0);

		glEnd();
}

void drawCtrlPolygon()
{
	int i;
	glColor3f(0.0, 1.0, 0.5);
	//printf("ncpts %d\n", ncpts);

	glBegin(GL_LINE_STRIP);
	for (i = 0; i < ncpts; i++) {
			//printf("i %d\n", i);
			glVertex3fv(cpts[i]);
	}
    glEnd();
}


/* This function implements the Cox deBoor algorithm.  */
float CoxdeBoor(int i, int p, float t)
{
	float	left, right;

	if (p==1) {
		if (knot[i] < knot[i+1] && knot[i] <= t && t < knot[i+1])
			return( 1.0 );
		else
			return( 0.0 );
	}
	else {
		if (knot[i+p-1] - knot[i] != 0.0)
			left = CoxdeBoor(i, p-1, t)*(t - knot[i])/
				(knot[i+p-1] - knot[i]);
		else
			left = 0.0;

		if (knot[i+p] - knot[i+1] != 0.0)
			right = CoxdeBoor(i+1, p-1, t)*(knot[i+p] - t)/
				(knot[i+p] - knot[i+1]);
		else
			right = 0.0;

		return( left + right );
	}
}


/* Tjis function compute and then draw the current B-Spline curve. */
void drawBsplineCurve()
{
	int		i;
	int		m = ncpts - 1;
	int		num_knots;
	float	t, B0, B1, B2, B3, x, y;


	// Compute the knot vector
		for (i=0; i<=3; i++)
			knot[i] = 0.0;

		for (i=4; i<=m; i++)
			knot[i] = i - 3.0;

		for (i=m+1;  i<=m+4; i++)
			knot[i] = m - 2.0;

		num_knots = m+4;
		printf("num_knots %d\n", num_knots);
		printf("knots: ");
		for (i=0;  i<=m+4; i++)
			printf("%.1f ", knot[i]);

		printf("\n");

        numBpts = -1;

	// Compute the store the points along the B-spline curve

		for (i=3; i < num_knots-3; i++) {
			for (t = knot[i]; t < knot[i+1]; t += 0.2) {
				B0 = CoxdeBoor(i, 4, t);
				B1 = CoxdeBoor(i-1, 4, t);
				B2 = CoxdeBoor(i-2, 4, t);
				B3 = CoxdeBoor(i-3, 4, t);

				x = cpts[i][0] * B0 +
				    cpts[i-1][0] * B1 +
				    cpts[i-2][0] * B2 +
				    cpts[i-3][0] * B3;

				y = cpts[i][1] * B0 +
				    cpts[i-1][1] * B1 +
				    cpts[i-2][1] * B2 +
				    cpts[i-3][1] * B3;

				numBpts++;
				B_points[ numBpts][0] = x;
				B_points[ numBpts][1] = y;
			}
		}

		// Store the last point of the B-spline curve
		numBpts++;
		B_points[numBpts][0] = cpts[ncpts-1][0];
		B_points[numBpts][1] = cpts[ncpts-1][1];

		glColor3f(1.0, 0.4, 0.4);

		glBegin(GL_LINE_STRIP);
		for (i = 0; i <= numBpts; i++) {
			glVertex3f(B_points[i][0],B_points[i][1], 0.0 );
		}
		glEnd();
}


/*******************************************************************************************************/
// drawing a wireframe suface by creating triangle mesh
void wireframeSurface(){
    float theta = 0, a,b,c;
    numCurves = 0;

    while(theta < 360){
        glColor3f(1, 0, 1);
        glBegin(GL_LINE_STRIP);
        //cout<<"number of points "<<numBpts<<endl;
		for (int i = 0; i <= numBpts; i++) {
            a = B_points[i][0] * cos(theta*PI/180);
            b = B_points[i][1];
            c = -B_points[i][0] * sin(theta*PI/180);
			glVertex3f(a, b, c);
			//cout<<"x = "<< a<< "y = "<< b<<"z = "	<<c<<endl;
			triangleMesh[i][numCurves].x = a;
			triangleMesh[i][numCurves].y = b;
			triangleMesh[i][numCurves].z = c;
        }
        glEnd();
        theta += 20;
        numCurves++;
    }
}


point calNormal(point v1, point v2, point v3){
    float ax, ay, az, bx, by, bz, norm;
    point normal;

    ax = v2.x - v1.x;
    ay = v2.y - v1.y;
    az = v2.z - v1.z;
    bx = v3.x - v1.x;
    by = v3.y - v1.y;
    bz = v3.z - v1.z;

    normal.x = ay*bz - by*az;
    normal.y = bx*az - ax*bz;
    normal.z = ax*by - bx*ay;
    norm = sqrt(normal.x * normal.x + normal.y * normal.y + normal.z * normal.z);

    if (norm != 0.0) {
        normal.x = normal.x/norm;
        normal.y = normal.y/norm;
        normal.z = normal.z/norm;
    }
    return normal;
}

void drawTriangleMesh(){
    int i,j;

   // if(textureMapOn == 1)glColor3f(1.0, 1.0, 1.0);
    SetMaterial(mat_specular, mat_diffuse, mat_ambient, mat_shininess);

    GLfloat texInd_x = 0.0, texInd_y=0.0, diff_x, diff_y;
    diff_x = 1.0 / numCurves;
    diff_y = 1.0 / numBpts;
    point normal;


    for(j=0; j<numCurves; j++){
        for(i=0; i<numBpts; i++){
            //glColor3f(0.8, 0.2, 0.6);

            //if(j>=numCurves/2)
            normal = calNormal(triangleMesh[i][j], triangleMesh[i+1][j],triangleMesh[i][(j+1)%numCurves]);
            //else
            //normal = calNormal(triangleMesh[i][(j+1)%numCurves], triangleMesh[i+1][j],triangleMesh[i][j]);

            if(shadeOn == -1 )glBegin(GL_LINE_STRIP);
            else glBegin(GL_TRIANGLES);

            glNormal3f(normal.x, normal.y, normal.z);
            glTexCoord2f(texInd_x, texInd_y);
            glVertex3f(triangleMesh[i][j].x, triangleMesh[i][j].y, triangleMesh[i][j].z);

            glNormal3f(normal.x, normal.y, normal.z);
            glTexCoord2f(texInd_x, texInd_y + diff_y);
            glVertex3f(triangleMesh[i+1][j].x, triangleMesh[i+1][j].y, triangleMesh[i+1][j].z);

            glNormal3f(normal.x, normal.y, normal.z);
            glTexCoord2f(texInd_x + diff_x, texInd_y);
            glVertex3f(triangleMesh[i][(j+1)%numCurves].x, triangleMesh[i][(j+1)%numCurves].y, triangleMesh[i][(j+1)%numCurves].z);
            glEnd();

           // glColor3f(0.8, 0.6, 0.2);
            //if(j>=numCurves/2)
            normal = calNormal(triangleMesh[i][(j+1)%numCurves], triangleMesh[i+1][j],triangleMesh[i+1][(j+1)%numCurves]);
            //else
            //normal = calNormal(triangleMesh[i+1][(j+1)%numCurves], triangleMesh[i+1][j],triangleMesh[i][(j+1)%numCurves]);

            if(shadeOn == -1 )glBegin(GL_LINE_STRIP);
            else glBegin(GL_TRIANGLES);

            glNormal3f(normal.x, normal.y, normal.z);
            glTexCoord2f(texInd_x + diff_x, texInd_y);
            glVertex3f(triangleMesh[i][(j+1)%numCurves].x, triangleMesh[i][(j+1)%numCurves].y, triangleMesh[i][(j+1)%numCurves].z);

            glNormal3f(normal.x, normal.y, normal.z);
            glTexCoord2f(texInd_x, texInd_y + diff_y);
            glVertex3f(triangleMesh[i+1][j].x, triangleMesh[i+1][j].y, triangleMesh[i+1][j].z);

            glNormal3f(normal.x, normal.y, normal.z);
            glTexCoord2f(texInd_x + diff_x, texInd_y + diff_y);
            glVertex3f(triangleMesh[i+1][(j+1)%numCurves].x, triangleMesh[i+1][(j+1)%numCurves].y, triangleMesh[i+1][(j+1)%numCurves].z);
            glEnd();

            texInd_y += diff_y;
            //cout<<texInd_x<< "  " << texInd_y<<endl;
        }
        texInd_y= 0.0;
        texInd_x += diff_x;
    }
    //cout<<texInd_x<< "  " << texInd_y<<endl;
}

/*******************************************************************************************************/




/*******************************************************************************************************/
// Mouse Motion

float lastPos[3] = {0.0F, 0.0F, 0.0F};

void mouseMotion(int x, int y)
{
    float  dx, dy, dz;

    float wx = (20.0 * x) / (float)(width - 1) - 10.0;
    float wy = (20.0 * (height - 1 - y)) / (float)(height - 1) - 10.0;
    lastPos[0] = curPos[0];
    lastPos[1] = curPos[1];
    lastPos[2] = curPos[2];

    curPos[0] = wx;
    curPos[1] = wy;
    curPos[2] = 0;

    if(selectPointOn == 1 && movePointOn == 1){
        cpts[selectedPoint][0] += (curPos[0] - lastPos[0]);
        cpts[selectedPoint][1] += (curPos[1] - lastPos[1]);
        cpts[selectedPoint][2] += (curPos[2] - lastPos[2]);
        glutPostRedisplay();
    }
}
// End of mouse motion
/********************************************************************************************************/


static void mouse(int button, int state, int x, int y)
{
    float wx, wy;

    /* We are only interested in left clicks */

    if (button != GLUT_LEFT_BUTTON || state != GLUT_DOWN)
        return;

    /* Translate back to our coordinate system */
    wx = (20.0 * x) / (float)(width - 1) - 10.0;
    wy = (20.0 * (height - 1 - y)) / (float)(height - 1) - 10.0;
    curPos[0] = wx;
    curPos[1] = wy;
    curPos[2] = 0;

    if(selectPointOn == 1){
        for(int i = 0; i < ncpts; i++){
            if(fabs(cpts[i][0] - wx) < 0.5 && fabs(cpts[i][1] - wy) < 0.5){
                selectedPoint = i;
            }
        }
        glutPostRedisplay();
    }
    /* See if we have room for any more control points */
    if (ctrlPointOn == 1 && ncpts != MAX_CPTS ){
        /* Save the point */
        cpts[ncpts][0] = wx;
        cpts[ncpts][1] = wy;
        cpts[ncpts][2] = 0.0;
        ncpts++;
        glutPostRedisplay();
    }
}

void spinView()
{

/* Idle callback, spin cube 2 degrees about selected axis */

	theta[axis] += 0.02;
	if( theta[axis] > 360.0 ) theta[axis] -= 360.0;
	glutPostRedisplay();
}


void SetMaterial(GLfloat spec[], GLfloat amb[], GLfloat diff[], GLfloat shin[]){
     glShadeModel(GL_SMOOTH);
     glMaterialfv(GL_FRONT, GL_SPECULAR, spec);
     glMaterialfv(GL_FRONT, GL_SHININESS, shin);
     glMaterialfv(GL_FRONT, GL_AMBIENT, amb);
     glMaterialfv(GL_FRONT, GL_DIFFUSE, diff);
}


void lighting(){

	GLfloat lightpos[] = {5, 5, 10, 1.0};
    glLightfv(GL_LIGHT0, GL_POSITION, lightpos);

    glColorMaterial ( GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE ) ;
    glEnable ( GL_COLOR_MATERIAL ) ;

    glShadeModel(GL_SMOOTH);
	glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT0);
	glEnable(GL_DEPTH_TEST);
}

void display(void)
{
    int i;
    glClear(GL_COLOR_BUFFER_BIT);
    glLoadIdentity();

	glRotatef(theta[0], 1.0, 0.0, 0.0);
	glRotatef(theta[1], 0.0, 1.0, 0.0);
	glRotatef(theta[2], 0.0, 0.0, 1.0);


	// Draw the grid if needed
	if (gridOn)
		drawGrid();

	// Draw the control polygon if needed
	if (ctrlPolygonOn == 1) {
		//printf("cpolygon %d\n", cpolygonOn);
		drawCtrlPolygon();
	}

	// Draw the B-spline curve if needed
	if (B_SplineCurveOn == 1 && ncpts>3 )
		drawBsplineCurve();

    if(wireframeOn == 1){
        wireframeSurface();
        drawTriangleMesh();
    }

	// Draw the control points
    glPointSize(5.0);
    glBegin(GL_POINTS);
    for (i = 0; i < ncpts; i++){
        if(i == selectedPoint) glColor3f(0.0, 0.0, 1.0);
        else glColor3f(0.0, 0.0, 0.0);
        glVertex3fv(cpts[i]);
    }
    glEnd();

    glFlush();
    glutSwapBuffers();
}

void keyboard(unsigned char key, int x, int y)
{

    switch (key)
    {
        case 'q': case 'Q':
            exit(0);
            break;
        case 'c': case 'C':
			ncpts = 0;
			glutPostRedisplay();
            break;
        case 's': case 'S':
            glutIdleFunc(spinView);
            break;
        case '1':
            axis = 0;
            break;
        case '2':
            axis = 1;
            break;
    }
}

/* set texture enabled or disabled*/
void userSettings(){
    if(textureMapOn == 1){
        glEnable(GL_TEXTURE_2D);
    }
    else {
        glDisable(GL_TEXTURE_2D);
    }
}


/*****************************************************************************************************************/
// Menu setting
typedef struct menuEntryStruct {
    char *label;
    int key;
} menuEntryStruct;

static menuEntryStruct mainMenu[] = {
    "Ctrl points On/Off ", 		    '0',
    "Ctrl polygon On/Off", 		    '1',
    "B-Spline Curve On/Off", 	    '2',
    "Select Ctrl point",            'a',
    "Delete Ctrl point",    		'd',
    "Move Ctrl point",              'm',
    "Save Control polygon", 	   	's',
    "Retrieve Control polygon",     'r',
    "Clear window",                	'c',
    "Draw Wireframe Surface",   	'w',
    "Shade Surface",            	'p',
    "Texture Surface",           	't',

    "quit", 			27,
};

int mainMenuEntries = sizeof(mainMenu)/sizeof(menuEntryStruct);

void selectMain(int choice)
{
    userEventAction(mainMenu[choice].key);
}

void setMenuEntries(bool init)
{
    int i, sub;

    if (init) {
	glutCreateMenu(selectMain);
	for (i=0; i < mainMenuEntries; i++) {
	    glutAddMenuEntry(mainMenu[i].label, i);
	}
	glutAttachMenu(GLUT_RIGHT_BUTTON);
    }
}

void userEventAction(int key) {
    switch(key) {
    case '0':
        ctrlPointOn = -ctrlPointOn;
        //selectPoint = -selectPoint;
        break;
    case '1':
        ctrlPolygonOn = -ctrlPolygonOn;
        break;
    case '2':
        B_SplineCurveOn = -B_SplineCurveOn;
        break;
    case 'a':
        selectPointOn = -selectPointOn;
        ctrlPointOn = -1;
        break;
    case 'd':
        deletePointOn = -deletePointOn;
        if(selectPointOn == 1 && selectedPoint > -1){
            for(int i = selectedPoint; i < ncpts - 1; i++){
                cpts[i][0] = cpts[i+1][0];
                cpts[i][1] = cpts[i+1][1];
                cpts[i][2] = cpts[i+1][2];
            }
            selectedPoint = -1;
        }
        break;
    case 'm':
        movePointOn = -movePointOn;
        break;
    case 's':
        //save the control points
        freopen ("bspline.txt","w",stdout);
        for(int i = 0; i < ncpts; i++){
            cout<<cpts[i][0]<<" "<<cpts[i][1]<<" "<<cpts[i][2]<< endl;
        }
        fclose(stdout);
        break;
    case 'r':
        // retrieve the control points
        freopen ("bspline.txt","r",stdin);
        float a,b,c;
        ncpts = 0;
        while(cin>>a>>b>>c){
            cpts[ncpts][0] = a;
            cpts[ncpts][1] = b;
            cpts[ncpts][2] = c;
            ncpts++;
        }
        fclose(stdin);
        break;
    case 'c':
        // clear the display window, delete all points
        ncpts = 0;
        numBpts = 0;
        selectPointOn = -1;
        textureMapOn = -1;
        shadeOn = -1;
        ctrlPointOn = 1;
        ctrlPolygonOn = 1;
        B_SplineCurveOn = 1;
        break;
    case 'w':
        // draw wireframe surface
        wireframeOn = -wireframeOn;
        ctrlPointOn = -1;
        break;
    case 'p':
        //shade surface
        shadeOn = -shadeOn;
        break;
    case 't':
        // texture mapping
        textureMapOn = - textureMapOn;
       break;
    case 27:
        exit(0);
    default:
        break;
    }
    userSettings();
    glutPostRedisplay();
}


/* This routine handles window resizes */
void reshape(int w, int h)
{
    width = w;
    height = h;

    /* Set the transformations */
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glOrtho(-10.0, 10.0, -10.0, 10.0, -10.0, 10.0);
    glMatrixMode(GL_MODELVIEW);
    glViewport(0, 0, w, h);
}



void init(){
/*
    Initialize default parameters.
*/
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glOrtho(-10.0, 10.0, -10.0, 10.0, -10.0, 10.0);
    glMatrixMode(GL_MODELVIEW);
    setMenuEntries(true);
   // gluLookAt(0.5, 1, 5, 0, 0, 0, 0,1 , 0);
}

#define BPP 4

GLuint tgaToTexture( char* filename, int* outWidth, int* outHeight ) {

	// open the file
	FILE* file = fopen( filename, "rb" );
	if( file == NULL ) {
		fprintf( stderr, "Could not open file %s for reading.\n", filename );
		return 0;
	}

	// skip first two bytes of data we don't need
	fseek( file, 2, SEEK_CUR );

	// read in the image type.  For our purposes the image type should
	// be either a 2 or a 3. This means that it's uncompressed (10 means compression)
	unsigned char imageTypeCode;
	fread( &imageTypeCode, 1, 1, file );
	if( imageTypeCode != 2 && imageTypeCode != 3 ) {
		fclose( file );
		fprintf( stderr, "File %s is an unsupported TGA type: %d\n", filename, imageTypeCode );
		return 0;
	}

	// skip 9 bytes of data we don't need
	fseek( file, 9, SEEK_CUR );

	// read image dimensions
	int imageWidth = 0;
	int imageHeight = 0;
	int bitCount = 0;
	fread( &imageWidth, sizeof( short ), 1, file );
	fread( &imageHeight, sizeof( short ), 1, file );
	fread( &bitCount, sizeof( unsigned char ), 1, file );

	// allocate memory for image data. It's basically a huge array of bytes
	unsigned char* bytes = (unsigned char*) malloc( (size_t) (imageWidth * imageHeight * BPP) );

	// read in data
	if( bitCount == 32 ) { // If it's a 32-bit file
		for( int i = 0; i != imageWidth * imageHeight; ++i ) {
			bytes[ i * BPP + 3 ] = fgetc( file );  // alpha
			bytes[ i * BPP + 2 ] = fgetc( file );  // blue
			bytes[ i * BPP + 1 ] = fgetc( file );  // green
			bytes[ i * BPP + 0 ] = fgetc( file );  // red
		}	// These above 4 lines are confiured for Paint.NET TGA files. I don't know why/how they're
			// different, or if there's a color code somewhere in the file that tells me how they're aligned
	} else {
		for( int i = 0; i != imageWidth * imageHeight; ++i ) {
			bytes[ i * BPP + 0 ] = fgetc( file ); // red
			bytes[ i * BPP + 1 ] = fgetc( file ); // green
			bytes[ i * BPP + 2 ] = fgetc( file ); // blue
			bytes[ i * BPP + 3 ] = 255;           // alpha
		}	// if its not 32-bit, then we're assuming its 24-bit (TGA only supports these formats, I believe)
			// since the last byte of a 32-bit pixel is alpha, for a 24-byte pixel, we just set the alpha value
			// to full opacity.
	}

	fclose( file );

	// generate an OpenGL texture from the bytes that we've loaded.
	GLuint tex;
	glGenTextures( 1, &tex );

	// build MipMaps, too.
	glBindTexture( GL_TEXTURE_2D, tex );
	gluBuild2DMipmaps( GL_TEXTURE_2D, GL_RGBA, imageWidth, imageHeight, GL_RGBA, GL_UNSIGNED_BYTE, bytes );

	// also return the width and the height to the provided variables. As I said in the header, providing NULL
	// to both of these dimensions doesn't alter function behavior. This is why.
	if( outWidth ) {
		*outWidth = imageWidth;
	}
	if( outHeight ) {
		*outHeight = imageHeight;
	}

	// let the user know the file's been loaded correctly
	//cout << "File " << fname << " (" << imageWidth << "x" << imageHeight << ") loaded successfully!" << endl;
	printf("File %s (%d x %d) loaded\n", filename, imageWidth, imageHeight);
	return tex;
}


void textureMap(){


    GLuint txtMap;
    int txtWidth, txtHeight;


    glEnable(GL_DEPTH_TEST);

    txtMap = tgaToTexture("Tex_img_Marufa.tga", &txtWidth, &txtHeight);

    glTexParameterf(GL_TEXTURE_2D,GL_TEXTURE_WRAP_S,GL_REPEAT);
    glTexParameterf(GL_TEXTURE_2D,GL_TEXTURE_WRAP_T,GL_REPEAT);
    glTexParameterf(GL_TEXTURE_2D,GL_TEXTURE_MAG_FILTER,GL_NEAREST);
    glTexParameterf(GL_TEXTURE_2D,GL_TEXTURE_MIN_FILTER,GL_NEAREST);

	// Replace the shaded object's surface with the texture pattern
	glTexEnvi(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE);


}

main(int argc, char **argv)
{
    /* Intialize the program */
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_RGB);
    glutInitWindowSize(width, height);
    glutCreateWindow("B-spline Surface");

    /* Register the callbacks */
    glutDisplayFunc(display);
    glutMouseFunc(mouse);
    glutMotionFunc(mouseMotion);
    glutKeyboardFunc(keyboard);
    glutReshapeFunc(reshape);

    init();
    lighting();
    textureMap();
    glClearColor(1.0, 1.0, 1.0, 1.0);
    glutMainLoop();
}
