// Project#4.cpp
// OS : Ubuntu 14.0
// Language: C++


/*
    Marufa Rahmi
    Std ID: W1128039
    Course: COEN 290
*/

/* This program is the 4th project of COEN 290
This program renders four spheres, three planes and a triangle mesh
using recursive ray tracing and Textures mapping.
*/

#include<iostream>
#include<stdio.h>
#include<math.h>
#include<GL/glut.h>
#include<GL/gl.h>

using namespace std;


#define	PI	2*acos(0)

/* Max image size allowed. */
#define MAX_SIZE 512

/*  Define some structures.  */

struct	points	{
    float   x, y, z;
};

typedef struct	rgb_struct	{
    float   r, g, b;
} rgb;

struct Plane{
    float a,b,c,d;
    float x_min, y_min, z_min, x_max, y_max, z_max;
    float r_ind;
};

struct Sphere{
    float cx,cy,cz;
    float r;
    float r_ind;
};

struct triangle{
    points v1, v2, v3;
    points normal;
};

/*  Viewing parameters.  */

struct	points	from, At, up;
float	VXR, VXL, VYB, VYT;
float	ax, ay, az, bx, by, bz, cx, cy, cz;
float	viewangle, angle, tanv2;
float	xinterval, yinterval;

/*  Illumination parameters.  */

points	light;
rgb	il, ia;
rgb	ka1, kd1, ks1;
rgb	ka2, kd2, ks2;
rgb	ka3, kd3, ks3;
rgb	ka4, kd4, ks4;
rgb	ka5, kd5, ks5;
rgb	ka6, kd6, ks6;
rgb	ka7, kd7, ks7;
rgb	tka1, tkd1, tka2, tkd2, tka3, tkd3;
rgb	tka4, tkd4, tka5, tkd5, tka6, tkd6;
rgb tka7, tkd7;

int	phong1, phong2, phong3, phong4, phong5, phong6;

/*  Image parameters.  */

int		xmax_pixel, ymax_pixel;

/* Image buffer.  A more efficient approach is to use one single array (texture_RGB)
   rather than using the three rays.  */

float *texture_R;
float *texture_G;
float *texture_B;
float noise_tabl[65][65];


/*  Object parameters.  */

Plane p_arr[20];
int numPlane = 0;
Sphere s_arr[10];
int numSphere = 0;

/*triangle mesh parameters*/
float surf_x_center;
float surf_y_center;
float surf_z_center;
int minx ,maxx, miny,maxy,minz,maxz;
triangle tri_arr[10000];
int numtri;


/* Image output file. */

FILE	*outpfile;


/*******************************************************************************
  Title:	Read_Information
  Purpose:	This function reads in the information about the objects (three spheres
            and three planes) to be rendered.  The information is
			stored in an ASCII text file called "in_MR.dat".

*******************************************************************************/

int Read_Information()
{
	string str;

	if ( freopen("in_MR.dat","r", stdin) == NULL) {
		printf("ERROR: Could not open in_MR.dat for read!\n");
		return(0);
	}

/*  Read in viewing information.  */

    while(cin>>str){

        if(str == "\n")continue;

        else if(str == "From" || str == "from"){
            cin>>from.x>>from.y>>from.z;
        }
        else if(str == "At" || str == "at"){
            cin>>At.x>>At.y>>At.z;
        }
        else if(str == "Up" || str == "up"){
            cin>>up.x>>up.y>>up.z;
        }
        else if(str == "ViewAngle" || str == "viewangle"){
            cin>>viewangle;
            angle = viewangle * PI/180.0;
            tanv2 = tan(angle/2.0);
        }
        else if(str == "Viewport" || str == "viewport"){
            cin>>VXL>>VXR>>VYB>>VYT;
        }
        else if(str == "Light" || str == "light"){
            cin>>light.x>>light.y>>light.z;
        }
        else if(str == "Plane" || str == "plane"){
            cin>>p_arr[numPlane].a>>p_arr[numPlane].b>>p_arr[numPlane].c;
            cin>>p_arr[numPlane].d;
            cin>>p_arr[numPlane].x_min>>p_arr[numPlane].y_min>>p_arr[numPlane].z_min;
            cin>>p_arr[numPlane].x_max>>p_arr[numPlane].y_max>>p_arr[numPlane].z_max;
            cin>>p_arr[numPlane].r_ind;
            numPlane++;
        }
        else if(str == "Sphere" || str == "sphere"){
            cin>>s_arr[numSphere].cx>>s_arr[numSphere].cy>>s_arr[numSphere].cz>>s_arr[numSphere].r;
            cin>>s_arr[numSphere].r_ind;
            numSphere++;
        }
        else if(str == "ImageSize" || str == "imagesize"){
            cin>>xmax_pixel>>ymax_pixel;
            if (xmax_pixel > MAX_SIZE || ymax_pixel > MAX_SIZE) {
                printf("Error: Exceeded max image size %d x %d\n", xmax_pixel, ymax_pixel);
            	printf("Reset to max image size: %d x %d\n", MAX_SIZE, MAX_SIZE);
            	xmax_pixel = MAX_SIZE-1;
            	ymax_pixel = MAX_SIZE - 1;
            }
        }
        else if(str == "Mesh"){

            cin>>surf_x_center>>surf_y_center>>surf_z_center;
            cin>>minx>>maxx>>miny>>maxy>>minz>>maxz;
            cin>>numtri;
            for(int i=0; i<numtri; i++){
                cin>>tri_arr[i].v1.x>>tri_arr[i].v1.y>>tri_arr[i].v1.z;
                cin>>tri_arr[i].v2.x>>tri_arr[i].v2.y>>tri_arr[i].v2.z;
                cin>>tri_arr[i].v3.x>>tri_arr[i].v3.y>>tri_arr[i].v3.z;
                cin>>tri_arr[i].normal.x>>tri_arr[i].normal.y>>tri_arr[i].normal.z;
            }
        }
    }
    fclose(stdin);

    printf("From: %f %f %f\n", from.x, from.y, from.z);
	printf("At: %f %f %f\n", At.x, At.y, At.z);
	printf("Up: %f %f %f \n View Angle: %f \n", up.x, up.y, up.z, viewangle);
	cout<<"x="<<surf_x_center<<endl;
    cout<<surf_y_center<<endl;
    cout<<surf_z_center<<endl;
    cout<<minx<<endl;
    cout<<maxx<<endl;
    cout<<miny<<endl;
    cout<<maxy<<endl;
    cout<<minz<<endl;
    cout<<maxz<<endl;
    cout<<"num "<<numtri<<endl;


/*  Open an output file to store the intensity values of the output image.  */

	if ((outpfile = fopen("image_MR.out","wb")) == NULL) {
		printf("ERROR:  cannot open image.out for write.\n");
		return(0);
	}

/*  Allocate memory for the image buffer.  */

	texture_R = new float [xmax_pixel * ymax_pixel ];
	texture_G = new float [xmax_pixel * ymax_pixel ];
	texture_B = new float [xmax_pixel * ymax_pixel ];
	printf("image_buf allocated.  Image size %d x %d\n", xmax_pixel, ymax_pixel);

	return(1);
}


/*******************************************************************************
  Title:	Normalize
  Purpose:	This function normalizes the given vector.

*******************************************************************************/

void Normalize(float *x,float *y,float *z)
{
	float	norm;

	norm = sqrt( *x * *x + *y * *y + *z * *z );
	if (norm != 0.0) {
		*x = *x / norm;
		*y = *y / norm;
		*z = *z / norm;
	}
}

/*******************************************************************************
  Title:	Power
  Purpose:	This function computes the power of the given base and
		    exponent.

*******************************************************************************/

float 	Power(float base,int exp)
{
	int	i;
	float	value;

	value = 1.0;
	for (i=1; i<=exp; i++)
		value *= base;

	return( value );
}


/*******************************************************************************
  Title:	Compute_M
  Purpose:	This function computes the transformation matrix to be used
		    in the perspective viewing model.

*******************************************************************************/

void Compute_M()
{

/*  Compute the line-of-sight vector, c.  */

	cx = At.x - from.x;
	cy = At.y - from.y;
	cz = At.z - from.z;
	Normalize(&cx, &cy, &cz);

/*  Compute the cross product of vector c and the up vector.  */

	ax = cy*up.z - up.y*cz;
	ay = up.x*cz - cx*up.z;
	az = cx*up.y - up.x*cy;
	Normalize(&ax, &ay, &az);

/*  Compute the cross product of vector a and c.  */

	bx = ay*cz - cy*az;
	by = cx*az - ax*cz;
	bz = ax*cy - cx*ay;
}

/*******************************************************************************
  Title:	Setup_Parameters
  Purpose:	This function sets up the necessary parameters for
		    performing the ray trace.  It first computes the
		    transformation matrix for the perspective viewing model, then
		    sets up the default illumination parameters.

*******************************************************************************/

void Setup_Parameters()
{

/*  Compute the transformation matrix for converting world coordinates to eye
    coordinates.  */

	Compute_M();

/*  Normalized the given directional light vector.

	Note:  DO NOT normalize this vector, if the position of the light source is given.
	       The light vector (L) would need to computed for each intersction point,
		   then normalized in Compute_Color().   */

	Normalize(&light.x, &light.y, &light.z);
	printf("light %f %f %f\n", light.x, light.y, light.z);

/*  Set up the conversion factors for converting from pixel coordinates to
    view port coordinates.  */

	xinterval = (VXR - VXL) / xmax_pixel;
	yinterval = (VYT - VYB) / ymax_pixel;

/*  Set up default illumination (Phong lighting) parameters.  */

	il.r = 1.0;	il.g = 1.0;	il.b = 1.0;
	ia.r = 1.0;	ia.g = 1.0;	ia.b = 1.0;

/*  Phone lighting parameters for the three spheres.  */

	ka1.r = 0.3;	ka1.g = 0.0;	ka1.b = 0.0;
	kd1.r = 0.7;	kd1.g = 0.0;	kd1.b = 0.0;
	ks1.r = 1.0;	ks1.g = 1.0;	ks1.b = 1.0;
	tka1.r = 0.2;	tka1.g = 0.0;	tka1.b = 0.0;
	tkd1.r = 0.2;	tkd1.g = 0.0;	tkd1.b = 0.0;
	phong1 = 60;

	ka2.r = 0.0;	ka2.g = 0.3;	ka2.b = 0.0;
	kd2.r = 0.0;	kd2.g = 0.7;	kd2.b = 0.0;
	ks2.r = 1.0;	ks2.g = 1.0;	ks2.b = 1.0;
	tka2.r = 0.0;	tka2.g = 0.2;	tka2.b = 0.0;
	tkd2.r = 0.0;	tkd2.g = 0.2;	tkd2.b = 0.0;
	phong2 = 90;

	ka3.r = 0.1;	ka3.g = 0.3;	ka3.b = 0.2;
	kd3.r = 0.0;	kd3.g = 0.2;	kd3.b = 0.3;
	ks3.r = 1.0;	ks3.g = 1.0;	ks3.b = 1.0;
	tka3.r = 0.0;	tka3.g = 0.0;	tka3.b = 0.2;
	tkd3.r = 0.0;	tkd3.g = 0.0;	tkd3.b = 0.2;
	phong3 = 120;

/*  Phone lighting parameters for the three planes (not shown).  */

	ka4.r = 0.1;	ka4.g = 0.1;	ka4.b = 0.0;
	kd4.r = 0.7;	kd4.g = 0.7;	kd4.b = 0.0;
	ks4.r = 1.0;	ks4.g = 1.0;	ks4.b = 1.0;
	tka4.r = 0.1;	tka4.g = 0.0;	tka4.b = 0.0;
	tkd4.r = 0.7;	tkd4.g = 0.0;	tkd4.b = 0.0;
	phong4 = 120;

	ka5.r = 0.1;	ka5.g = 0.0;	ka5.b = 0.1;
	kd5.r = 0.7;	kd5.g = 0.0;	kd5.b = 0.7;
	ks5.r = 1.0;	ks5.g = 1.0;	ks5.b = 1.0;
	tka5.r = 0.0;	tka5.g = 0.0;	tka5.b = 0.1;
	tkd5.r = 0.0;	tkd5.g = 0.0;	tkd5.b = 0.7;
	phong5 = 120;


	ka6.r = 1.0;	ka6.g = 1.0;	ka6.b = 1.0;
	kd6.r = 0.8;	kd6.g = 0.8;	kd6.b = 0.8;
	ks6.r = 0.8;	ks6.g = 0.8;	ks6.b = 0.8;
	tka6.r = 0.0;	tka6.g = 0.0;	tka6.b = 0.0;
	tkd6.r = 0.0;	tkd6.g = 0.0;	tkd6.b = 0.0;
	phong6 = 120;



	ka7.r = 0.24725;	ka7.g = 0.1995;	ka7.b = 0.0745;
	kd7.r = 0.75164;	kd7.g = 0.60648;	kd7.b = 0.22648;
	ks7.r = 0.628281;	ks7.g = 0.555802;	ks7.b = 0.366065;
	tka7.r = 0.0;	tka7.g = 0.0;	tka7.b = 0.0;
	tkd7.r = 0.0;	tkd7.g = 0.0;	tkd7.b = 0.0;
}

int TABLE_SIZE_1 = 60;
int TABLE_SIZE = 60;

/*******************************************************************************
  Title:	calc_noise
  Purpose:	This function computes the noise for Bump map.

*******************************************************************************/

float calc_noise(float iu, float iv, int direction)
{
    int i, j, x, y, left, right;
    float noise, u, v, w, u1, v1, w1;
    float val;

    //Filling noise table with random noise
    for (i = 0; i < 65; i++)
        for (j = 0; j < 65; j++)
        {
            val = rand() % 255;
            noise_tabl[i][j] = val / 256;
        }

    i = (int)iu;
    j = (int)iv;
    x = i%TABLE_SIZE;
    y = j%TABLE_SIZE;

    if (direction == 1)
    {
        if (x <= 0)
            left = 0;
        else
            left = x - 1;

        if (x >= TABLE_SIZE_1)
            right = TABLE_SIZE_1;
        else
            right = x + 1;
        //Getting the value of noise
        noise = (noise_tabl[right][y] - noise_tabl[left][y]) / 2.0;
    }

    else {
        if (y <= 0)
            left = 0;
        else
            left = y - 1;

        if (y >= TABLE_SIZE_1)
            right = TABLE_SIZE_1;
        else
            right = y + 1;
        //Getting the value of random noise
        noise = (noise_tabl[x][right] - noise_tabl[x][left]) / 2.0;
    }
    return(noise);

}

/*******************************************************************************
  Title:	Bump_Map
  Purpose:	This function computes the new Normal vector for Bump map texture.

*******************************************************************************/

void Bump_Map(float x, float y, float z, float xc, float yc, float zc, float r, float *nx, float *ny, float *nz, int methd)
{
    //Variables
    float xp, yp, zp, iu, iv, xu, yu, zu, xv, yv, zv;
    float fu, fv, a, dx, dy, dz, u, v, nnx, nny, nnz;

    xp = (x - xc) / r;
    yp = (y - yc) / r;
    zp = (z - zc) / r;

    u = asin(zp);
    v = atan2(yp, xp);

    //Getting values to generate the texture
    iu = (u + PI) / (2 * PI) * (TABLE_SIZE_1);
    iv = (v + (PI / 2)) / PI * (TABLE_SIZE_1);

    xu = -r*cos(v)*sin(u);
    xv = -r*sin(v)*cos(u);
    yu = r*cos(v)*cos(u);
    yv = -r*sin(v)*sin(u);
    zu = 0.0;
    zv = r*cos(v);

    fu = calc_noise(iu, iv, 1);
    fv = calc_noise(iu, iv, 2);

    if(methd == 1){
        dx = fu*xu + fv*xv;
        dy = fu*yu + fv*yv;
        dz = fu*zu + fv*zv;
    }

    else{
        nnx = *nx;
        nny = *ny;
        nnz = *nz;
        Normalize(&nnx, &nny, &nnz);

    //Displacing normals
        dx = fv* (yu*nnz - nny*zu);
        dy = fv* (nnx*zu - xu*nnz);
        dz = fv* (xu*nny - nnx*yu);

        dx += fu* (yv*nnz - nny*zv);
        dy += fu* (nnx*zv - xv*nnz);
        dz += fu* (xv*nny - nnx*yv);
    }
    Normalize(&dx, &dy, &dz);

    a = sqrt(fu*fu + fv*fv);
    dx *= a;
    dy *= a;
    dz *= a;

    //Finding new normals after displacement
    *nx += dx;
    *ny += dy;
    *nz += dz;
}



/*******************************************************************************
  Title:	Check_Sphere
  Purpose:	This function determines if the give ray intercepts the given
		    sphere.
*******************************************************************************/

void Check_Sphere(float px,float py,float pz,float dx,float dy,float dz,float xc,float yc,float zc,float r,
float *t1, float *t2)
{
	float	a, b, c, xdiff, ydiff, zdiff, discr;

	xdiff = px-xc;
	ydiff = py-yc;
	zdiff = pz-zc;
	a = xdiff*xdiff + ydiff*ydiff + zdiff*zdiff - r*r;
	b = 2.0*( dx*xdiff + dy*ydiff + dz*zdiff );
	c = dx*dx + dy*dy + dz*dz;

/*  Check if there are any intersections.  */

	discr = b*b - 4.0*a*c;
	if (discr < 0.0) {
		*t1 = -1.0;
		*t2 = -1.0;
	}
	else if (discr == 0.0) {
		*t1 = -b / (2.0*c);
		*t2 = -1.0;
	}
	else {
		discr = sqrt(discr);
		*t1 = (-b + discr) / (2.0*c);
		*t2 = (-b - discr) / (2.0*c);
	}
}


/*******************************************************************************
  Title:	Check_Plane
  Purpose:	This function checks if the given ray intercepts the given
		    plane.
*******************************************************************************/

void Check_Plane(float px,float py,float pz,float dx,float dy,float dz,float a,float b,float c,float d,float *t1)
{
	*t1 = (-a*px - b*py - c*pz - d) / (a*dx + b*dy + c*dz);
}

/*******************************************************************************
  Title:	Check_Shadow
  Purpose:	This function checks if there is any object between light and
   the point
*******************************************************************************/

int Check_Shadow(float px, float py, float pz){

    float i, t1, t2;

    for(int i = 0; i<numSphere; i++){
        Check_Sphere(px, py, pz, light.x, light.y, light.z, s_arr[i].cx, s_arr[i].cy, s_arr[i].cz , s_arr[i].r, &t1, &t2);
        //cout<<px<<" "<<py<<" "<<pz<<" "<<disc<<" "<< s_arr[i].cx
         //<<" "<<s_arr[i].cy<<" "<<s_arr[i].cz<<" "<<s_arr[i].r<<endl;
        if(t1>= 0.0001 || t2>=0.0001)return 1;
    }
//    for(i=0; i<numPlane; i++){
//
//        Check_Plane(px, px, py, light.x, light.y, light.z, p_arr[i].a, p_arr[i].b, p_arr[i].c, p_arr[i].d, &t1);
//
//        if (t1 >= 0.0001) return 1;
//    }
    return 0;
}

/*******************************************************************************
  Title:	dot
  Purpose:	This function computes dot product of two vectors

*******************************************************************************/

float dot(float v0[3], float v1[3])
{
	return( v0[0]*v1[0] + v0[1]*v1[1] + v0[2]*v1[2] );
}

/*******************************************************************************
  Title:	Check_Triangle
  Purpose:   Given the triangle number (triNumber), the point (px,py,pz), the ray
    described by (dx,dy,dz), and the equation of the plane described by (a,b,c,d).
    Return t1.  Set t1 = -1 if the ray does not intersect the plane or the
    intersection point is not inside the triangle triNumber.

*******************************************************************************/

void Check_Triangle(int triNumber, float px,float py,float pz,float dx,float dy,float dz,
float a,float b,float c,float d,float *t1)
{
	int	 vert1, vert2, vert3;
	float dot00, dot01, dot02, dot11, dot12, invDenom;
	float t, ipx, ipy, ipz, u, v, p[3][3], v0[3], v1[3], v2[3];

	t = (-a*px - b*py - c*pz - d) / (a*dx + b*dy + c*dz);
	if (t > 0.00001) {

		// Check if the intersection point is inside the triangle.

			ipx = px + t*dx;
			ipy = py + t*dy;
			ipz = pz + t*dz;


			p[0][0] = tri_arr[triNumber].v1.x;
			p[0][1] = tri_arr[triNumber].v1.y;
			p[0][2] = tri_arr[triNumber].v1.z;

			p[1][0] = tri_arr[triNumber].v2.x;
			p[1][1] = tri_arr[triNumber].v2.y;
			p[1][2] = tri_arr[triNumber].v2.z;

			p[2][0] = tri_arr[triNumber].v3.x;
			p[2][1] = tri_arr[triNumber].v3.y;
			p[2][2] = tri_arr[triNumber].v3.z;

			v0[0] = p[2][0] - p[0][0];
			v0[1] = p[2][1] - p[0][1];
			v0[2] = p[2][2] - p[0][2];

			v1[0] = p[1][0] - p[0][0];
			v1[1] = p[1][1] - p[0][1];
			v1[2] = p[1][2] - p[0][2];

			v2[0] = ipx - p[0][0];
			v2[1] = ipy - p[0][1];
			v2[2] = ipz - p[0][2];

			dot00 = dot(v0, v0);
			dot01 = dot(v0, v1);
			dot02 = dot(v0, v2);
			dot11 = dot(v1, v1);
			dot12 = dot(v1, v2);

			invDenom = 1.0/( dot00*dot11 - dot01*dot01);

			u = (dot11 * dot02 - dot01 * dot12) * invDenom;
			v = (dot00 * dot12 - dot01 * dot02) * invDenom;

			if ( u>=0.0 && v >= 0.0 && u+v <= 1.0)
				*t1 = t;
			else
				*t1 = -1;
	}
	else
		*t1 = -1;
}

/*******************************************************************************
  Title:	Compute_Intersection
  Purpose:	This function computes the intersection of ray with an
		    object.  The intersection point is given by a parametric value
		    t, where ray = p + d*t, d = the direction of the ray, and p is
		    the starting point of the ray.

*******************************************************************************/

void Compute_Intersection(float px,float py,float pz,float dx,float dy, float dz,float t,float *newx,float *newy,float *newz)
{
	*newx = px + t*dx;
	*newy = py + t*dy;
	*newz = pz + t*dz;
}


/*******************************************************************************
  Title:	Compute_Color
  Purpose:	This function computes the intensity of the color for the
		    given location based on the Phong lighting model.

*******************************************************************************/

void Compute_Color(int shadow_flag, float ipx,float ipy,float  ipz,float  nx,float  ny,float  nz,
				   rgb ia,rgb ka,rgb kd, rgb ks,int n,float *r,float *g, float *b)
{
	float	vx, vy, vz, rx, ry, rz;
	float	ndotl, vdotr, cosalphapower;

/*  Compute the view vector.  */

	vx = from.x - ipx;
	vy = from.y - ipy;
	vz = from.z - ipz;
	Normalize(&vx, &vy, &vz);

/*  Compute the R (reflection) vector.  */

	ndotl = nx*light.x + ny*light.y + nz*light.z;
	rx = 2.0*ndotl*nx - light.x;
	ry = 2.0*ndotl*ny - light.y;
	rz = 2.0*ndotl*nz - light.z;

/* Compute the V (view) vector. */

	vdotr = vx*rx + vy*ry + vz*rz;

/* Compute Ia * Ka.  */

	*r = ia.r * ka.r;
	*g = ia.g * ka.g;
	*b = ia.b * ka.b;

/* Compute diffuse reflection. */

	if (ndotl >= 0.0 && shadow_flag==0) {

		/*  diffuse reflection = kd * N dot L * Il  */
		*r = *r + kd.r*ndotl*il.r;
		*g = *g + kd.g*ndotl*il.g;
		*b = *b + kd.b*ndotl*il.b;

		if (vdotr >= 0.0) {

			/*  specular reflection = ks * cos(alpha)**K^n * Il */
			cosalphapower = Power(vdotr, n);

			*r = *r + ks.r*cosalphapower*il.r;
			*g = *g + ks.g*cosalphapower*il.g;
			*b = *b + ks.b*cosalphapower*il.b;
		}
	}

/*  Make sure that the color is within range.  */
	if (*r > 1.0) *r = 1.0;
	if (*g > 1.0) *g = 1.0;
	if (*b > 1.0) *b = 1.0;
}





/*******************************************************************************
  Title:	get_reflected_ray
  Purpose:	This function generates reflected ray of a given ray.

*******************************************************************************/

void get_reflected_ray(float px, float py, float pz, float dx, float dy, float dz, float nx, float ny, float nz,
    float *rlx, float *rly, float *rlz){
    /*
    c1 = -dot_product( N, V )
    Rl = V + (2 * N * c1 )
    */
    //Normalize(&nx, &ny, &nz);
    float ix = dx - px;
    float iy = dy - py;
    float iz = dz - pz;
    float c1 = -(nx*dx + ny*dy + nz*dz);
    *rlx = dx + ( 2 * nx * c1);
    *rly = dy + ( 2 * ny * c1);
    *rlz = dz + ( 2 * nz * c1);

}

/*******************************************************************************
  Title:	get_refracted_ray
  Purpose:	This function generates refracted ray of a given ray.

*******************************************************************************/

void get_refracted_ray(float dx, float dy, float dz, float nx, float ny, float nz,
    float *rrx, float *rry, float *rrz, float n2, float n1){
    /*
    n1 = index of refraction of original medium
    n2 = index of refraction of new medium
    n = n1 / n2
    c2 = sqrt( 1 - n2 * (1 - c12) )
    Rr = (n * V) + (n * c1 - c2) * N */
    float n, c1, c2, rr;
    //n1 = 1;
    //n2 = 1.33;

    n = n1 / n2;
    c1 = - ((nx*dx)+(ny*dy)+(nz*dz));
    c2 = sqrt( 1 - (n * n) * (1 - (c1 * c1)));
    *rrx = (n * dx) + (n * c1 - c2) * nx;
    *rry = (n * dy) + (n * c1 - c2) * ny;
    *rrz = (n * dz) + (n * c1 - c2) * nz;
}


/*******************************************************************************
  Title:	Ray_Tracer_recur
  Purpose:	This is a recursive function to compute color for every pixel with
            reflection and refraction.
*******************************************************************************/

rgb Ray_Tracer_recur(int level, points frm, float dx, float dy, float dz, float r_ind, int obj_type, int obj_id){

    int	    obj, obj_num, shadow_flag;
	int	    texture, buf_ptr;
	float	nx, ny, nz;
	float	t_min, t1, t2, t, ipx, ipy, ipz;
	float   rlx, rly, rlz, rrx, rry, rrz;
	float	r, g, b;
	int     i,j;
	rgb     newcolor, reflect_color, refract_color;
	float   refraction_index = r_ind;

    t_min = 999.0;
    obj_num = 0;
    obj = 0;
	texture = 0;

/*  Check if the current ray intercepts spheres  */

    for(i = 0; i<numSphere; i++){
        if(obj_type == 1 && obj_id == i)continue;
        Check_Sphere(frm.x, frm.y, frm.z, dx, dy, dz, s_arr[i].cx, s_arr[i].cy, s_arr[i].cz ,s_arr[i].r, &t1, &t2);

        if (t1>=0.0) {
            t_min = t1;
            obj = 1;
            obj_num = i;
            refraction_index = s_arr[i].r_ind;
            Compute_Intersection(frm.x, frm.y, frm.z, dx, dy, dz, t1, &ipx, &ipy, &ipz);
        }

        if (t2>=0.0 && t2<t_min) {
            t_min = t2;
            obj = 1;
            obj_num = i;
            refraction_index = s_arr[i].r_ind;
            Compute_Intersection(frm.x, frm.y, frm.z, dx, dy, dz, t2, &ipx, &ipy, &ipz);
        }

    }

/*  Check if the current ray intercepts any of the planes  */

    for(i=0; i<numPlane; i++){
        if(obj_type == 2 && obj_id == i)continue;
        Check_Plane(frm.x, frm.y, frm.z, dx, dy, dz, p_arr[i].a, p_arr[i].b, p_arr[i].c, p_arr[i].d, &t1);

        if (t1 >= 0.0 && t1<t_min) {
        /*  Check if the intersection point is inside the min/max values. */

            Compute_Intersection(frm.x, frm.y, frm.z, dx, dy, dz, t1, &ipx, &ipy, &ipz);

            if (ipx >= p_arr[i].x_min && ipx <= p_arr[i].x_max &&
                ipy >= p_arr[i].y_min  && ipy <= p_arr[i].y_max &&
                ipz >=  p_arr[i].z_min && ipz <= p_arr[i].z_max ) {

                t_min = t1;
                obj_num = i;
                obj = 2;
                refraction_index = p_arr[i].r_ind;
            }
        }
    }

    //test to see if the ray intersects the bounding sphere of the B-Spline surface
    float r_bs;
    if((maxx - surf_x_center)> (maxy - surf_y_center))r_bs = maxx - surf_x_center;
    else r_bs = maxy - surf_y_center;
    Check_Sphere(frm.x, frm.y, frm.z, dx, dy, dz, surf_x_center, surf_y_center, surf_z_center,
    r_bs, &t1, &t2);

// Perform the following ray intersection tests if the current ray intersects the bounding sphere of the B-spline surface mesh.

    if((t1 >= 0.0 && t1 < t_min)||(t2 >= 0.0 && t2< t_min)){
        if(obj_type != 3){
            float a,d,c,d0;
            for (i=0; i<numtri; i++) {

                a = tri_arr[i].normal.x;
                d = tri_arr[i].normal.y;
                c = tri_arr[i].normal.z;

                d0 = -a * tri_arr[i].v1.x - d * tri_arr[i].v1.y - c * tri_arr[i].v1.z;
		// Check if the ray intersects the plane containing the triangle.

                Check_Triangle(i, frm.x, frm.y, frm.z, dx, dy, dz, a, d, c, d0, &t);
                if (t > 0.0 && t < t_min) {

                    Compute_Intersection(frm.x, frm.y, frm.z,dx, dy, dz, t, &ipx, &ipy, &ipz);
                    t_min = t;
                    obj_num = i;
                    obj = 3;
                }
            }
        }
	}

    if(fabs(frm.x - ipx) < 0.00001 && fabs(frm.y - ipy) < 0.00001 && fabs(frm.z - ipz) < 0.00001)
        obj = 0;

/*  Compute the intensity to use at the current pixel.  */

        switch (obj) {
/*  The current ray does not intersect any of the objects.  */

            case 0 :
                r = 0.0;
                g = 0.2;
                b = 0.2;
                break;

/*  The current ray intercept spheres.  */

            case 1 :
                nx = ipx - s_arr[obj_num].cx;
                ny = ipy - s_arr[obj_num].cy;
                nz = ipz - s_arr[obj_num].cz;
                Normalize(&nx, &ny, &nz);

                shadow_flag = 0;
                shadow_flag = Check_Shadow(ipx, ipy, ipz );
                texture = 0;

                if(obj_num == 0){
                    Compute_Color(shadow_flag, ipx, ipy, ipz, nx, ny, nz, ia, ka1, kd1, ks1, phong1, &r, &g, &b);
                }
                if(obj_num == 1){
                    Compute_Color(shadow_flag, ipx, ipy, ipz, nx, ny, nz, ia, ka2, kd2, ks2, phong2, &r, &g, &b);
                }
                if(obj_num == 2){
                    Compute_Color(shadow_flag, ipx, ipy, ipz, nx, ny, nz, ia, ka3, kd3, ks3, phong3, &r, &g, &b);
                }
                if(obj_num == 3){
                    if(rand()%3 == 1)Bump_Map(ipx, ipy, ipz, dx, dy, dz, 1, &nx, &ny, &nz,1);
                    Compute_Color(shadow_flag, ipx, ipy, ipz, nx, ny, nz, ia, ka7, kd7, ks7, phong6, &r, &g, &b);
                }
                break;


            case 2 :
                nx = p_arr[obj_num].a;
                ny = p_arr[obj_num].b;
                nz = p_arr[obj_num].c;

                shadow_flag = 0;
                shadow_flag = Check_Shadow( ipx, ipy, ipz );

                if(obj_num == 0){

                    if((((int)ipx % 4) < 2 && ((int)ipz % 4) < 2)|| (((int)ipx % 4) >= 2 && ((int)ipz % 4) >= 2)){
                        texture = 0;
                    }
                    else texture = 1;

                    if (texture == 1){
                            Compute_Color(shadow_flag, ipx, ipy, ipz, nx, ny, nz, ia, tka5, tkd5, ks5, phong6, &r, &g, &b);
                    }
                    else{
                        Compute_Color(shadow_flag, ipx, ipy, ipz, nx, ny, nz, ia, ka6, kd6, ks6, phong6, &r, &g, &b);
                    }
                }
                if(obj_num == 1){
                    //calculating wood grain
                    float radius1, radius2, ang;
                    float u, v, w;
                    int grain;
                    u = (p_arr[1].x_max - p_arr[1].x_min)/2 + ipx;
                    v = (p_arr[1].y_max - p_arr[1].y_min)/2 + ipy;
                    w = (p_arr[1].z_max - p_arr[1].z_min)/2 + ipz;
                    radius1 = sqrt(u*u + v * v);
                    ang = atan(u/v);
                    radius2 = radius1 + 2 * sin(20 * ang) + u / 15;
                    grain = (int)radius2 % 4;
                    if(grain < 3 )
                        Compute_Color(shadow_flag, ipx, ipy, ipz, nx, ny, nz, ia, ka4, kd4, ks4, phong4, &r, &g, &b);
                    else
                        Compute_Color(shadow_flag, ipx, ipy, ipz, nx, ny, nz, ia, tka6, tkd6, ks6, phong4, &r, &g, &b);

                }
                if(obj_num  == 2){
                    if(rand() % 2 == 1){
                        Bump_Map(ipx, ipy, ipz, dx, dy, dz, 1, &nx, &ny, &nz,1);
                        Compute_Color(shadow_flag, ipx, ipy, ipz, nx, ny, nz, ia, ka3, kd3, ks3, phong3, &r, &g, &b);
                    }
                    else
                        Compute_Color(shadow_flag, ipx, ipy, ipz, nx, ny, nz, ia, ka3, kd3, ks3, phong3, &r, &g, &b);
                }
                break;
            case 3:
                    //color for B-Spline
                nx = tri_arr[obj_num].normal.x;
                ny = tri_arr[obj_num].normal.y;
                nz = tri_arr[obj_num].normal.z;
                Normalize(&nx, &ny, &nz);
                shadow_flag = 0;
                shadow_flag = Check_Shadow( ipx, ipy, ipz );

                Compute_Color(shadow_flag, ipx, ipy, ipz, nx, ny, nz, ia, ka1, kd1, ks1, phong6, &r, &g, &b);
                break;
			}

        /* generating reflected and refracted ray*/
            get_reflected_ray(ipx,ipy,ipz, dx, dy, dz, nx, ny, nz, &rlx, &rly, &rlz);
            if(r_ind > 1) get_refracted_ray(dx, dy, dz, nx, ny, nz, &rrx, &rry, &rrz, r_ind, 1);
            else get_refracted_ray(dx, dy, dz, nx, ny, nz, &rrx, &rry, &rrz, r_ind, refraction_index);

			points newFrom;
            newFrom.x = ipx;
            newFrom.y = ipy;
            newFrom.z = ipz;


        /* Base case: if depth is 3 or no intersection return*/

        if(level >= 4 || obj == 0){
                newcolor.r = r ;
                newcolor.g = g ;
                newcolor.b = b ;
                return newcolor;
        }
        level += 1;


        if( obj == 1 && obj_num == 2){
			refract_color = Ray_Tracer_recur(level, newFrom, rrx, rry, rrz, refraction_index,obj, obj_num);
			newcolor.r = r*0.1 +  refract_color.r/(1.1*level);
			newcolor.g = g*0.1 +  refract_color.g/(1.1*level);
			newcolor.b = b*0.1 +  refract_color.r/(1.1*level);
        }
        else if(obj == 2 && obj_num == 2){
            newcolor.r = r ;
			newcolor.g = g ;
			newcolor.b = b;
        }
//        else if(obj == 2 && obj_num == 1){
//			reflect_color = Ray_Tracer_recur(level, newFrom, rlx, rly, rlz, 1, obj, obj_num);
//
//			newcolor.r = r * 0.9 + reflect_color.r/(2*level) ;
//			newcolor.g = g * 0.9 + reflect_color.g/(2*level);
//			newcolor.b = b * 0.9 + reflect_color.b/(2*level);
//        }
//        else if(obj == 2 && obj_num == 0){
//			reflect_color = Ray_Tracer_recur(level, newFrom, rlx, rly, rlz, 1, obj, obj_num);
//
//			newcolor.r = r * 0.5 + reflect_color.r/(1.2*level) ;
//			newcolor.g = g * 0.5 + reflect_color.g/(1.2*level);
//			newcolor.b = b * 0.5 + reflect_color.b/(1.2*level);
//        }
        else{
			reflect_color = Ray_Tracer_recur(level, newFrom, rlx, rly, rlz, 1, obj, obj_num);

			newcolor.r = r + reflect_color.r/(1.2*level) ;
			newcolor.g = g + reflect_color.g/(1.2*level);
			newcolor.b = b + reflect_color.b/(1.2*level);
        }
    return newcolor;
}

/*******************************************************************************
  Title:	Ray_Generate
  Purpose:	This function generates the primary ray from camera to every pixel
            of the window and call recursive ray tracer

*******************************************************************************/

void Ray_Generate(){

	int	    xp, yp, obj, obj_num, shadow_flag;
	int	    texture, buf_ptr;
	float	xv, yv, dx, dy, dz, nx, ny, nz;
	float	t_min, t1, t2, ipx, ipy, ipz;
	float   rlx, rly, rlz, rrx, rry, rrz;
	float	r, g, b;
	float   u, v;
	int     i,j,k;
	rgb     newcolor;


/*  Generate a ray for each pixel in the desired image.  */

	buf_ptr = 0;
	for (xp=0; xp<xmax_pixel; xp++) {
		u = (float)xp/xmax_pixel;

		for (yp=0; yp<ymax_pixel; yp++) {
			v = (float)yp/ymax_pixel;

/*  Compute the corresponding view port coordinates.  */

			xv = VXL + xp * xinterval;
			yv = VYB + yp * yinterval;

/*  Compute the direction of the current ray from the "From" point to the
    current position on the image.  */

			dx = ax*xv*tanv2 + bx*yv*tanv2 + cx;
			dy = ay*xv*tanv2 + by*yv*tanv2 + cy;
			dz = az*xv*tanv2 + bz*yv*tanv2 + cz;


			newcolor = Ray_Tracer_recur(1, from, dx, dy, dz, 1, -1, -1);
			/* Save the computed color intensity to the image buffer. */

			texture_R[xp + xmax_pixel * yp] = newcolor.r;
			texture_G[xp + xmax_pixel * yp] = newcolor.g;
			texture_B[xp + xmax_pixel * yp] = newcolor.b;

		}
	}

/*  Write the image to the output file.  */
//
//	printf("Writing to image...\n");
//  fwrite(&xmax_pixel, sizeof(int), 1, outpfile);
//	fwrite(&ymax_pixel, sizeof(int), 1, outpfile);
//
//	fwrite(texture_R, sizeof(float), xmax_pixel*ymax_pixel, outpfile);
//	fwrite(texture_G, sizeof(float), xmax_pixel*ymax_pixel, outpfile);
//	fwrite(texture_B, sizeof(float), xmax_pixel*ymax_pixel, outpfile);

//	fclose(outpfile);
}



/* Initialize the projection matrix.  */

void myinit(void)
{
/* attributes */

      glClearColor(0.0, 0.0, 0.0, 1.0); /* white background */

/* set up viewing */
/* 512 x 512 window with origin lower left */

      glMatrixMode(GL_PROJECTION);
      glLoadIdentity();
      gluOrtho2D(0.0, 512.0, 0.0, 512.0);
      glMatrixMode(GL_MODELVIEW);
}

/* Display the ray traced image.   A more efficient method is
   to use glDrawPixels(). */

void display( void )
{
	int s, t;
	float  r, g, b;

	glClear(GL_COLOR_BUFFER_BIT);  /*clear the window */

	for(t = 0; t < ymax_pixel; t++) {
		for(s = 0; s < xmax_pixel; s++) {

			r = texture_R[s + xmax_pixel * t];
			g = texture_G[s + xmax_pixel * t];
			b = texture_B[s + xmax_pixel * t];

			glColor3f(r, g, b);
			glBegin(GL_POINTS);
               glVertex2f(s,t);
			glEnd();
		}
	 }

     glFlush(); /* clear buffers */
}


int main(int argc, char**argv)
{

/* Standard GLUT initialization */

    glutInit(&argc,argv);
    glutInitDisplayMode (GLUT_SINGLE | GLUT_RGB);   /* default, not needed */
    glutInitWindowSize(500,500);                    /* 500 x 500 pixel window */
    glutInitWindowPosition(0,0);                    /* place window top left on display */
    glutCreateWindow("Ray Trace");                  /* window title */
    glutDisplayFunc(display);                       /* display callback invoked when window opened */

	Read_Information();
	Setup_Parameters();
	Ray_Generate();

    myinit();       /* set attributes */
    glutMainLoop(); /* enter event loop */

	return(0);
}
