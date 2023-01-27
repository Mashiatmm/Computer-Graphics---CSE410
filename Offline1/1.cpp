#include<stdio.h>
#include<stdlib.h>
#include<math.h>

#include <windows.h>
#include <GL/glut.h>

#define pi (2*acos(0.0))

double cameraHeight;
double cameraAngle;
int drawgrid;
int drawaxes;
double angle;

struct point
{
	double x,y,z;
};

point pos;
point u, r, l;
double mov_angle;

double glob_rx = 10;
double glob_ry = 10;
double glob_rz = 10;
double glob_radius = 6  ;



void rotate_axis(point* unit_vec1, point* unit_vec2){
    double radian_angle = mov_angle * pi / 180;
    unit_vec1->x = unit_vec1->x * cos(radian_angle) + unit_vec2->x * sin(radian_angle);
    unit_vec1->y = unit_vec1->y * cos(radian_angle) + unit_vec2->y * sin(radian_angle);
    unit_vec1->z = unit_vec1->z * cos(radian_angle) + unit_vec2->z * sin(radian_angle);

    unit_vec2->x = unit_vec2->x * cos(radian_angle) - unit_vec1->x * sin(radian_angle);
    unit_vec2->y = unit_vec2->y * cos(radian_angle) - unit_vec1->y * sin(radian_angle);
    unit_vec2->z = unit_vec2->z * cos(radian_angle) - unit_vec1->z * sin(radian_angle);
}

void drawAxes()
{
	if(drawaxes==1)
	{
		glColor3f(1.0, 1.0, 1.0);
		glBegin(GL_LINES);{
			glVertex3f( 100,0,0);
			glVertex3f(-100,0,0);

			glVertex3f(0,-100,0);
			glVertex3f(0, 100,0);

			glVertex3f(0,0, 100);
			glVertex3f(0,0,-100);
		}glEnd();
	}
}


void drawGrid()
{
	int i;
	if(drawgrid==1)
	{
		glColor3f(0.6, 0.6, 0.6);	//grey
		glBegin(GL_LINES);{
			for(i=-8;i<=8;i++){

				if(i==0)
					continue;	//SKIP the MAIN axes

				//lines parallel to Y-axis
				glVertex3f(i*10, -90, 0);
				glVertex3f(i*10,  90, 0);

				//lines parallel to X-axis
				glVertex3f(-90, i*10, 0);
				glVertex3f( 90, i*10, 0);
			}
		}glEnd();
	}
}

void drawSquare(double a)
{
    glColor3f(1.0,1.0,1.0);
	glBegin(GL_QUADS);{
		glVertex3f( a, a, 0);
		glVertex3f( a,-a, 0);
		glVertex3f(-a,-a, 0);
		glVertex3f(-a, a, 0);
	}glEnd();
}

void drawSquares(double ax, double ay, double az, double radius)
{
    glColor3f(1.0, 1.0, 1.0);
	glBegin(GL_QUADS);{
		glVertex3f( ax, ay, az + radius);
		glVertex3f( ax,-ay, az + radius);
		glVertex3f(-ax,-ay, az + radius);
		glVertex3f(-ax, ay, az + radius);
	}glEnd();

	glBegin(GL_QUADS);{
		glVertex3f( ax + radius, ay, az );
		glVertex3f( ax + radius, ay, -az );
		glVertex3f( ax + radius, -ay, -az );
		glVertex3f( ax + radius, -ay, az );
	}glEnd();

	glBegin(GL_QUADS);{
		glVertex3f( ax, ay + radius, az );
		glVertex3f( ax, ay + radius, -az );
		glVertex3f( -ax, ay + radius, -az );
		glVertex3f( -ax, ay + radius, az );
	}glEnd();
}


void drawCircle(double radius,int segments)
{
    int i;
    struct point points[100];
    glColor3f(0.7,0.7,0.7);
    //generate points
    for(i=0;i<=segments;i++)
    {
        points[i].x=radius*cos(((double)i/(double)segments)*2*pi);
        points[i].y=radius*sin(((double)i/(double)segments)*2*pi);
    }
    //draw segments using generated points
    for(i=0;i<segments;i++)
    {
        glBegin(GL_LINES);
        {
			glVertex3f(points[i].x,points[i].y,0);
			glVertex3f(points[i+1].x,points[i+1].y,0);
        }
        glEnd();
    }
}

void drawOneForthCircle(double radius,int segments)
{
    int i;
    struct point points[100];
    //generate points
    for(i=0;i<=segments;i++)
    {
        points[i].x=radius*cos(((double)i/(double)segments)* ( pi / 2 ) );
        points[i].y=radius*sin(((double)i/(double)segments)* ( pi / 2 ) );
        points[i].z = 0;
    }
    //draw segments using generated points
    for(i=0;i<segments;i++)
    {
        glBegin(GL_LINES);
        {
			glVertex3f(points[i].x,points[i].y, points[i].z);
			glVertex3f(points[i+1].x,points[i+1].y, points[i+1].z);
        }
        glEnd();
    }
}

void drawCone(double radius,double height,int segments)
{
    int i;
    double shade;
    struct point points[100];
    //generate points
    for(i=0;i<=segments;i++)
    {
        points[i].x=radius*cos(((double)i/(double)segments)*2*pi);
        points[i].y=radius*sin(((double)i/(double)segments)*2*pi);
    }
    //draw triangles using generated points
    for(i=0;i<segments;i++)
    {
        //create shading effect
        if(i<segments/2)shade=2*(double)i/(double)segments;
        else shade=2*(1.0-(double)i/(double)segments);
        glColor3f(shade,shade,shade);

        glBegin(GL_TRIANGLES);
        {
            glVertex3f(0,0,height);
			glVertex3f(points[i].x,points[i].y,0);
			glVertex3f(points[i+1].x,points[i+1].y,0);
        }
        glEnd();
    }
}


void drawSphere(double radius,int slices,int stacks)
{
	struct point points[100][100];
	int i,j;
	double h,r;
	//generate points
	for(i=0;i<=stacks;i++)
	{
		h=radius*sin(((double)i/(double)stacks)*(pi/2));
		r=radius*cos(((double)i/(double)stacks)*(pi/2));
		for(j=0;j<=slices;j++)
		{
			points[i][j].x=r*cos(((double)j/(double)slices)*2*pi);
			points[i][j].y=r*sin(((double)j/(double)slices)*2*pi);
			points[i][j].z=h;
		}
	}
	//draw quads using generated points
	for(i=0;i<stacks;i++)
	{
        glColor3f((double)i/(double)stacks,(double)i/(double)stacks,(double)i/(double)stacks);
		for(j=0;j<slices;j++)
		{
			glBegin(GL_QUADS);{
			    //upper hemisphere
				glVertex3f(points[i][j].x,points[i][j].y,points[i][j].z);
				glVertex3f(points[i][j+1].x,points[i][j+1].y,points[i][j+1].z);
				glVertex3f(points[i+1][j+1].x,points[i+1][j+1].y,points[i+1][j+1].z);
				glVertex3f(points[i+1][j].x,points[i+1][j].y,points[i+1][j].z);
                //lower hemisphere
                glVertex3f(points[i][j].x,points[i][j].y,-points[i][j].z);
				glVertex3f(points[i][j+1].x,points[i][j+1].y,-points[i][j+1].z);
				glVertex3f(points[i+1][j+1].x,points[i+1][j+1].y,-points[i+1][j+1].z);
				glVertex3f(points[i+1][j].x,points[i+1][j].y,-points[i+1][j].z);
			}glEnd();


		}
	}
}


void drawOneEightSphere(double radius,int slices,int stacks, bool upper)
{
    glColor3f(1, 0, 0);
	struct point points[100][100];
	int i,j;
	double h,r;
	//generate points
	for(i=0;i<=stacks;i++)
	{
		h=radius*sin(((double)i/(double)stacks)*(pi/2));
		r=radius*cos(((double)i/(double)stacks)*(pi/2));
		if(upper == false){
            h = -h;
        }
		for(j=0;j<=slices;j++)
		{
			points[i][j].x=r*cos(((double)j/(double)slices)*(pi / 2) );
			points[i][j].y=r*sin(((double)j/(double)slices)*(pi / 2) );
			points[i][j].z=h;

		}
	}
	//draw quads using generated points
	for(i=0;i<stacks;i++)
	{
        //glColor3f((double)i/(double)stacks,(double)i/(double)stacks,(double)i/(double)stacks);
        //glColor3f(1, 0, 0);
		for(j=0;j<slices;j++)
		{
			glBegin(GL_QUADS);{
			    //upper hemisphere

                glVertex3f(points[i][j].x,points[i][j].y,points[i][j].z);
                glVertex3f(points[i][j+1].x,points[i][j+1].y,points[i][j+1].z);
                glVertex3f(points[i+1][j+1].x,points[i+1][j+1].y,points[i+1][j+1].z);
                glVertex3f(points[i+1][j].x,points[i+1][j].y,points[i+1][j].z);
			}glEnd();

			 glBegin(GL_TRIANGLES);
            {
                glVertex3f(0,0, points[i][j].z);
                glVertex3f(points[i][j].x, points[i][j].y, points[i][j].z);
                glVertex3f(points[i][j+1].x, points[i][j+1].y, points[i][j+1].z);
            }
            glEnd();
		}
	}
}


void drawCylinder(double radius,int slices, double height)
{
    glPushMatrix();
    glColor3f(0, 1, 0);
    int stacks = int( height * 10 );
    double h = ( 1.0 / stacks ) * height;
    for(int i = 0 ; i < stacks; i++){
        glTranslatef(0, 0, h);
        drawOneForthCircle(radius, slices);
    }
    glPopMatrix();


}



void drawSS()
{
    double rx = glob_rx;
    double ry = glob_ry;
    double rz = glob_rz;
    double radius = glob_radius;
    int slice = 12;
    int stacks = 20;
    bool upper = true;
    double height = rz*2;
    double square_height = 0;
    int ver_dir = 1;
    int below_dir = 2;
    int up_dir = 3;
    glColor3f(1, 0, 0);
    //drawSphere(30,12,20);
    for(int i = 0 ; i < 2 ; i++){



        if(i%2 == 1){
            upper = false;
            rz = -rz;
        }


        // Drawing Squares
        glPushMatrix();
        // up down squares
        double mov_z  = rz + radius;
        if(rz < 0){
            mov_z = mov_z  - (2*radius);
        }
        glTranslatef(0, 0, mov_z);
        drawSquare(rx);
        glTranslatef(0, 0, -mov_z);

        // Front back squares
        glPushMatrix();
        double mov_y = ry + radius;
        if( i == 1){
            mov_y = - mov_y;
        }
        glTranslatef(0, mov_y, 0);
        glRotatef(90, 1, 0, 0);
        drawSquare(rx);
        glPopMatrix();

        // Left Right squares
        double mov_x = rx + radius;
        if( i == 1){
            mov_x = - mov_x;
        }
        glTranslatef(mov_x, 0, 0);
        glRotatef(90, 0, 1, 0);
        drawSquare(rx);

        glPopMatrix();


        glPushMatrix();
        glTranslatef(rx, ry, rz);
        drawOneEightSphere(radius, slice, stacks, upper);
        if(i == 1){
           drawCylinder(radius, slice, height);

            glPushMatrix();
            glTranslatef(-rx*2, 0, 0);
            glRotatef(90, 0, 1, 0);
            drawCylinder(radius, slice, height);
            glPopMatrix();

            glTranslatef(0, 0, -rz*2);
            glRotatef(-90, 0, 1, 0);
            drawCylinder(radius, slice, height);
        }

        glPopMatrix();



        glPushMatrix();
        glTranslatef(-rx,  ry, rz);
        glRotatef(90,0,0,1);
        drawOneEightSphere(radius, slice, stacks, upper);
        if(i == 1){
            drawCylinder(radius, slice, height);

            glPushMatrix();
            glTranslatef(-rx*2, 0, 0);
            glRotatef(90, 0, 1, 0);
            drawCylinder(radius, slice, height);
            glPopMatrix();

            glTranslatef(0, 0, -rz*2);
            glRotatef(-90, 0, 1, 0);
            drawCylinder(radius, slice, height);
        }
        glPopMatrix();




        glPushMatrix();
        glTranslatef(-rx, -ry, rz);
        glRotatef(180,0,0,1);
        drawOneEightSphere(radius, slice, stacks, upper);
        if(i == 1){
            drawCylinder(radius, slice, height);

            glPushMatrix();
            glTranslatef(-rx*2, 0, 0);
            glRotatef(90, 0, 1, 0);
            drawCylinder(radius, slice, height);
            glPopMatrix();

            glTranslatef(0, 0, -rz*2);
            glRotatef(-90, 0, 1, 0);
            drawCylinder(radius, slice, height);
        }

        glPopMatrix();

        glPushMatrix();
        glTranslatef(rx,-ry, rz);
        glRotatef(270,0,0,1);
        drawOneEightSphere(radius, slice, stacks, upper);
        if(i == 1){

            drawCylinder(radius, slice, height);

            glPushMatrix();
            glTranslatef(-rx*2, 0, 0);
            glRotatef(90, 0, 1, 0);
            drawCylinder(radius, slice, height);
            glPopMatrix();

            glTranslatef(0, 0, -rz*2);
            glRotatef(-90, 0, 1, 0);
            drawCylinder(radius, slice, height);

        }
        glPopMatrix();



    }

}

void keyboardListener(unsigned char key, int x,int y){
	switch(key){

		case '1':
			//drawgrid=1-drawgrid;
            mov_angle = - mov_angle;
            rotate_axis(&l, &r);
            mov_angle = - mov_angle;
			break;
        case '2':
            rotate_axis(&l, &r);
            break;
        case '3':
            rotate_axis(&l, &u);
            break;
        case '4':
            mov_angle = - mov_angle;
            rotate_axis(&l, &u);
            mov_angle = - mov_angle;
            break;
        case '5':
            rotate_axis(&u, &r);
            break;
        case '6':
            mov_angle = - mov_angle;
            rotate_axis(&u, &r);
            mov_angle = - mov_angle;
            break;

		default:
			break;
	}
}


void specialKeyListener(int key, int x,int y){
	switch(key){
		case GLUT_KEY_DOWN:		//down arrow key
			cameraHeight -= 3.0;
			pos.x -= 5.0*l.x;
			pos.y -= 5.0*l.y;
			pos.z -= 5.0*l.z;
			break;
		case GLUT_KEY_UP:		// up arrow key
			cameraHeight += 3.0;
			pos.x += 5.0*l.x;
			pos.y += 5.0*l.y;
			pos.z += 5.0*l.z;
			break;

		case GLUT_KEY_RIGHT:
			cameraAngle += 0.03;
			pos.x += 5.0*r.x;
			pos.y += 5.0*r.y;
			pos.z += 5.0*r.z;
			break;
		case GLUT_KEY_LEFT:
			cameraAngle -= 0.03;
			pos.x -= 5.0*r.x;
			pos.y -= 5.0*r.y;
			pos.z -= 5.0*r.z;
			break;

		case GLUT_KEY_PAGE_UP:
		    pos.x += 5.0*u.x;
			pos.y += 5.0*u.y;
			pos.z += 5.0*u.z;
			break;
		case GLUT_KEY_PAGE_DOWN:
		    pos.x -= 5.0*u.x;
			pos.y -= 5.0*u.y;
			pos.z -= 5.0*u.z;
			break;

		case GLUT_KEY_INSERT:
			break;

		case GLUT_KEY_HOME:
		    if(glob_rx >= 1){
                glob_radius += 1;
                glob_rx -= 1;
                glob_ry -= 1;
                glob_rz -= 1;
		    }

			break;
		case GLUT_KEY_END:
		    if(glob_radius >= 1){
                 glob_radius -= 1;
                glob_rx += 1;
                glob_ry += 1;
                glob_rz += 1;
		    }

			break;

		default:
			break;
	}
}


void mouseListener(int button, int state, int x, int y){	//x, y is the x-y of the screen (2D)
	switch(button){
		case GLUT_LEFT_BUTTON:
			if(state == GLUT_DOWN){		// 2 times?? in ONE click? -- solution is checking DOWN or UP
				drawaxes=1-drawaxes;
			}
			break;

		case GLUT_RIGHT_BUTTON:
			//........
			break;

		case GLUT_MIDDLE_BUTTON:
			//........
			break;

		default:
			break;
	}
}



void display(){

	//clear the display
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glClearColor(0,0,0,0);	//color black
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	/********************
	/ set-up camera here
	********************/
	//load the correct matrix -- MODEL-VIEW matrix
	glMatrixMode(GL_MODELVIEW);

	//initialize the matrix
	glLoadIdentity();

	//now give three info
	//1. where is the camera (viewer)?
	//2. where is the camera looking?
	//3. Which direction is the camera's UP direction?

	//gluLookAt(100,100,100,	0,0,0,	0,0,1);
	//gluLookAt(200*cos(cameraAngle), 200*sin(cameraAngle), cameraHeight,		0,0,0,		0,0,1);
	//gluLookAt(pos.x, pos.y, pos.z,	pos.x + 0,pos.y + 0,pos.z + 0,	0,1,0);
	//printf("%f %f %f\n", pos.x, pos.y, pos.z);
	gluLookAt(pos.x, pos.y, pos.z,     pos.x + l.x, pos.y + l.y, pos.z + l.z,	u.x, u.y , u.z);
	//gluLookAt(pos.x, pos.y, pos.z,    0, 0, 0,	0, 1.0, 0);


	//again select MODEL-VIEW
	glMatrixMode(GL_MODELVIEW);


	/****************************
	/ Add your objects from here
	****************************/
	//add objects

	drawAxes();
	drawGrid();

    //glColor3f(1,0,0);
    //drawSquare(10);

    drawSS();

    //drawCircle(30,24);

    //drawCone(20,50,24);

	//drawSphere(30,24,20);




	//ADD this line in the end --- if you use double buffer (i.e. GL_DOUBLE)
	glutSwapBuffers();
}


void animate(){
	angle+=0.05;
	//codes for any changes in Models, Camera
	glutPostRedisplay();
}

void init(){
	//codes for initialization
	drawgrid=0;
	drawaxes=1;
	cameraHeight=150.0;
	cameraAngle=1.0;
	angle=0;

	u.x = 0;
	u.y = 0;
	u.z = 1;
	r.x = - 0.70710678118;
	r.y = 0.70710678118;
	r.z = 0;
	l.x = - 0.70710678118;
	l.y = - 0.70710678118;
	l.z = 0;

	pos.x = 50;
	pos.y = 50;
	pos.z = 0;
	mov_angle = 5.0;






	//clear the screen
	glClearColor(0,0,0,0);

	/************************
	/ set-up projection here
	************************/
	//load the PROJECTION matrix
	glMatrixMode(GL_PROJECTION);

	//initialize the matrix
	glLoadIdentity();

	//give PERSPECTIVE parameters
	gluPerspective(80,	1,	1,	1000.0);
	//field of view in the Y (vertically)
	//aspect ratio that determines the field of view in the X direction (horizontally)
	//near distance
	//far distance
}

int main(int argc, char **argv){
	glutInit(&argc,argv);
	glutInitWindowSize(500, 500);
	glutInitWindowPosition(0, 0);
	glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGB);	//Depth, Double buffer, RGB color

	glutCreateWindow("My OpenGL Program");

	init();

	glEnable(GL_DEPTH_TEST);	//enable Depth Testing

	glutDisplayFunc(display);	//display callback function
	glutIdleFunc(animate);		//what you want to do in the idle time (when no drawing is occuring)

	glutKeyboardFunc(keyboardListener);
	glutSpecialFunc(specialKeyListener);
	glutMouseFunc(mouseListener);

	glutMainLoop();		//The main loop of OpenGL

	return 0;
}
