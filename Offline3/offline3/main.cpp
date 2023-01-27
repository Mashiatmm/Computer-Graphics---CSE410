/*
 * GLUT Shapes Demo
 *
 * Written by Nigel Stewart November 2003
 *
 * This program is test harness for the sphere, cone
 * and torus shapes in GLUT.
 *
 * Spinning wireframe and smooth shaded shapes are
 * displayed until the ESC or q key is pressed.  The
 * number of geometry stacks and slices can be adjusted
 * using the + and - keys.
 */

#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

#include <stdlib.h>
#include <sstream>
#include <fstream>
#include "1705005_classes.h"
#include "bitmap_image.hpp"


static int slices = 16;
static int stacks = 16;

extern vector<Object*> Objects;
extern vector<PointLight*> pointLights;
extern vector<SpotLight*> spotlights;


extern int level_recursion;
int pixels;

double cameraHeight;
double cameraAngle;
int drawgrid;
int drawaxes;
double angle;

int windowWidth = 500;
int windowHeight = 500;

int image_count = 11;

Point eye;
Point u, r, l;
double mov_angle;

double glob_rx = 10;
double glob_ry = 10;
double glob_rz = 10;
double glob_radius = 6  ;



void rotate_axis(Point* unit_vec1, Point* unit_vec2){
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




/* GLUT callback Handlers */

static void resize(int width, int height)
{
    const float ar = (float) width / (float) height;

    glViewport(0, 0, width, height);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glFrustum(-ar, ar, -1.0, 1.0, 2.0, 100.0);

    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity() ;
}



static void key(unsigned char key, int x, int y)
{
    switch (key)
    {
        case 27 :
        case 'q':
            exit(0);
            break;

        case '+':
            slices++;
            stacks++;
            break;

        case '-':
            if (slices>3 && stacks>3)
            {
                slices--;
                stacks--;
            }
            break;
    }

    glutPostRedisplay();
}

static void idle(void)
{
    glutPostRedisplay();
}

const GLfloat light_ambient[]  = { 0.0f, 0.0f, 0.0f, 1.0f };
const GLfloat light_diffuse[]  = { 1.0f, 1.0f, 1.0f, 1.0f };
const GLfloat light_specular[] = { 1.0f, 1.0f, 1.0f, 1.0f };
const GLfloat light_position[] = { 2.0f, 5.0f, 5.0f, 0.0f };

const GLfloat mat_ambient[]    = { 0.7f, 0.7f, 0.7f, 1.0f };
const GLfloat mat_diffuse[]    = { 0.8f, 0.8f, 0.8f, 1.0f };
const GLfloat mat_specular[]   = { 1.0f, 1.0f, 1.0f, 1.0f };
const GLfloat high_shininess[] = { 100.0f };


void loadData(){
    cout<<"hello"<<endl;

    ifstream fin;
    string line;
    // by default open mode = ios::in mode
    fin.open("E:/4-1/CSE410/Offline3/offline3/demo.txt");
    fin>>level_recursion;
    fin>>pixels;

    int count_obj;
    fin>>count_obj;

    double color[3];
    double coeff[4];
    int shine;

    for(int i = 0 ; i < count_obj ; i++){
        string obj_type;
        fin >> obj_type;
        cout << obj_type << endl;
        if(obj_type == "triangle"){

            Object *temp;

            Point vertices[3];
            for(int j = 0 ; j < 3 ; j++){
                fin >> vertices[j].x >> vertices[j].y >> vertices[j].z;
            }
            temp = new Triangle(vertices[0], vertices[1], vertices[2]);

            fin >> color[0] >> color[1] >> color[2];
            temp->setColor(color);

            fin >> coeff[0] >> coeff[1] >> coeff[2] >> coeff[3];
            temp->setCoefficients(coeff);

            fin >> shine;
            temp->setShine(shine);

            Objects.push_back(temp);

        }
        else if(obj_type == "sphere"){
            Object *temp;

            Point center;
            fin >> center.x >> center.y >> center.z;
            double radius;
            fin >> radius;
            temp = new Sphere(center, radius);

            fin >> color[0] >> color[1] >> color[2];
            temp->setColor(color);

            fin >> coeff[0] >> coeff[1] >> coeff[2] >> coeff[3];
            temp->setCoefficients(coeff);

            fin >> shine;
            temp->setShine(shine);

            Objects.push_back(temp);


        }
        else if(obj_type == "general"){

            string line;
            getline(fin, line);
            getline(fin, line);
            // cout << line << endl;
            istringstream ss(line);
            string token;
            vector<double> degree_coeff;
            while(getline(ss, token, ' ')) {
                // cout << "token  : ";
                // cout << token << '\n';
                degree_coeff.push_back( stod(token) );
            }

            Object* temp = new General(degree_coeff);
            Point ref_point;
            double length, width, height;
            double color[3];
            double coeff[4];
            int shine;
            fin >> ref_point.x >> ref_point.y >> ref_point.z >> length >> width >> height;
            fin >> color[0] >> color[1] >> color[2];
            fin >> coeff[0] >> coeff[1] >> coeff[2] >> coeff[3];
            fin >> shine;

            temp->reference_point = ref_point;
            temp->height = height;
            temp->width = width;
            temp->length = length;
            temp->setColor(color);
            temp->setCoefficients(coeff);
            temp->setShine(shine);

            Objects.push_back(temp);
        }


    }

    // POINT LIGHT
    int point_light_count = 0;
    fin >> point_light_count;
    cout << point_light_count << endl;
    for(int j = 0 ; j < point_light_count ; j++ ){
        Point pos;
        fin >> pos.x >> pos.y >> pos.z;
        double color[3];
        fin >> color[0] >> color[1] >> color[2];
        PointLight* pl = new PointLight(pos);
        pl->setColor(color);

        pointLights.push_back(pl);
    }

    // SPOT LIGHTS
    int spot_light_count = 0;
    fin >> spot_light_count;
    cout << spot_light_count << endl;
    for(int i = 0 ; i < spot_light_count ; i++ ){
        SpotLight* sl;
        Point pos;

        Point dir;
        double angle;

        fin >> pos.x >> pos.y >> pos.z ;
        fin >> color[0] >> color[1] >> color[2];
        fin >> dir.x >> dir.y >> dir.z;
        fin >> angle;

        PointLight pl(pos);
        pl.setColor(color);

        sl = new SpotLight(pl, dir, angle);
        spotlights.push_back(sl);


    }

    // Close the file
    fin.close();

    Object* floor_tile = new Floor(1000, 20); // 1000 - floorwidth, 20 - tilewidth
    for(int i = 0; i < 3 ; i++ ){
        color[i] = 1;
    }
    floor_tile->setColor(color);
    Objects.push_back(floor_tile);
}

void printObjects(){
   for (auto & element : Objects) {
        element->print_object();
    }
}

void printLights(){
    for (auto & element : pointLights) {
        element->PrintLight();
    }

    for (auto & element : spotlights) {
        element->PrintLight();
    }
}
/* Program entry point */

void Capture(){
    bitmap_image image(pixels, pixels);
    double epsilon = 0.00001;

    for(int i=0;i<pixels;i++){
        for(int j=0;j<pixels;j++){
            image.set_pixel(i, j, 0, 0, 0);
        }

    }

    int imageWidth = pixels;
    int imageHeight = pixels;

    double planeDistance = ( windowHeight / 2.0 ) / ( tan(cameraAngle * pi / ( 2.0 * 180 ) ) );
    Point topLeft;
    topLeft.x = eye.x + ( l.x * planeDistance ) - ( r.x * windowWidth / 2 ) + ( u.x * windowHeight / 2) ;
    topLeft.y = eye.y + ( l.y * planeDistance ) - ( r.y * windowWidth / 2 ) + ( u.y * windowHeight / 2) ;
    topLeft.z = eye.z + ( l.z * planeDistance ) - ( r.z * windowWidth / 2 ) + ( u.z * windowHeight / 2) ;

    double du = windowWidth * 1.0 / imageWidth;
    double dv = windowHeight * 1.0 / imageHeight;

    topLeft.x = topLeft.x + (r.x) * ( 0.5 * du ) - (u.x) * ( 0.5 * dv );
    topLeft.y = topLeft.y + (r.y) * ( 0.5 * du ) - (u.y) * ( 0.5 * dv );
    topLeft.z = topLeft.z + (r.z) * ( 0.5 * du ) - (u.z) * ( 0.5 * dv );

    int nearest;
    double t, tMin;
    for(int i = 0 ; i < imageWidth ; i++ ){
        for(int j = 0 ; j < imageHeight ; j++ ){
            Point curPixel;
            curPixel.x = topLeft.x + ( i * r.x * du ) - ( j * u.x * dv  ) ;
            curPixel.y = topLeft.y + ( i * r.y * du ) - ( j * u.y * dv  );
            curPixel.z = topLeft.z + ( i * r.z * du ) - ( j * u.z * dv  );
            // curPixel.printPoint();

            Ray* ray = new Ray(eye, curPixel);
            // ray.printRay();

            double *color_ray = new double[3];

            Object *nearest = NULL;
            int near_idx;
            double t_min = INT_MAX;
            int count = 0;

             for (auto & element : Objects) {
                // element->intersect();
                double* color_temp = new double[3];
                double t = element->intersect(ray, color_temp, 0);
                // double t = Objects[3]-> intersect(ray, color_temp, 0);
                if( ( t < t_min ) && ( t > epsilon)){
                    // cout << "Near idx : " << near_idx << ", t : " << t << endl;
                    t_min = t;
                    near_idx = count;
                }
                count += 1;
            }

            // cout << t_min << endl;

            if(t_min != INT_MAX ){
                // cout << " Near index : " << endl;
                // cout << near_idx << endl;
                nearest = Objects[near_idx];
                t_min = nearest->intersect(ray, color_ray, 1);
                // cout << color_ray[0] << " " << color_ray[1] << " " << color_ray[2] << endl;
                image.set_pixel(i, j, round( color_ray[0] * 255 ) , round( color_ray[1] * 255 ), round( color_ray[2] * 255) );

            }

        }
    }


    image.save_image("E:/4-1/CSE410/Offline3/offline3/Output_" + to_string( image_count ) + ".bmp");
    image_count += 1;

}

void drawSS()
{
    // cout << "Inside ss " << endl;
    // Objects[3]->draw();

    for (auto & element : Objects) {
        glPushMatrix();
        glTranslatef(element->reference_point.x, element->reference_point.y, element->reference_point.z);
        element->draw();
        glPopMatrix();

    }

    for( auto &pl : pointLights){
        glPushMatrix();
        pl->draw();
        glPopMatrix();
    }

    for( auto &sl : spotlights){
        glPushMatrix();
        sl->draw();
        glPopMatrix();
    }
    /*int slices = 20;
    int direction = 1;
    int height = 5;



    glTranslatef(translate_x, translate_y, radius);
    glRotatef( mov_theta, 0, 0, 1);
    glRotatef( mov_angle, 0, 1, 0);


    glTranslatef( 0, 0, -radius);
    drawWheel(slices, height);*/


}

void keyboardListener(unsigned char key, int x,int y){
	switch(key){

        case '0':
            Capture();
            break;

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
			eye.x -= 5.0*l.x;
			eye.y -= 5.0*l.y;
			eye.z -= 5.0*l.z;
			break;
		case GLUT_KEY_UP:		// up arrow key
			cameraHeight += 3.0;
			eye.x += 5.0*l.x;
			eye.y += 5.0*l.y;
			eye.z += 5.0*l.z;
			break;

		case GLUT_KEY_RIGHT:
			cameraAngle += 0.03;
			eye.x += 5.0*r.x;
			eye.y += 5.0*r.y;
			eye.z += 5.0*r.z;
			break;
		case GLUT_KEY_LEFT:
			cameraAngle -= 0.03;
			eye.x -= 5.0*r.x;
			eye.y -= 5.0*r.y;
			eye.z -= 5.0*r.z;
			break;

		case GLUT_KEY_PAGE_UP:
		    eye.x += 5.0*u.x;
			eye.y += 5.0*u.y;
			eye.z += 5.0*u.z;
			break;
		case GLUT_KEY_PAGE_DOWN:
		    eye.x -= 5.0*u.x;
			eye.y -= 5.0*u.y;
			eye.z -= 5.0*u.z;
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
	//gluLookAt(eye.x, eye.y, eye.z,	eye.x + 0,eye.y + 0,eye.z + 0,	0,1,0);
	//printf("%f %f %f\n", eye.x, eye.y, eye.z);
	gluLookAt(eye.x, eye.y, eye.z,     eye.x + l.x, eye.y + l.y, eye.z + l.z,	u.x, u.y , u.z);
	//gluLookAt(eye.x, eye.y, eye.z,    0, 0, 0,	0, 1.0, 0);


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
	cameraAngle= 80.0;
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

	eye.x = 100;
	eye.y = 100;
	eye.z = 0;
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
	gluPerspective(cameraAngle,	1,	1,	1000.0);
	//field of view in the Y (vertically)
	//aspect ratio that determines the field of view in the X direction (horizontally)
	//near distance
	//far distance
}

void freeMemory(){

    for (int i = 0; i < Objects.size(); ++i) {
        delete Objects[i]; // Calls ~object and deallocates *tmp[i]
    }
    Objects.clear();

    for (int i = 0; i < pointLights.size(); ++i) {
        delete pointLights[i]; // Calls ~object and deallocates *tmp[i]
    }
    pointLights.clear();

    for (int i = 0; i < spotlights.size(); ++i) {
        delete spotlights[i]; // Calls ~object and deallocates *tmp[i]
    }
    spotlights.clear();

    cout << "Freeing Memory " << endl;

}


int main(int argc, char *argv[])
{
    atexit(freeMemory);

    loadData();
    printObjects();
    printLights();

    glutInit(&argc,argv);
	glutInitWindowSize(windowWidth, windowHeight);
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

    return EXIT_SUCCESS;
}
