#include <iostream>
#include <vector>
#include<math.h>

#define pi (2*acos(0.0))

using namespace std;

class Point{
public:
    double x;
    double y;
    double z;

    Point(){}

    Point(double x, double y, double z){
        this->x = x;
        this->y = y;
        this->z = z;
    }

    double length(){
        return sqrt( pow(x, 2) + pow(y, 2) + pow(z, 2) );
    }

    double normalize(){
        double len = length();
        this->x = x * 1.0 / len;
        this->y = y * 1.0 / len;
        this->z = z * 1.0 / len;
    }

    void printPoint(){
        cout << " Point : " << x << " " << y << " " << z << endl;
        return;
    }
};

class Ray{
public:
    Point start;
    Point dir;

    Ray(Point eye, Point cur_pixel){
        start = eye;
        dir.x = cur_pixel.x - eye.x;
        dir.y = cur_pixel.y - eye.y;
        dir.z = cur_pixel.z - eye.z;

        double length = sqrt( ( dir.x * dir.x ) + ( dir.y * dir.y ) + ( dir.z * dir.z ) );
        dir.x = dir.x / length;
        dir.y = dir.y / length;
        dir.z = dir.z / length;
    }

    Ray(){}

    void setDir(Point dir){
        double length = sqrt( ( dir.x * dir.x ) + ( dir.y * dir.y ) + ( dir.z * dir.z ) );
        this->dir.x = dir.x / length;
        this->dir.y = dir.y / length;
        this->dir.z = dir.z / length;
    }

    void setStart(Point start){
        this->start = start;
    }

    void printRay(){
        cout << "Ray start  " << endl;
        start.printPoint() ;
        cout << "Ray dir  " << endl;
        dir.printPoint() ;
    }
};

class PointLight{
public:
    Point light_pos;
    double color[3];

    PointLight(Point pos){
        this->light_pos = pos;
    }

    void setColor(double c[]){
        this->color[0] = c[0];
        this->color[1] = c[1];
        this->color[2] = c[2];
    }

    void draw(){
        glColor3f(color[0], color[1], color[2]);
        glBegin(GL_QUADS);{
            //upper hemisphere
            glVertex3f(light_pos.x, light_pos.y, light_pos.z);
            glVertex3f(light_pos.x + 5, light_pos.y, light_pos.z);
            glVertex3f(light_pos.x + 5, light_pos.y + 5, light_pos.z);
            glVertex3f(light_pos.x, light_pos.y + 5, light_pos.z);
        }glEnd();

    }

    void PrintLight(){
        cout << "PointLight " << endl;
        light_pos.printPoint();
        cout << "Color : " << color[0] << " " << color[1] << " " << color[2] << endl;
    }
};

class SpotLight{
public:
    PointLight* pl;
    Point light_direction;
    double cutoff_angle;

    SpotLight(PointLight pl, Point light_dir, double angle){
        this-> pl = new PointLight(pl.light_pos);
        this->pl->setColor(pl.color);
        this->light_direction = light_dir;
        this->cutoff_angle = angle;
    }

    void draw(){
        pl->draw();

    }

    void PrintLight(){
        cout << "Spotlight " << endl;
        pl->PrintLight();
        light_direction.printPoint();
        cout << "Cutoff Angle : " << cutoff_angle << endl;
    }
};


class Object{
public:
    Point reference_point;
    double height, width, length;
    double color[3];
    double coEfficients[4];
    double shine;

    Object(){
    }

    virtual void draw(){}

    virtual double intersect(Ray *r, double* color, int level){
        return INT_MAX;
    }

    void setLightingColor(double* c){
        c[0] = color[0] * coEfficients[0];
        c[1] = color[1] * coEfficients[0];
        c[2] = color[2] * coEfficients[0];
    }

    void setColor(double c[]){
        this->color[0] = c[0];
        this->color[1] = c[1];
        this->color[2] = c[2];
    }

    void setShine(double shine){
        this->shine = shine;
    }

    void setCoefficients(double coeff[]){
        this->coEfficients[0] = coeff[0];
        this->coEfficients[1] = coeff[1];
        this->coEfficients[2] = coeff[2];
        this->coEfficients[3] = coeff[3];
    }

    double dot_product(Point first, Point second){
        double result = ( first.x * second.x ) + ( first.y * second.y ) + ( first.z * second.z );
        return result;
    }

    void calculateReflection(double* colorRay, double* color_in){
        color_in[0] += colorRay[0] * coEfficients[3];
        color_in[1] += colorRay[1] * coEfficients[3];
        color_in[2] += colorRay[1] * coEfficients[3];

        if( color_in[0] > 1.0 ) color_in[0] = 1;
        if( color_in[1] > 1.0 ) color_in[1] = 1;
        if( color_in[2] > 1.0 ) color_in[2] = 1;
    }

    void calculateSpecularDiffuse(Point normal, Ray* ray_point_light, double* color_in, Point intersection_point, Ray* r, PointLight* pl){

        // double diffuse_coeff = this->coEfficients[1] * cos_theta;
        double cos_theta = dot_product(normal, ray_point_light->dir);
        double lambert_value = -1.0 * cos_theta; // both unit vector so no need to divide by length to get costheta
        Point reflectedRayDir( - (2.0 * cos_theta * normal.x ) + ray_point_light->dir.x, - (2.0 * cos_theta * normal.y ) + ray_point_light->dir.y, - (2.0 * cos_theta * normal.z ) + ray_point_light->dir.z);
        reflectedRayDir.normalize();
        Ray reflectedRay;
        reflectedRay.setStart( intersection_point );
        reflectedRay.setDir(reflectedRayDir);
        double phong_value = -1.0 * dot_product(reflectedRay.dir, r->dir); // both unit vector so no need to divide by length for getting costheta

        double const_diffuse = coEfficients[1] * max(lambert_value, 0.0);
        double const_specular = coEfficients[2] * pow( max(phong_value, 0.0), shine);

        color_in[0] += color[0] * pl->color[0] * const_diffuse;
        color_in[1] += color[1] * pl->color[1] * const_diffuse;
        color_in[2] += color[2] * pl->color[2] * const_diffuse;

        color_in[0] += color[0] * pl->color[0] * const_specular;
        color_in[1] += color[1] * pl->color[1] * const_specular;
        color_in[2] += color[2] * pl->color[2] * const_specular;


        if( color_in[0] > 1.0 ) color_in[0] = 1;
        if( color_in[1] > 1.0 ) color_in[1] = 1;
        if( color_in[2] > 1.0 ) color_in[2] = 1;

    }


    virtual void print_object(){
        cout << "Color : " << color[0] << " " << color[1] << " " << color[2] << endl;
        cout << "Coefficients: " << coEfficients[0] << " " << coEfficients[1] << " " << coEfficients[2] << " " << coEfficients[3] << endl;
        cout << "Shine : " << shine << endl;
    }

};


vector<Object*> Objects;
vector<PointLight*> pointLights;
vector<SpotLight*> spotlights;
int level_recursion;



class Sphere : public Object{
public:
    Sphere(Point center, double radius){
        reference_point = center;
        length = radius;
    }

    void draw(){
        // printf("Draw circle\n");
        Point points[100][100];
        int i,j;
        int slices = 12;
        int stacks = 20;
        double h,r;
        double radius = length;
       // cout << "radius : " << radius << endl;
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
            glColor3f(color[0], color[1], color[2]);
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
        // cout<<"drawing is over" << endl;
    }

    double intersect(Ray *r, double* color_in, int level){
        double epsilon = 0.00001;
        Point r_start_translate( r->start.x - reference_point.x, r->start.y - reference_point.y, r->start.z - reference_point.z);


        double ld = - ( r->dir.x * r_start_translate.x ) - ( r->dir.y * r_start_translate.y ) - ( r->dir.z * r_start_translate.z );
        if(ld < 0 ){
            return INT_MAX;
        }

        double r_start_len = ( r_start_translate.x * r_start_translate.x ) + ( r_start_translate.y * r_start_translate.y ) + ( r_start_translate.z * r_start_translate.z );
        double d_squared = ( r_start_len ) - ( ld * ld );
        if(d_squared > ( length * length ) ){
            return INT_MAX;
        }

        double l_prime = sqrt( ( length*length ) - d_squared );
        double t1 = ld + l_prime;
        double t2 = ld - l_prime;



        // if(level == 1) cout << "intersecting sphere on  " << min(t1, t2) << endl;

        double t = NULL;

        if(r_start_len < (length * length) ){
            t = t1;
        }
        else if(r_start_len > ( length * length ) ){
            t = t2;
        }
        else if(t == NULL ){
            t = min(t1, t2);
        }
        // t = min(t1, t2);
        if( level == 0 ){
            return t;
        }

        // level != 0
        // cout << "level " << level << endl;
        // cout << pointLights.size() << endl;
        setLightingColor(color_in);

        Point intersect_point( r->start.x + ( t * r->dir.x ), r->start.y + ( t * r->dir.y ), r->start.z + ( t * r->dir.z ) );
        Point normal( intersect_point.x - reference_point.x , intersect_point.y - reference_point.y, intersect_point.z - reference_point.z );
        normal.normalize();

        // intersect_point.printPoint() ;

        vector<PointLight*> all_pl = pointLights;


        for( auto &sl : spotlights ){
            PointLight* pl = sl->pl;
            // cout << "hello " << endl;
            Ray* ray_point_light = new Ray(pl->light_pos, intersect_point);
            // ray_point_light->printRay();

            // check if crosses cutoff angle
            double alpha = dot_product(ray_point_light->dir, sl->light_direction);
            double light_dir_len = sl->light_direction.length();
            double ray_len = ray_point_light->dir.length();
            alpha = acos( alpha * 1.0 / (light_dir_len * ray_len) ) * 180.0 / pi;
            // cout << "angle : " << alpha << endl;
            if(alpha > sl->cutoff_angle ){
                // cout << "not spotlight added " << endl;
                continue;
            }
            else{
                all_pl.push_back(pl);
            }
        }

        for( auto &pl : all_pl ){

            // cout << "hello " << endl;
            Ray* ray_point_light = new Ray(pl->light_pos, intersect_point);
            // ray_point_light->printRay();
            double min_t_pl = INT_MAX;
            for( auto &obj : Objects ){
                double* dummyColor = new double[3];
                // cout << "here" << endl;
                double t_pl = obj->intersect(ray_point_light, dummyColor, 0);
                // cout << "t_pl : " << t_pl << endl;
                if( ( t_pl > epsilon ) && ( t_pl < min_t_pl ) ){
                    min_t_pl = t_pl;
                }
            }
            if(min_t_pl != INT_MAX ){
                // cout << "inside here "<< endl;
                // cout << "Min t pl : " << min_t_pl << endl;
                Point shadow_intersect_point( ray_point_light->start.x + ( min_t_pl * ray_point_light->dir.x ), ray_point_light->start.y + ( min_t_pl * ray_point_light->dir.y ), ray_point_light->start.z + ( min_t_pl * ray_point_light->dir.z ) );

                // check if the intersection point is in shadow
                double dis_shadow_ip = sqrt( pow(shadow_intersect_point.x - ray_point_light->start.x, 2) + pow(shadow_intersect_point.y - ray_point_light->start.y, 2) + pow(shadow_intersect_point.z - ray_point_light->start.z, 2) );
                double dis_ip = sqrt( pow(intersect_point.x - ray_point_light->start.x, 2) + pow(intersect_point.y - ray_point_light->start.y, 2) + pow(intersect_point.z - ray_point_light->start.z, 2) );


                // cout << "dis shadow  : " << dis_shadow_ip << endl;
                // cout << "dis ip : " << dis_ip - epsilon << endl;
                if( dis_shadow_ip < dis_ip - epsilon ){
                    continue;
                }

                // not in shadow, calculate specular_diffuse components
                // cout << "hello" << endl;
                calculateSpecularDiffuse(normal, ray_point_light, color_in, intersect_point, r, pl);

            }


        }

        if(level >= level_recursion ) return t;

        double dot_ray_n = dot_product(normal, r->dir);
        Point reflectedRayDir( - (2.0 * dot_ray_n * normal.x ) + r->dir.x, - (2.0 * dot_ray_n * normal.y ) + r->dir.y, - (2.0 * dot_ray_n * normal.z ) + r->dir.z);
        reflectedRayDir.normalize();

        Point reflect_initial(intersect_point.x + reflectedRayDir.x , intersect_point.y + reflectedRayDir.y, intersect_point.z + reflectedRayDir.z);
        Ray* reflectedRay = new Ray();

        reflectedRay->setStart(reflect_initial);
        reflectedRay->setDir(reflectedRayDir);

        // ray.printRay();

        double *color_ray = new double[3];

        Object *nearest = NULL;
        int near_idx;
        double t_min = INT_MAX;
        int count = 0;

         for (auto & element : Objects) {
            // element->intersect();
            double* color_temp = new double[3];
            double t = element->intersect(reflectedRay, color_temp, 0);
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
            // cout << near_idx << endl;
            nearest = Objects[near_idx];
            t_min = nearest->intersect(reflectedRay, color_ray, level + 1);
            // cout << color_ray[0] << " " << color_ray[1] << " " << color_ray[2] << endl;
        }

        calculateReflection(color_ray, color_in);

        return t;


    }

    void print_object(){
        cout << "Sphere " << endl;
        cout << "Radius : " << length << endl;
        cout << "Center : " << reference_point.x << reference_point.y << reference_point.z <<endl;
        Object::print_object();

    }
};


class Triangle : public Object{
public:
    Point first, second, third;

    Triangle(Point first, Point second, Point third){
        this->first =  first;
        this->second = second;
        this->third = third;
    }

    void draw(){
        // printf("Draw Triangle\n");

        glColor3f(color[0], color[1], color[2]);
         glBegin(GL_TRIANGLES);
        {
            glVertex3f(first.x, first.y, first.z);
            glVertex3f(second.x, second.y, second.z);
            glVertex3f(third.x, third.y, third.z);
        }
        glEnd();

    }

    void print_object(){
        cout << "First : " << first.x << first.y << first.z <<endl;
        cout << "Second : " << second.x << second.y << second.z <<endl;
        cout << "Third : " << third.x << third.y << third.z <<endl;
        Object::print_object();

    }

    double intersect(Ray *r, double* color_in, int level){

        double epsilon = 0.00001;

        setLightingColor(color_in);

        double a1 = first.x - second.x;
        double a2 = first.y - second.y;
        double a3 = first.z - second.z;

        double b1 = third.x - second.x;
        double b2 = third.y - second.y;
        double b3 = third.z - second.z;

        double c1 = - r->dir.x;
        double c2 = - r->dir.y;
        double c3 = - r->dir.z;

        double d1 = r->start.x - second.x ;
        double d2 = r->start.y - second.y ;
        double d3 = r->start.z - second.z ;

        double D = a1 * (b2*c3 - c2*b3) + b1 * (c2*a3 - c3*a2) + c1 * (a2*b3 - a3*b2);
        double D1 = d1 * (b2*c3 - c2*b3) + b1 * (c2*d3 - c3*d2) + c1 * (d2*b3 - d3*b2);
        double D2 = a1 * (d2*c3 - c2*d3) + d1 * (c2*a3 - c3*a2) + c1 * (a2*d3 - a3*d2);
        double D3 = a1 * (b2*d3 - d2*b3) + b1 * (d2*a3 - d3*a2) + d1 * (a2*b3 - a3*b2);

        if( D != 0 ){
            double k1 = D1 / D;
            double k2 = D2 / D;
            double t = D3 / D;
            // cout << (k1 + k2 ) << endl;
            if((k1 > 0 ) && (k2 > 0) && ( k1 + k2 <= 1 ) ){
                if(level == 0 ) return t;

                // level != 0
                Point intersect_point(r->start.x + ( t * r->dir.x ), r->start.y + ( t * r->dir.y ), r->start.z + ( t * r->dir.z ));
                Point normal;
                normal.x = ( a2 * b3 ) - ( b2 * a3 );
                normal.y = ( a3 * b1 ) - ( a1 * b3 );
                normal.z = ( a1 * b2 ) - ( a2 * b1 );
                normal.normalize();

                vector<PointLight*> all_pl = pointLights;


                for( auto &sl : spotlights ){
                    PointLight* pl = sl->pl;
                    // cout << "hello " << endl;
                    Ray* ray_point_light = new Ray(pl->light_pos, intersect_point);
                    // ray_point_light->printRay();

                    // check if crosses cutoff angle
                    double alpha = dot_product(ray_point_light->dir, sl->light_direction);
                    double light_dir_len = sl->light_direction.length();
                    double ray_len = ray_point_light->dir.length();
                    alpha = acos( alpha * 1.0 / (light_dir_len * ray_len) ) * 180.0 / pi;
                    // cout << "angle : " << alpha << endl;
                    if(alpha > sl->cutoff_angle ){
                        // cout << "not spotlight added " << endl;
                        continue;
                    }
                    else{
                        all_pl.push_back(pl);
                    }
                }

                for( auto &pl : all_pl ){
                    // cout << "hello " << endl;
                    Ray* ray_point_light = new Ray(pl->light_pos, intersect_point);
                    // ray_point_light->printRay();
                    double min_t_pl = INT_MAX;
                    for( auto &obj : Objects ){
                        double* dummyColor = new double[3];
                        // cout << "here" << endl;
                        double t_pl = obj->intersect(ray_point_light, dummyColor, 0);
                        // cout << "t_pl : " << t_pl << endl;
                        if( ( t_pl > epsilon ) && ( t_pl < min_t_pl ) ){
                            min_t_pl = t_pl;
                        }
                    }
                    if(min_t_pl != INT_MAX ){
                        // cout << "inside here "<< endl;
                        // cout << "Min t pl : " << min_t_pl << endl;
                        Point shadow_intersect_point( ray_point_light->start.x + ( min_t_pl * ray_point_light->dir.x ), ray_point_light->start.y + ( min_t_pl * ray_point_light->dir.y ), ray_point_light->start.z + ( min_t_pl * ray_point_light->dir.z ) );

                        // check if the intersection point is in shadow
                        double dis_shadow_ip = sqrt( pow(shadow_intersect_point.x - ray_point_light->start.x, 2) + pow(shadow_intersect_point.y - ray_point_light->start.y, 2) + pow(shadow_intersect_point.z - ray_point_light->start.z, 2) );
                        double dis_ip = sqrt( pow(intersect_point.x - ray_point_light->start.x, 2) + pow(intersect_point.y - ray_point_light->start.y, 2) + pow(intersect_point.z - ray_point_light->start.z, 2) );


                        // cout << "dis shadow  : " << dis_shadow_ip << endl;
                        // cout << "dis ip : " << dis_ip - epsilon << endl;
                        if( dis_shadow_ip < dis_ip - epsilon ){
                            continue;
                        }

                        // not in shadow, calculate specular_diffuse components
                        // cout << "hello" << endl;
                        calculateSpecularDiffuse(normal, ray_point_light, color_in, intersect_point, r, pl);

                    }

                }

                if(level >= level_recursion ) return t;

                double dot_ray_n = dot_product(normal, r->dir);
                Point reflectedRayDir( - (2.0 * dot_ray_n * normal.x ) + r->dir.x, - (2.0 * dot_ray_n * normal.y ) + r->dir.y, - (2.0 * dot_ray_n * normal.z ) + r->dir.z);
                reflectedRayDir.normalize();

                Point reflect_initial(intersect_point.x + reflectedRayDir.x , intersect_point.y + reflectedRayDir.y, intersect_point.z + reflectedRayDir.z);
                Ray* reflectedRay = new Ray();

                reflectedRay->setStart(reflect_initial);
                reflectedRay->setDir(reflectedRayDir);

                // ray.printRay();

                double *color_ray = new double[3];

                Object *nearest = NULL;
                int near_idx;
                double t_min = INT_MAX;
                int count = 0;

                 for (auto & element : Objects) {
                    // element->intersect();
                    double* color_temp = new double[3];
                    double t = element->intersect(reflectedRay, color_temp, 0);
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
                    // cout << near_idx << endl;
                    nearest = Objects[near_idx];
                    t_min = nearest->intersect(reflectedRay, color_ray, level + 1);
                    // cout << color_ray[0] << " " << color_ray[1] << " " << color_ray[2] << endl;
                }

                calculateReflection(color_ray, color_in);

                return t;

            }
        }

        return INT_MAX;

    }

};


class General : public Object{
public:
    vector<double> degree_coeff;

    General( vector<double> coeff){
        this->degree_coeff = coeff;
        /*cout<<"length : " << degree_coeff.size() << endl;
        for(double element : degree_coeff){
            cout<<element << " ";
        }
        cout<<endl;*/
    }

    void draw(){
    }

    double intersect(Ray *r, double* color_in, int level){

        double epsilon = 0.00001;

        setLightingColor(color_in);


        Point r_start(r->start.x - reference_point.x, r->start.y - reference_point.y, r->start.z - reference_point.z);



        /*r_start.x =  r->start.x ;
        r_start.y =  r->start.y ;
        r_start.z =  r->start.z ;*/

        double a = ( degree_coeff[0] * pow( r->dir.x , 2 ) ) + ( degree_coeff[1] * pow(r->dir.y, 2) ) + (degree_coeff[2] * pow(r->dir.z, 2) ) \
                            + ( degree_coeff[3] * r->dir.x * r->dir.y ) + ( degree_coeff[4] * r->dir.y * r->dir.z ) + ( degree_coeff[5] * r->dir.x * r->dir.z );
        double b = ( 2 * degree_coeff[0] * r_start.x * r->dir.x ) + ( 2* degree_coeff[1] * r_start.y * r->dir.y ) + ( 2* degree_coeff[2] * r_start.z * r->dir.z ) \
                    + ( degree_coeff[3] * ( ( r->dir.x * r_start.y ) + ( r_start.x * r->dir.y ) ) ) + ( degree_coeff[4] * ( ( r->dir.y * r_start.z ) + ( r_start.y * r->dir.z ) ) ) \
                    + ( degree_coeff[5] * ( ( r->dir.z * r_start.x ) + ( r_start.z * r->dir.x ) ) ) + ( degree_coeff[6] * r->dir.x ) + ( degree_coeff[7] * r->dir.y ) + ( degree_coeff[8] * r->dir.z );

        double c = ( degree_coeff[0] * pow( r_start.x , 2 ) ) + ( degree_coeff[1] * pow(r_start.y, 2) ) + (degree_coeff[2] * pow(r_start.z, 2) ) \
                        + ( degree_coeff[3] * r_start.x * r_start.y ) + ( degree_coeff[4] * r_start.y * r_start.z ) + ( degree_coeff[5] * r_start.x * r_start.z ) \
                        + ( degree_coeff[6] * r_start.x ) + ( degree_coeff[7] * r_start.y ) + ( degree_coeff[8] * r_start.z ) + degree_coeff[9];

        double root_squared = ( b * b ) - ( 4 * a * c );
        if(root_squared < 0 ){
            return INT_MAX;
        }
        root_squared = sqrt( root_squared );
        double t1 = ( - b + root_squared ) / ( 2 * a );
        double t2 = ( - b - root_squared ) / ( 2 * a );

        double t = INT_MAX;
        bool not_valid = false;
        Point min_intersect;
        if( t1 > 0 ){
            Point intersect_point(r->start.x + ( t1 * r->dir.x ), r->start.y + ( t1 * r->dir.y ), r->start.z + ( t1 * r->dir.z ));


            if(length != 0 ){
                if( ( intersect_point.x < reference_point.x ) || ( intersect_point.x > ( reference_point.x + length )) ){
                    not_valid = true;
                }
            }


            if(width != 0 ){
                if( ( intersect_point.y < reference_point.y ) || ( intersect_point.y > ( reference_point.y + width )) ){
                    not_valid = true;
                }
            }

            if(height != 0 ){
                if( ( intersect_point.z < reference_point.z ) || ( intersect_point.z > ( reference_point.z + height )) ){
                    not_valid = true;
                }
            }

            if( not_valid == false ){
                min_intersect = intersect_point;
                t = min(t, t1);
            }

        }
        not_valid = false;

         if( t2 > 0 ){
            Point intersect_point( r->start.x + ( t2 * r->dir.x ), r->start.y + ( t2 * r->dir.y ), r->start.z + ( t2 * r->dir.z ) );


            if(length != 0 ){
                if( ( intersect_point.x < reference_point.x ) || ( intersect_point.x > ( reference_point.x + length )) ){
                    not_valid = true;
                }
            }


            if(width != 0 ){
                if( ( intersect_point.y < reference_point.y ) || ( intersect_point.y > ( reference_point.y + width )) ){
                    not_valid = true;
                }
            }

            if(height != 0 ){
                if( ( intersect_point.z < reference_point.z ) || ( intersect_point.z > ( reference_point.z + height )) ){
                    not_valid = true;
                }
            }

            if( not_valid == false ){
                if( t != NULL ){
                    t = min(t, t2);
                    if( t == t2 ){
                        min_intersect = intersect_point;
                    }
                }
                else{
                    min_intersect = intersect_point;
                    t = t2;
                }
            }

        }


       if(level == 0 ) return t;

       Point intersect_point = min_intersect;

       Point normal;
       normal.x = ( 2 * degree_coeff[0] * min_intersect.x ) + ( degree_coeff[3] * min_intersect.y ) + ( degree_coeff[5] * min_intersect.z ) + degree_coeff[6];
       normal.y = ( 2 * degree_coeff[1] * min_intersect.y ) + ( degree_coeff[3] * min_intersect.x ) + ( degree_coeff[4] * min_intersect.y ) + degree_coeff[7];
       normal.z = ( 2 * degree_coeff[2] * min_intersect.z ) + ( degree_coeff[4] * min_intersect.y ) + ( degree_coeff[5] * min_intersect.x ) + degree_coeff[8];
       normal.normalize();

       vector<PointLight*> all_pl = pointLights;


        for( auto &sl : spotlights ){
            PointLight* pl = sl->pl;
            // cout << "hello " << endl;
            Ray* ray_point_light = new Ray(pl->light_pos, intersect_point);
            // ray_point_light->printRay();

            // check if crosses cutoff angle
            double alpha = dot_product(ray_point_light->dir, sl->light_direction);
            double light_dir_len = sl->light_direction.length();
            double ray_len = ray_point_light->dir.length();
            alpha = acos( alpha * 1.0 / (light_dir_len * ray_len) ) * 180.0 / pi;
            // cout << "angle : " << alpha << endl;
            if(alpha > sl->cutoff_angle ){
                // cout << "not spotlight added " << endl;
                continue;
            }
            else{
                all_pl.push_back(pl);
            }
        }

        for( auto &pl : all_pl ){
            // cout << "hello " << endl;
            Ray* ray_point_light = new Ray(pl->light_pos, intersect_point);
            // ray_point_light->printRay();
            double min_t_pl = INT_MAX;
            for( auto &obj : Objects ){
                double* dummyColor = new double[3];
                // cout << "here" << endl;
                double t_pl = obj->intersect(ray_point_light, dummyColor, 0);
                // cout << "t_pl : " << t_pl << endl;
                if( ( t_pl > epsilon ) && ( t_pl < min_t_pl ) ){
                    min_t_pl = t_pl;
                }
            }
            if(min_t_pl != INT_MAX ){
                // cout << "inside here "<< endl;
                // cout << "Min t pl : " << min_t_pl << endl;
                Point shadow_intersect_point( ray_point_light->start.x + ( min_t_pl * ray_point_light->dir.x ), ray_point_light->start.y + ( min_t_pl * ray_point_light->dir.y ), ray_point_light->start.z + ( min_t_pl * ray_point_light->dir.z ) );

                // check if the intersection point is in shadow
                double dis_shadow_ip = sqrt( pow(shadow_intersect_point.x - ray_point_light->start.x, 2) + pow(shadow_intersect_point.y - ray_point_light->start.y, 2) + pow(shadow_intersect_point.z - ray_point_light->start.z, 2) );
                double dis_ip = sqrt( pow(intersect_point.x - ray_point_light->start.x, 2) + pow(intersect_point.y - ray_point_light->start.y, 2) + pow(intersect_point.z - ray_point_light->start.z, 2) );


                // cout << "dis shadow  : " << dis_shadow_ip << endl;
                // cout << "dis ip : " << dis_ip - epsilon << endl;
                if( dis_shadow_ip < dis_ip - epsilon ){
                    continue;
                }

                // not in shadow, calculate specular_diffuse components
                // cout << "hello" << endl;
                calculateSpecularDiffuse(normal, ray_point_light, color_in, intersect_point, r, pl);

            }

        }

        if(level >= level_recursion ) return t;

        double dot_ray_n = dot_product(normal, r->dir);
        Point reflectedRayDir( - (2.0 * dot_ray_n * normal.x ) + r->dir.x, - (2.0 * dot_ray_n * normal.y ) + r->dir.y, - (2.0 * dot_ray_n * normal.z ) + r->dir.z);
        reflectedRayDir.normalize();

        Point reflect_initial(intersect_point.x + reflectedRayDir.x , intersect_point.y + reflectedRayDir.y, intersect_point.z + reflectedRayDir.z);
        Ray* reflectedRay = new Ray();

        reflectedRay->setStart(reflect_initial);
        reflectedRay->setDir(reflectedRayDir);

        // ray.printRay();

        double *color_ray = new double[3];

        Object *nearest = NULL;
        int near_idx;
        double t_min = INT_MAX;
        int count = 0;

         for (auto & element : Objects) {
            // element->intersect();
            double* color_temp = new double[3];
            double t = element->intersect(reflectedRay, color_temp, 0);
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
            // cout << near_idx << endl;
            nearest = Objects[near_idx];
            t_min = nearest->intersect(reflectedRay, color_ray, level + 1);
            // cout << color_ray[0] << " " << color_ray[1] << " " << color_ray[2] << endl;
        }

        calculateReflection(color_ray, color_in);

       return t;

    }



};




class Floor : public Object{
public:
    double floorwidth;

    Floor(double floorwidth, double tilewidth){
        this->floorwidth = floorwidth;
        reference_point.x = - floorwidth / 2.0;
        reference_point.y = - floorwidth / 2.0;
        reference_point.z = 0;
        length = tilewidth;
        setCoefficients(new double[4]{0.2, 0.2, 0.4, 0.2});
        setShine(0.5);
        setColor(new double[3]{1, 1, 1});
    }

    void print_object(){
        cout << "Floor : " << endl;
        cout << "Floorwidth: " << floorwidth << endl;
        cout << "tilewidth : " << length << endl;
        Object::print_object();
    }

    void draw(){
        int i;

        Point cur_point(0, 0, 0);

        double row_color[3] = {1, 1, 1};
        // cout << "y : " << cur_point.y << " " << reference_point.y << endl;
        while( cur_point.y <= floorwidth ){

            double col_color[3] = {row_color[0], row_color[1], row_color[2]};
            while(cur_point.x <= floorwidth){
                // cout << "cur point y : " << cur_point.y << endl;
                glColor3f(col_color[0], col_color[1], col_color[2]);
                glBegin(GL_QUADS);{

                glVertex3f(cur_point.x, cur_point.y, 0);
                glVertex3f(cur_point.x + length, cur_point.y, 0);
                glVertex3f(cur_point.x + length, cur_point.y + length, 0);
                glVertex3f(cur_point.x, cur_point.y + length, 0);

                }glEnd();

                cur_point.x = cur_point.x + length;

                col_color[0] = 1 - col_color[0];
                col_color[1] = 1 - col_color[1];
                col_color[2] = 1 - col_color[2];

            }

            cur_point.y = cur_point.y + length;
            cur_point.x = 0;
            row_color[0] = 1 - row_color[0];
            row_color[1] = 1 - row_color[1];
            row_color[2] = 1 - row_color[2];

            // cout << "updated cur point y : " << cur_point.y << endl;

        }



    }

     double intersect(Ray *r, double* color_in, int level){

        Point normal(0.0, 0.0, 1.0);
        if(r->start.z < 0 ){
            normal.z = - normal.z;
        }
        double dot_d_n = dot_product(normal, r->dir);

        if( dot_d_n == 0 ){
            return INT_MAX;
        }

        double dot_p_n = reference_point.z * r->dir.z;
        double dot_r0_n = dot_product(r->start, normal);

        double t = ( dot_p_n - dot_r0_n ) * 1.0 / dot_d_n;

        if( (t < 0.0 ) && ( t > INT_MAX ) ){
            return INT_MAX;
        }

        Point intersection_point( r->start.x + ( t * r->dir.x ),r->start.y + ( t * r->dir.y ), r->start.z + ( t * r->dir.z ) );

        double temp_color[3];
        temp_color[0] = 1;
        temp_color[1] = 1;
        temp_color[2] = 1;

        if( ( intersection_point.y > reference_point.y ) && ( intersection_point.y < - reference_point.y ) ){
                int col_idx = int( (intersection_point.y - reference_point.y) / length );
                if( ( intersection_point.x > reference_point.x ) && ( intersection_point.x < - reference_point.x ) )
                {
                    int row_idx = int( (intersection_point.x - reference_point.x) / length );
                    if( ( row_idx + col_idx ) % 2 == 0){
                        temp_color[0] = 0;
                        temp_color[1] = 0;
                        temp_color[2] = 0;
                    }



                    if(level == 0) return t;


                    color[0] = temp_color[0];
                    color[1] = temp_color[1];
                    color[2] = temp_color[2];

                    setLightingColor(color_in); // set ambient color
                    // level != 0
                    double epsilon = 0.00001;

                    Point intersect_point = intersection_point;

                    vector<PointLight*> all_pl = pointLights;


                    for( auto &sl : spotlights ){
                        PointLight* pl = sl->pl;
                        // cout << "hello " << endl;
                        Ray* ray_point_light = new Ray(pl->light_pos, intersect_point);
                        // ray_point_light->printRay();

                        // check if crosses cutoff angle
                        double alpha = dot_product(ray_point_light->dir, sl->light_direction);
                        double light_dir_len = sl->light_direction.length();
                        double ray_len = ray_point_light->dir.length();
                        alpha = acos( alpha * 1.0 / (light_dir_len * ray_len) ) * 180.0 / pi;
                        // cout << "angle : " << alpha << endl;
                        if(alpha > sl->cutoff_angle ){
                            //cout << "not spotlight added " << endl;
                            continue;
                        }
                        else{
                            // cout << "spotlight added " << endl;
                            all_pl.push_back(pl);
                        }
                    }

                    for( auto &pl : all_pl ){
                        // cout << "hello " << endl;
                        Ray* ray_point_light = new Ray(pl->light_pos, intersect_point);
                        if(ray_point_light->start.z < 0.0 ){
                            normal.z = - 1.0;
                        }
                        // ray_point_light->printRay();
                        double min_t_pl = INT_MAX;
                        for( auto &obj : Objects ){
                            double* dummyColor = new double[3];
                            // cout << "here" << endl;
                            double t_pl = obj->intersect(ray_point_light, dummyColor, 0);
                            // cout << "t_pl : " << t_pl << endl;
                            if( ( t_pl > epsilon ) && ( t_pl < min_t_pl ) ){
                                min_t_pl = t_pl;
                            }
                        }
                        if(min_t_pl != INT_MAX ){
                            // cout << "inside here "<< endl;
                            // cout << "Min t pl : " << min_t_pl << endl;
                            Point shadow_intersect_point( ray_point_light->start.x + ( min_t_pl * ray_point_light->dir.x ), ray_point_light->start.y + ( min_t_pl * ray_point_light->dir.y ), ray_point_light->start.z + ( min_t_pl * ray_point_light->dir.z ) );

                            // check if the intersection point is in shadow
                            double dis_shadow_ip = sqrt( pow(shadow_intersect_point.x - ray_point_light->start.x, 2) + pow(shadow_intersect_point.y - ray_point_light->start.y, 2) + pow(shadow_intersect_point.z - ray_point_light->start.z, 2) );
                            double dis_ip = sqrt( pow(intersect_point.x - ray_point_light->start.x, 2) + pow(intersect_point.y - ray_point_light->start.y, 2) + pow(intersect_point.z - ray_point_light->start.z, 2) );


                            // cout << "dis shadow  : " << dis_shadow_ip << endl;
                            // cout << "dis ip : " << dis_ip - epsilon << endl;
                            if( dis_shadow_ip < dis_ip - epsilon ){
                                continue;
                            }

                            // not in shadow, calculate specular_diffuse components
                            // cout << "hello" << endl;
                            calculateSpecularDiffuse(normal, ray_point_light, color_in, intersect_point, r, pl);

                        }

                    }


                    if(level >= level_recursion ) return t;

                    double dot_ray_n = dot_product(normal, r->dir);
                    Point reflectedRayDir( - (2.0 * dot_ray_n * normal.x ) + r->dir.x, - (2.0 * dot_ray_n * normal.y ) + r->dir.y, - (2.0 * dot_ray_n * normal.z ) + r->dir.z);
                    reflectedRayDir.normalize();

                    Point reflect_initial(intersect_point.x + reflectedRayDir.x , intersect_point.y + reflectedRayDir.y, intersect_point.z + reflectedRayDir.z);
                    Ray* reflectedRay = new Ray();

                    reflectedRay->setStart(reflect_initial);
                    reflectedRay->setDir(reflectedRayDir);

                    // ray.printRay();

                    double *color_ray = new double[3];

                    Object *nearest = NULL;
                    int near_idx;
                    double t_min = INT_MAX;
                    int count = 0;

                     for (auto & element : Objects) {
                        // element->intersect();
                        double* color_temp = new double[3];
                        double t = element->intersect(reflectedRay, color_temp, 0);
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
                        // cout << near_idx << endl;
                        nearest = Objects[near_idx];
                        t_min = nearest->intersect(reflectedRay, color_ray, level + 1);
                        // cout << color_ray[0] << " " << color_ray[1] << " " << color_ray[2] << endl;
                    }

                    calculateReflection(color_ray, color_in);

                   return t;

                }

        }

        return INT_MAX;
    }

};




