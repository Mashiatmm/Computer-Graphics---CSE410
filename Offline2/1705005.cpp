#include "bitmap_image.hpp"
#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include <stack>
#include <cmath>
#include <iomanip>
#include <vector>
#include <time.h>

using namespace std;

struct Point{
    double x;
    double y;
    double z;
};

struct Triangle{
    Point point[3];
    double color[3];
};

double eye[3];
double look[3];
double up[3];
double fovY, aspectRatio, near, far;

double screen_Width, screen_Height, frontZ, rearZ, left_limit, right_limit, top_limit, bottom_limit;
double dx, dy, Top_Y, Left_X;
bitmap_image image;

double pi = 2 * acos(0.0);
vector< vector <double> > z_buffer;
vector<Triangle> triangles;
int tri_counter = 0;

ofstream outfile1;
ofstream outfile2;
ofstream outfile3;

ifstream infile;
ifstream configFile;
ifstream stagefile;


void readGlobalVar(){

    string line;

    getline( infile, line );
    istringstream iss(line);
    iss >> eye[0] >> eye[1] >> eye[2];
    //cout<< eyeX << eyeY << eyeZ <<endl;

    getline( infile, line );
    istringstream isl(line);
    isl >> look[0] >> look[1] >> look[2];
    //cout<< lookX << lookY << lookZ <<endl;

    getline( infile, line );
    istringstream isf(line);
    isf >> up[0] >> up[1] >> up[2];
    //cout<< upX << upY << upZ <<endl;

    getline( infile, line );
    istringstream isr(line);
    isr >> fovY >> aspectRatio >> near >> far;
    isr.clear();

    getline( configFile, line );
    isr.str(line);
    isr >> screen_Width >> screen_Height;
    isr.clear();

    getline( configFile, line );
    isr.str(line);
    isr >> left_limit;
    right_limit = - left_limit;
    isr.clear();

    getline( configFile, line );
    isr.str(line);
    isr >> bottom_limit;
    top_limit = - bottom_limit;
    isr.clear();

    getline( configFile, line );
    isr.str(line);
    isr >> frontZ >> rearZ;
    isr.clear();


}

double** IdentityMatrix(){
    double** identity = new double*[4];
    for( int i = 0 ; i < 4 ; i++ ){
        identity[i] = new double[4];
        for( int j = 0 ; j < 4 ; j++ ){
            if(i == j){
                identity[i][j] = 1;
            }
            else{
                identity[i][j] = 0;
            }
        }
    }
    return identity;
}



void print_matrix(double** matrix, int row, int column, ofstream& outfile){
    for( int i = 0 ; i < row ; i++ ){
        for( int j = 0 ; j < column ; j++ ){
            outfile << std::fixed << setprecision(7) << matrix[i][j] << " ";
        }
        outfile<<endl;
    }
    outfile<<endl;
}



double* transformPoint(double** transform_matrix, double* P ){
    //print_matrix(transform_matrix, 4, 4);
    //cout<<endl<<endl;
    double* P_prime = new double[4];
    double sum = 0;
    for(int i = 0 ; i < 4 ; i++){
        sum = 0;
        for(int j = 0 ; j < 4 ; j++){
            sum += transform_matrix[i][j] * P[j];
        }
        P_prime[i] = sum;
    }

    for(int i = 0 ; i < 4 ; i++ ){
        P_prime[i] = P_prime[i] / P_prime[3];
    }

    return P_prime;
}



double* R(double x[3], double a[3], double angle){
    double* temp = new double[3];

    double cross[3];
    double dot = a[0]*x[0] + a[1]*x[1] + a[2]*x[2];
    double rad_angle = angle * pi / 180.0;
    double cos_val = cos(rad_angle);
    double sin_val = sin(rad_angle);
    //outfile1<<"cos val " << cos_val << endl;
    //outfile1<<"sin val " << sin_val << endl;

    for(int  i = 0 ; i < 3 ; i++ ){
        cross[i] = a[ (i + 1) % 3 ] * x[ (i + 2) % 3 ] - x[ (i + 1) % 3 ] * a[ (i + 2) % 3 ];
        temp[i] = ( cos_val * x[i] ) + ( ( 1 - cos_val)  * dot * a[i] ) + ( sin_val * cross[i] ) ;

        //outfile1 << dot << " " << cross[i] << " " << temp[i] << endl;
    }

    return temp;

}

double** transform_matrix(double** transform_matrix, double** command_matrix){
    double** matrix = IdentityMatrix();
    double sum = 0;

    for(int i = 0 ; i < 4 ; i++){
        matrix[i] = new double[4];
        for(int j = 0 ; j < 4 ; j++){
            sum = 0;
            for( int k = 0 ; k < 4 ; k++){
                sum += transform_matrix[i][k] * command_matrix[k][j];
            }
            matrix[i][j] = sum;
        }

    }
    delete command_matrix;
    return matrix;
}

double** generate_rotation( double angle, double ax, double ay, double az){

    // normalize a
   //cout<<"generate rotation"<< endl;
    //cout << angle << " " << ax << " " << ay << " " << az << endl;
    double square_root_a = sqrt( ( ax*ax ) + ( ay*ay ) + ( az*az ) );
    ax = ax / square_root_a;
    ay = ay / square_root_a;
    az = az / square_root_a;
    // cout << angle << " " << ax << " " << ay << " " << az << endl;
    double a[3] = {ax, ay, az};
    double** matrix = IdentityMatrix();
    double sum = 0;

    for( int i = 0 ; i < 3 ; i++ ){
        double x[3] = {0, 0, 0};
        x[i] = 1;
        double* c = R(x, a, angle);

        for( int j = 0 ; j < 3 ; j++ ){
            matrix[j][i] = c[j];
        }

    }
    return matrix;
}

double** view_transformation(){
    double l[3], u[3], r[3];
    double** R = IdentityMatrix();

    for(int i = 0 ; i < 3 ; i++ ){
        l[i] = look[i] - eye[i];
    }

    double square_root_l = sqrt( ( l[0]*l[0] ) + ( l[1]*l[1] ) + ( l[2]*l[2] ) );
    l[0] = l[0] / square_root_l;
    l[1] = l[1] / square_root_l;
    l[2] = l[2] / square_root_l;

    for(int i = 0 ; i < 3 ; i++ ){
        r[i] = l[ (i + 1) % 3 ] * up[ (i + 2) % 3 ] - up[ (i + 1) % 3 ] * l[ (i + 2) % 3 ];
        R[2][i] = - l[i];
    }

    double square_root_r = sqrt( ( r[0]*r[0] ) + ( r[1]*r[1] ) + ( r[2]*r[2] ) );
    r[0] = r[0] / square_root_r;
    r[1] = r[1] / square_root_r;
    r[2] = r[2] / square_root_r;

    for(int i = 0 ; i < 3 ; i++ ){
        u[i] = r[ (i + 1) % 3 ] * l[ (i + 2) % 3 ] - l[ (i + 1) % 3 ] * r[ (i + 2) % 3 ];;
        R[1][i] = u[i];
        R[0][i] = r[i];
    }

    double** T = IdentityMatrix();
    for(int i = 0 ; i < 3 ; i++ ){
        T[i][3] = -eye[i];
    }

    return transform_matrix(R, T);

}

double** projection_transform(){
    double fovX = fovY * aspectRatio;
    double tan_angle = tan(fovY / 2.0 *  pi / 180.0);

    double t = near * tan_angle;
    tan_angle = tan(fovX / 2.0 * pi / 180.0);

    double r = near * tan_angle;

    double** projection_matrix = IdentityMatrix();

    projection_matrix[0][0] = near / r;
    projection_matrix[1][1] = near / t;
    projection_matrix[2][2] = -(far + near) / (far - near);
    projection_matrix[2][3] = -(2 * far * near) / (far - near);
    projection_matrix[3][2] = -1;
    projection_matrix[3][3] = 0;

    return projection_matrix;
}


//-------------------------------------------Task 4 Functions ------------------------------------------

void read_stagefile(){
    stagefile.open("stage3.txt");


    for(int i = 0 ; i < tri_counter ; i++)
    {
       Triangle tri;

       stagefile >> tri.point[0].x >> tri.point[0].y >> tri.point[0].z;
       stagefile >> tri.point[1].x >> tri.point[1].y >> tri.point[1].z;
       stagefile >> tri.point[2].x >> tri.point[2].y >> tri.point[2].z;

       for(int i=0;i<3;i++){
            tri.color[i] = rand()%256;
        }

        triangles.push_back(tri);

    }


    dx = ( right_limit - left_limit ) / screen_Width;
    dy = ( top_limit - bottom_limit ) / screen_Height;
    Top_Y = top_limit - ( dy / 2 );
    Left_X = left_limit + ( dx / 2 );

}

void initiate_z_buffer(){
    for (int i = 0; i < screen_Width; i++) {
        vector<double> v1;

        for (int j = 0; j < screen_Height; j++) {
            v1.push_back(rearZ);
        }


        z_buffer.push_back(v1);
    }

}



void print_triangles(){
    for (auto i = triangles.begin(); i != triangles.end(); ++i){
        cout << i->point[0].x << " " << i->point[0].y << " " << i->point[0].z << endl;
        cout << i->point[1].x << " " << i->point[1].y << " " << i->point[1].z << endl;
        cout << i->point[2].x << " " << i->point[2].y << " " << i->point[2].z << endl;
        cout << i->color[0]<< " " << i->color[1] << " " << i->color[2] << endl;

        cout << endl << endl;
    }

}

void print_tri(Triangle i){
    cout << i.point[0].x << " " << i.point[0].y << " " << i.point[0].z << endl;
    cout << i.point[1].x << " " << i.point[1].y << " " << i.point[1].z << endl;
    cout << i.point[2].x << " " << i.point[2].y << " " << i.point[2].z << endl;
    cout << i.color[0]<< " " << i.color[1] << " " << i.color[2] << endl;
}

void set_bg_color(){
     for(int i = 0 ; i < screen_Width ; i++){
        for(int j = 0 ; j < screen_Height ; j++){
            image.set_pixel(i, j, 0, 0, 0);
        }
    }
}

int find_intersect_col(Triangle tri, int scanline, int &left_col, int &right_col, double &xa, double &xb){
    int index = 0;

    double max_x = max({tri.point[0].x, tri.point[1].x, tri.point[2].x});
    double min_x = min({tri.point[0].x, tri.point[1].x, tri.point[2].x});

    double y = Top_Y - ( scanline * dy );
    double left_intersect_x = right_limit + 1;
    double right_intersect_x = left_limit - 1;
    double temp;

    double slope = ( tri.point[0].y - tri.point[1].y ) / ( tri.point[0].x - tri.point[1].x );
    double large_y = max(tri.point[0].y, tri.point[1].y);
    double small_y = min(tri.point[0].y, tri.point[1].y);

    if( ( y >= small_y ) && ( y <= large_y) ){
         temp = tri.point[0].x + ( ( y - tri.point[0].y ) / slope );
         left_intersect_x = min( temp, left_intersect_x);
         right_intersect_x = max( temp, right_intersect_x);
    }
    else if( ( slope == 0 ) && ( y == small_y) ){
        left_intersect_x = min(tri.point[0].x, tri.point[1].x);
        right_intersect_x = max(tri.point[0].x, tri.point[1].x);
    }
    else{
        index = 2;
    }

    slope = ( tri.point[1].y - tri.point[2].y ) / ( tri.point[1].x - tri.point[2].x );
    large_y = max(tri.point[2].y, tri.point[1].y);
    small_y = min(tri.point[2].y, tri.point[1].y);
    if(( y >= small_y ) && ( y <= large_y) ){
         temp = tri.point[1].x + ( ( y - tri.point[1].y ) / slope );
         left_intersect_x = min( temp, left_intersect_x);
         right_intersect_x = max( temp, right_intersect_x);
    }
    else if( ( slope == 0 ) && ( y == small_y) ){
        left_intersect_x = min(tri.point[2].x, tri.point[1].x);
        right_intersect_x = max(tri.point[2].x, tri.point[1].x);
    }
    else{
        index = 0;
    }

    slope = ( tri.point[0].y - tri.point[2].y ) / ( tri.point[0].x - tri.point[2].x );
    large_y = max(tri.point[2].y, tri.point[0].y);
    small_y = min(tri.point[2].y, tri.point[0].y);
    if(( y >= small_y ) && ( y <= large_y) ){
         temp = tri.point[0].x + ( ( y - tri.point[0].y ) / slope );
         left_intersect_x = min( temp, left_intersect_x);
         right_intersect_x = max( temp, right_intersect_x);
    }
    else if( ( slope == 0 ) && ( y == small_y) ){
        left_intersect_x = min(tri.point[0].x, tri.point[2].x);
        right_intersect_x = max(tri.point[0].x, tri.point[2].x);
    }
    else{
        index = 1;
    }

    // for clipping purpose
    left_intersect_x = max(left_limit, left_intersect_x);
    right_intersect_x = min(right_limit, right_intersect_x);

    double temp_x = ( left_intersect_x - left_limit ) / dx;
    left_col = floor(temp_x);

    temp_x = ( right_intersect_x - left_limit ) / dx;
    right_col = floor(temp_x);
    if( temp_x == right_col ){
        right_col -= 1;
    }

    xa = left_intersect_x;
    xb = right_intersect_x;

    return index;
}

void calc_za_zb(Triangle tri, int common_idx, int scanline, double &za, double &zb){
    double y = Top_Y - ( scanline * dy );
    int left_index = ( common_idx + 1 ) % 3;
    int right_index = ( common_idx + 2 ) % 3;

    if( tri.point[ ( common_idx + 1 ) % 3 ].x > tri.point[ ( common_idx + 2 ) % 3 ].x ){
        left_index = ( common_idx + 2 ) % 3;
        right_index = ( common_idx + 1 ) % 3;
    }

    double slope = (y - tri.point[common_idx].y) / ( tri.point[left_index ].y - tri.point[common_idx].y );
    za = tri.point[common_idx].z + slope * ( tri.point[ left_index ].z - tri.point[common_idx].z );

    slope = (y - tri.point[common_idx].y) / ( tri.point[ right_index ].y - tri.point[common_idx].y );
    zb = tri.point[common_idx].z + slope * ( tri.point[ right_index ].z - tri.point[common_idx].z );


}




void clip_scan(Triangle tri){
    print_tri(tri);
    // clipping
    double max_y = max({tri.point[0].y, tri.point[1].y, tri.point[2].y});
    double min_y = min({tri.point[0].y, tri.point[1].y, tri.point[2].y});


    int top_scanline = 0;
    if(top_limit > max_y){
        double temp = ( top_limit - max_y ) / dy ;
        top_scanline = floor( temp );
        // fully divisible, just on the row

    }
    int bottom_scanline = round( ( top_limit - bottom_limit ) / dy ) - 1;
    if( bottom_limit < min_y ){
        double temp =  ( top_limit - min_y ) / dy ;
        bottom_scanline = floor( temp );
        // fully divisible, just on the row
        if( bottom_scanline == temp ){
            bottom_scanline -= 1;
        }
    }

    // cout << "top scanline : " << top_scanline << " ; bottom scanline :  " << bottom_scanline << endl;
    for(int i = top_scanline; i <= bottom_scanline ; i++ ){

        int left_col, right_col;
        double za, zb, xa, xb;
        int common_idx = find_intersect_col( tri, i, left_col, right_col , xa, xb);

        double y = Top_Y - ( i * dy );


        // double za, zb;
        calc_za_zb(tri, common_idx, i, za, zb);

        // double xa =  Left_X + ( left_col * dx );
        // double xb =  Left_X + ( right_col * dx );

        double ys = Top_Y - ( i * dy );
        double xp;
        double zp = NULL;
        double constant = dx* (zb - za) / (xb - xa);

        // cout << "scanline : " << y << " ; left col : " << left_col << " ; right_col : " << right_col << endl;
        // cout << "zb : " << zb << " ; za : " << za << " ; constant : " << constant << endl;

        for( int j = left_col ; j <= right_col ; j++ ){
            if(zp == NULL ){
                xp =  Left_X + ( j * dx );
                // calc_z_val(tri, common_idx, i, j);
                double slope = (xp - xa) / (xb - xa);
                zp = za + ( slope * (zb - za) );
            }
            else{
                zp = zp + constant;
            }

            //  cout << " zp : " << zp << endl;
            if( ( zp >= frontZ ) && ( zp < z_buffer[i][j] ) ){
                z_buffer[i][j] = zp;

                // update pixel info
                image.set_pixel(j, i, tri.color[0], tri.color[1], tri.color[2]);
            }

        }
    }

    cout<<endl<<endl;
}

z_buffer_calc(){
    for (auto i = triangles.begin(); i != triangles.end(); ++i){
        clip_scan(*i);
    }
}


void save_z_file(){
    ofstream z_buffer_file;
    z_buffer_file.open("z_buffer.txt");

    for(auto i = z_buffer.begin() ; i<z_buffer.end() ; i++)
      {
         for(auto j = i->begin() ; j<i->end() ; j++){
            if(*j < rearZ ){
                z_buffer_file << std::fixed << setprecision(6) << *j <<"\t";
            }
         }

         z_buffer_file <<"\n";
         //similarly you can do other things
      }

    z_buffer_file.close();
}

int main () {

    // open a file in read mode.

    infile.open("scene.txt");
    configFile.open("config.txt");

    outfile1.open("stage1.txt");
    outfile2.open("stage2.txt");
    outfile3.open("stage3.txt");


    string line;

    readGlobalVar();

    double** view_transform_matrix = view_transformation();
    double** projection_matrix = projection_transform();

    string command;
    double** current_triangle = new double*[3];
    stack<double**> S;
    S.push( IdentityMatrix() );
    double** model_matrix = S.top();
    //print_matrix(pop_matrix);
    while ( getline( infile, command ) )
    {
        //cout<<command<<endl;

        if( command == "triangle" ){
            tri_counter += 1;
            for( int  i = 0 ; i < 3 ; i++ ){
                current_triangle[i] = new double[4];
                getline( infile, line );
                istringstream is_stream(line);
                is_stream >> current_triangle[i][0] >> current_triangle[i][1] >> current_triangle[i][2];
                current_triangle[i][3] = 1;
                current_triangle[i] = transformPoint(model_matrix, current_triangle[i] );
            }
            //print_matrix(model_matrix, 4, 4, outfile1);
            print_matrix(current_triangle, 3, 3, outfile1);

            for( int  i = 0 ; i < 3 ; i++ ){
                current_triangle[i] = transformPoint(view_transform_matrix, current_triangle[i] );
            }
            print_matrix(current_triangle, 3, 3, outfile2);

            for( int  i = 0 ; i < 3 ; i++ ){
                current_triangle[i] = transformPoint(projection_matrix, current_triangle[i] );
            }
            print_matrix(current_triangle, 3, 3, outfile3);

        }
        else if( command == "translate" ){
            double tx, ty, tz;
            getline( infile, line );
            istringstream is_stream(line);
            is_stream >> tx >> ty >> tz;
            double** translation_matrix = IdentityMatrix();
            translation_matrix[0][3] = tx;
            translation_matrix[1][3] = ty;
            translation_matrix[2][3] = tz;

            model_matrix = transform_matrix(model_matrix, translation_matrix);

        }
        else if( command == "scale" ){
            double sx, sy, sz;
            getline( infile, line );
            istringstream is_stream(line);
            is_stream >> sx >> sy >> sz;
            double** scale_matrix = IdentityMatrix();
            scale_matrix[0][0] = sx;
            scale_matrix[1][1] = sy;
            scale_matrix[2][2] = sz;

            model_matrix = transform_matrix(model_matrix, scale_matrix);

        }
        else if( command == "rotate" ){
            double angle, ax, ay, az;
            getline( infile, line );
            istringstream is_stream(line);
            is_stream >> angle >> ax >> ay >> az;

            double** rotation_matrix = generate_rotation(angle, ax, ay, az);
            model_matrix = transform_matrix(model_matrix, rotation_matrix);
            //outfile1<<"Rotation"<<endl;
            //print_matrix(rotation_matrix, 4, 4, outfile1);

        }
        else if( command == "push" ){
           S.push(model_matrix);
        }
        else if( command == "pop" ){
            model_matrix = S.top();
            S.pop();
        }
        else if( command == "end" ){
            break;
        }
    }

   // close the opened file.
   infile.close();
   outfile1.close();
   outfile2.close();
   outfile3.close();

   read_stagefile();
   print_triangles();
   initiate_z_buffer();
   image.setwidth_height(screen_Width, screen_Height);

   z_buffer_calc();
   image.save_image("out.bmp");
   save_z_file();

   stagefile.close();
   configFile.close();

   return 0;
}
