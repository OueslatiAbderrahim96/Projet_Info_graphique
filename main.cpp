// Raytracer.cpp : Defines the entry point for the console application.
#define _CRT_SECURE_NO_WARNINGS // for Visual Studio 2017 (maybe 2015 as well)

#include <iostream>
#include <stdio.h>
#include <vector>
#include <algorithm>
#include <limits>
#include <math.h>
#include <random>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

#define MATH_PI 3.1415922653589

std::default_random_engine engine;
std::uniform_real_distribution<double> distrib(0,1);

inline double sqr(double x){
    return x*x;
}

void save_image(const char* filename, const unsigned char* tableau, int w, int h) { // (0,0) is top-left corner

    FILE *f;

    int filesize = 54 + 3 * w*h;

    unsigned char bmpfileheader[14] = { 'B','M', 0,0,0,0, 0,0, 0,0, 54,0,0,0 };
    unsigned char bmpinfoheader[40] = { 40,0,0,0, 0,0,0,0, 0,0,0,0, 1,0, 24,0 };
    unsigned char bmppad[3] = { 0,0,0 };

    bmpfileheader[2] = (unsigned char)(filesize);
    bmpfileheader[3] = (unsigned char)(filesize >> 8);
    bmpfileheader[4] = (unsigned char)(filesize >> 16);
    bmpfileheader[5] = (unsigned char)(filesize >> 24);

    bmpinfoheader[4] = (unsigned char)(w);
    bmpinfoheader[5] = (unsigned char)(w >> 8);
    bmpinfoheader[6] = (unsigned char)(w >> 16);
    bmpinfoheader[7] = (unsigned char)(w >> 24);
    bmpinfoheader[8] = (unsigned char)(h);
    bmpinfoheader[9] = (unsigned char)(h >> 8);
    bmpinfoheader[10] = (unsigned char)(h >> 16);
    bmpinfoheader[11] = (unsigned char)(h >> 24);

    f = fopen(filename, "wb");
    fwrite(bmpfileheader, 1, 14, f);
    fwrite(bmpinfoheader, 1, 40, f);
    unsigned char *row = new unsigned char[w * 3];
    for (int i = 0; i < h; i++)
    {
        for (int j = 0; j < w; j++) {
            row[j * 3] = tableau[(w*(h - i - 1) * 3) + j * 3+2];
            row[j * 3+1] = tableau[(w*(h - i - 1) * 3) + j * 3+1];
            row[j * 3+2] = tableau[(w*(h - i - 1) * 3) + j * 3];
        }
        fwrite(row, 3, w, f);
        fwrite(bmppad, 1, (4 - (w * 3) % 4) % 4, f);
    }
    fclose(f);
    delete[] row;
}
 class Vector {
 public :
     Vector(double x=0, double y=0, double z=0){
         coords[0] = x;
         coords[1] = y;
         coords[2] = z;
     }
      double operator[](int i) const {
        return coords[i];
     }
     double operator[](int i){
        return coords[i];
     }
     double getNormSquare() const{
         return coords[0] * coords[0] + coords[1] * coords[1] + coords[2] * coords[2] ;
     }
     void operator*=(double a){
        coords[0]*=a;
        coords[1]*=a;
        coords[2]*=a;
     }
     void operator/=(double a){
        coords[0]/=a;
        coords[1]/=a;
        coords[2]/=a;
     }
     void normalize(){
         double aux = sqrt(getNormSquare());
         coords[0] /= aux;
         coords[1] /= aux;
         coords[2] /= aux;
     }
     double coords[3];
 };
 Vector operator+(const Vector& A, const Vector& B ){
     return Vector(A[0]+B[0], A[1]+B[1], A[2]+B[2]);
 }
 Vector operator-(const Vector& A, const Vector& B ){
     return Vector(A[0]-B[0], A[1]-B[1], A[2]-B[2]);
 }

 Vector operator*(double k, const Vector& B ){
     return Vector(k*B[0], k*B[1], k*B[2]);
 }
  Vector operator*(const Vector& B , double k){
     return Vector(k*B[0], k*B[1], k*B[2]);
 }
  Vector operator/(const Vector& B , double k){
     return Vector(B[0]/k, B[1]/k, B[2]/k);
 }

 double dot(const Vector& A, const Vector& B){
     return A[0] * B[0] + A[1] * B[1] + A[2] * B[2];
 }

 class Ray{
public :
     Ray(const Vector& C, Vector u):C(C), u(u){};
     Vector C, u;
 };

class Sphere{
public:
    Sphere(const Vector& O, double R, const Vector& albedo, bool mirro = false, bool transp = false): O(O),R(R), albedo(albedo),mirro(mirro),transp(transp){};


    bool intersect(const Ray& r, Vector& P, Vector& N, double& t) const {
        //2*t^2+b*t+c
        double a = 1;
        double b = 2 * dot(r.u, r.C - O);
        double c =(r.C-O).getNormSquare() -R*R;

        double delta = b*b - 4*a*c;
        if(delta<0) return false;
        double sqrtDelta = sqrt(delta);
        double t1 = (-b + sqrtDelta)/(2*a);
        if(t1<0) return false;

        double t0 = (-b - sqrtDelta)/(2*a);
        if(t0 > 0){
            t = t0;
        }else{
            t = t1;
        }
        P = r.C + t*r.u;
        N = (P - O);
	N.normalize();
        return true;
    }
    Vector O;
    double R;
    Vector albedo;
    bool mirro,transp;
};

class Scene{
public :
    Scene(){};

    void addSphere(Sphere* s){
        spheres.push_back(s);
    }

    bool intersect(const Ray& r, Vector& P, Vector& N, int& indice, double& t) const {
        bool has_inter = false;
        t = std::numeric_limits<double>::max();

        for(int i=0; i<spheres.size();i++){
           Vector Plocal, Nlocal;
           double tlocal;
           bool interlocal = spheres[i]->intersect(r,Plocal, Nlocal, tlocal);
           if(interlocal){
                has_inter = true;
                if(tlocal < t){
                    t = tlocal;
                    P = Plocal;
                    N = Nlocal;
                    indice = i;
                }
           }
        }
       return has_inter; 
    }

    Vector getColor(const Ray& r,int numrebond){
        if(numrebond<0) return Vector(0.,0.,0.);
        Vector P, N;
        int indice_sphere;
	double t;
        bool has_intersection = intersect(r, P, N, indice_sphere,t);

	Vector intensite_pixel=Vector(0.,0.,0.);
        if(has_intersection){

            if(spheres[indice_sphere]->mirro){
                Vector direction_miroir = r.u - 2*dot(r.u,N)*N;
                Ray rray_miroir(P+0.0001*N, direction_miroir);
                return getColor(rray_miroir, numrebond-1);
            }

            else {
		  if(spheres[indice_sphere]->transp){
		        double n1=1;
		        double n2=1.4;
		        Vector NforTransp(N);
		        if(dot(r.u,N)>0){
		            std::swap(n1,n2);
		            NforTransp = -1*N;
		        }
		        double delta = 1-sqr(n1/n2)*(1-sqr(dot(r.u,NforTransp)));
                         if(delta>0){
                               Vector direction_refracte = (n1/n2)*(r.u - dot(r.u, NforTransp)*NforTransp) - NforTransp*sqrt(delta);
                               Ray rayon_refracte(P - 0.0001*NforTransp,direction_refracte);
                               return getColor(rayon_refracte, numrebond-1);
                          }
		        
		    }
	            else{
                            Vector vlocal = (L-P);
                            vlocal.normalize();
			    Ray shadowRay(P+0.0001*N, vlocal);     
			    Vector Pprime, Nprime;
			    int indiceprime;
			    double tprime;
			    double d_lightlocal = (L-P).getNormSquare(); 
			    if(intersect(shadowRay,Pprime, Nprime, indiceprime,tprime) && tprime*tprime < d_lightlocal){
				 return Vector(0.,0.,0.);
			    }
			    else{
				return spheres[indice_sphere]->albedo * intensiteL * std::max(0., dot(vlocal,N)) / d_lightlocal;
			     }
		       }
        	}
	}
        return intensite_pixel;
    }
    std::vector<Sphere*> spheres;
    Vector L ;
    double intensiteL;
};

int main()
{
    int W = 1024;
    int H = 1024;

    Scene s;
    Sphere s1(Vector(-15., 0., -55.), 10, Vector(1.,0.,0.),false,true);
    Sphere s2(Vector(15., 0., -55.), 10, Vector(1.,0.,0.),true);
    Sphere ssol(Vector(0., -2000-20, 0.), 2000, Vector(1.,1.,1.));
    Sphere splafond(Vector(0., 2000+100., 0.), 2000, Vector(1.,1.,1.));
    Sphere smurgauche(Vector(-2000-50., 0., 0.), 2000, Vector(0.,1.,0.));
    Sphere smurdroit(Vector(2000+50., 0., 0.), 2000, Vector(0.,0.,1.));
    Sphere smurfond(Vector(0., 0., -2000-100.), 2000, Vector(0.,1.,1.));

    s.addSphere(&s1);
    s.addSphere(&s2);
    s.addSphere(&ssol);
    s.addSphere(&splafond);
    s.addSphere(&smurgauche);
    s.addSphere(&smurdroit);
    s.addSphere(&smurfond);

    s.L=Vector(15, 70, -30);
    s.intensiteL = 1000000000;

    double alpha = 60 * MATH_PI /180;
    double d = W / (2 * tan(alpha/2.));

    std::vector<unsigned char> img(W*H * 3, 0);

    for(int i = 0 ; i<H ; i++){
        for(int j = 0; j < W; j++){
                Vector u(j-W / 2, -i+H /2, -d);
                u.normalize();

                Ray r(Vector(0.,0.,0.), u);
                Vector I=s.getColor(r, 5);

                img[(i*W + j)*3 + 0] = std::min(255., pow(I[0], 0.45));
                img[(i*W + j)*3 + 1] = std::min(255., pow(I[1], 0.45));
                img[(i*W + j)*3 + 2] = std::min(255., pow(I[2], 0.45));
           }
        }
    

    //img[(10 * W + 50)*3] = 255;  // pixel at (x, y) = (50, 10) is red

    stbi_write_png("image.png", W, H,3, &img[0], 0);

    return 0;
}
