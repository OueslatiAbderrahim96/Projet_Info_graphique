// Raytracer.cpp : Defines the entry point for the console application.
#define _CRT_SECURE_NO_WARNINGS // for Visual Studio 2017 (maybe 2015 as well)
#include <omp.h>
#include <iostream>
#include <stdio.h>
#include <vector>
#include <algorithm>
#include <limits>
#include <math.h>
#include <random>
#include <string>
#include <map>
#include <cstdio>
#include <list>
//email : nicolas.bonneel@liris.cnrs.fr
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

#define MATH_PI 3.1415922653589

//les textures : perso.liris.cnrs.fr/nbonneel/tmp/texturesbmp.zip
using namespace std;
std::default_random_engine engine[8];
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
     Vector& operator+=(const Vector& b){
         coords[0]+=b[0];
         coords[1]+=b[1];
         coords[2]+=b[2];
         return *this;
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
 Vector operator*(const Vector& B , const Vector& A){
     return Vector(A[0]*B[0], A[1]*B[1], A[2]*B[2]);
 }
 Vector operator/(const Vector& B , double k){
     return Vector(B[0]/k, B[1]/k, B[2]/k);
 }


 double dot(const Vector& A, const Vector& B){
     return A[0] * B[0] + A[1] * B[1] + A[2] * B[2];
 }

 Vector cross(const Vector& A, const Vector& B){
     return Vector(A[1] * B[2] - A[2] * B[1], A[2] * B[0] - A[0] * B[2], A[0] * B[1] - A[1] * B[0]);
 }

Vector random_cos(const Vector& N){
     double r1 = distrib(engine[omp_get_thread_num()]);
     double r2 = distrib(engine[omp_get_thread_num()]);
     Vector direction_aleatoir_repere_local(cos(2 * MATH_PI * r1) * sqrt(1 - r2), sin(2 * MATH_PI * r1) * sqrt(1-r2),sqrt(r2));
     Vector aleatoir(distrib(engine[omp_get_thread_num()])-0.5,distrib(engine[omp_get_thread_num()])-0.5,distrib(engine[omp_get_thread_num()])-0.5);
     Vector tangent1 = cross(N,aleatoir);
     tangent1.normalize();
     Vector tangent2 = cross(tangent1, N);                            
     return direction_aleatoir_repere_local[2] * N + direction_aleatoir_repere_local[0] * tangent1 + direction_aleatoir_repere_local[1] * tangent2;
}

class Ray{
public :
     Ray(const Vector& C, Vector u):C(C), u(u){};
     Vector C, u;
 };

class Object{
 public :
        Object(){};
 	virtual bool intersect(const Ray& d, Vector& P, Vector& N, double& t) const =0;

        Vector albedo;
        bool mirro,transp;
};

class Triangle : public Object{
public :
const Vector A, B, C;
Triangle(const Vector& A, const Vector& B, const Vector& C, const Vector& albedo, bool mirro = false, bool transp = false) :A(A), B(B), C(C){
	this->albedo = albedo;
	this->mirro = mirro; 
	this->transp = transp;
        
};

bool intersect(const Ray& d, Vector& P, Vector& N, double& t) const {
	double alpha, beta, gamma;
	return intersect(d,P,N,t,alpha, beta, gamma);
}

bool intersect(const Ray& d, Vector& P, Vector& N, double& t,double &alpha,double &beta,double &gamma) const {
	N = -1*cross(B - A, C - A);
        //cout<<A[0]<<" : "<<A[1]<<" : "<<A[2]<<endl;
        N.normalize();
        t = dot(C-d.C,N) / dot(d.u,N);
        //cout<<"N = "<<N.coords[0]<<"t = "<<t<<endl;
        if(t<0) return false;
        P = d.C + t*d.u;
        Vector u = B-A;
        Vector v = C-A;
        Vector w = P-A;
        double m11 = u.getNormSquare();
        double m12 = dot(u,v);
        double m22 = v.getNormSquare();
        double detm = m11*m22 - m12*m12;
        double b11 = dot(w,u);
        double b21 = dot(w, v);
        double detb=b11*m22 - b21*m12;
        beta = detb/detm;
        double detg = m11*b21 - m12*b11;
        gamma = detg/detm;
        alpha = 1-beta-gamma;
        if(alpha<0 || alpha>1) return false;
        if(beta<0 || beta>1) return false;
        if(gamma<0 || gamma>1) return false;
                
        return true;       
}


};


class Sphere : public Object{
public:
    Sphere(const Vector& O, double R, const Vector& albedo, bool mirro = false, bool transp = false): O(O),R(R){this->albedo = albedo;
this->mirro = mirro; 
this->transp = transp;};


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
    
};


class TriangleIndices {
public:
	TriangleIndices(int vtxi = -1, int vtxj = -1, int vtxk = -1, int ni = -1, int nj = -1, int nk = -1, int uvi = -1, int uvj = -1, int uvk = -1) : vtxi(vtxi), vtxj(vtxj), vtxk(vtxk), uvi(uvi), uvj(uvj), uvk(uvk), ni(ni), nj(nj), nk(nk) {
	};
	int vtxi, vtxj, vtxk;
	int uvi, uvj, uvk;
	int ni, nj, nk;
	int faceGroup;
};

class BBox {
public:
   BBox(){};
   BBox(const Vector& bmin, const Vector& bmax):bmin(bmin),bmax(bmax){};

   bool intersect(const Ray& d) const{
	double t_1_x = (bmin[0] - d.C[0])/d.u[0];
	double t_2_x = (bmax[0] - d.C[0])/d.u[0];
	double t_min_x = std::min(t_1_x,t_2_x);
	double t_max_x = std::max(t_1_x,t_2_x);

	double t_1_y = (bmin[1] - d.C[1])/d.u[1];
	double t_2_y = (bmax[1] - d.C[1])/d.u[1];
	double t_min_y = std::min(t_1_y,t_2_y);
	double t_max_y = std::max(t_1_y,t_2_y);

	double t_1_z = (bmin[2] - d.C[2])/d.u[2];
	double t_2_z = (bmax[2] - d.C[2])/d.u[2];
	double t_min_z = std::min(t_1_z,t_2_z);
	double t_max_z = std::max(t_1_z,t_2_z);
        
        if((std::min(std::min(t_max_x,t_max_y),t_max_z) - std::max(std::max(t_min_x,t_min_y),t_min_z)) > 0)
		return true;
	return false;
   }
   Vector bmin, bmax;
};

class BVH{
public:
	
	int i0, i1;
	BBox bbox;

	BVH* fg; BVH* fd;
};

class Geometry :public Object{
private : 
	BVH bvh;
public:
	Geometry() {};
	Geometry(const char* obj, double scaling, const Vector& offset, const Vector& albedo, bool mirro=false,bool transp=false) {
                this->albedo = albedo;
                this->mirro = mirro;
                this->transp = transp;
		readOBJ(obj);
		for (int i = 0; i < vertices.size(); i++) {
			vertices[i] = vertices[i] * scaling + offset;
		}
	}
       
       bool intersect(const Ray& d, Vector& P, Vector& N, double& t) const {
		t=1E99;
		bool has_inter=false;

		if (!bvh.bbox.intersect(d)) 
			return false;

		std::list<const BVH*> listBVH;
		listBVH.push_front(&bvh);
		while(!listBVH.empty()){
			const BVH* currentNode = listBVH.front();
			listBVH.pop_front();
			if(currentNode->fg && currentNode->fg->bbox.intersect(d)){
				listBVH.push_back(currentNode->fg);
			}
			if(currentNode->fd && currentNode->fd->bbox.intersect(d)){
				listBVH.push_back(currentNode->fd);
			}
			if(!currentNode->fg){                                
				for(int i=currentNode->i0; i<currentNode->i1;i++){                                       
					Triangle tri1(vertices[indices[i].vtxi],vertices[indices[i].vtxj],vertices[indices[i].vtxk],albedo,mirro,transp);
/*cout<<(vertices[indices[i].vtxi].coords[0]/25)<<" : "<<(vertices[indices[i].vtxi].coords[1]+20)/25<<" : "<<(vertices[indices[i].vtxi].coords[2]+55)/25<<endl;*/
					Vector localP, localN ;
					double localt;
					double alpha,beta, gamma;
					if (tri1.intersect(d,localP,localN,localt,alpha,beta, gamma)){
						has_inter=true;
						if (localt<t){
							t =localt;
							P=localP;
							//N=localN;	
							N = normals[indices[i].ni] * alpha +normals[indices[i].nj] * beta + normals[indices[i].nk] * gamma;
							N.normalize();
						}
					}
				}
			}
                   
		}

		/*if (!bb.intersect(d)) 
			return false;
		t=1E99;
		bool has_inter=false;
		for (int i=0;i<indices.size();i++){
			Triangle tri(vertices[indices[i].vtxi],vertices[indices[i].vtxj],vertices[indices[i].vtxk],albedo,mirro,transp);
			Vector localP, localN ;
			double localt;
			if (tri.intersect(d,localP,localN,localt)){
				has_inter=true;
				if (localt<t){
					t =localt;
					P=localP;
					N=localN;		
				}
			}
		
		}*/
		return has_inter;
	}
      
	void readOBJ(const char* obj) {

		char matfile[255];
		char grp[255];

		FILE* f;
		f = fopen(obj, "r");

		std::map<std::string, int> groupNames;
		int curGroup = -1;
		while (!feof(f)) {
			char line[255];
			if (!fgets(line, 255, f)) break;

			std::string linetrim(line);
			linetrim.erase(linetrim.find_last_not_of(" \r\t") + 1);
			strcpy(line, linetrim.c_str());

			if (line[0] == 'u' && line[1] == 's') {
				sscanf(line, "usemtl %[^\n]\n", grp);
				if (groupNames.find(std::string(grp)) != groupNames.end()) {
					curGroup = groupNames[std::string(grp)];
				}
				else {
					curGroup = groupNames.size();
					groupNames[std::string(grp)] = curGroup;
				}
			}
			if (line[0] == 'm' && line[1] == 't' && line[2] == 'l') {
				sscanf(line, "mtllib %[^\n]\n", matfile);
			}
			if (line[0] == 'v' && line[1] == ' ') {
				Vector vec;
				Vector col;
				if (sscanf(line, "v %lf %lf %lf %lf %lf %lf\n", &vec.coords[0], &vec.coords[2], &vec.coords[1], &col.coords[0], &col.coords[1], &col.coords[2]) == 6) {
					vertices.push_back(vec);
					vertexcolors.push_back(col);
				}
				else {
					sscanf(line, "v %lf %lf %lf\n", &vec.coords[0], &vec.coords[2], &vec.coords[1]);  // helmet
																				 //vec[2] = -vec[2]; //car2
					vertices.push_back(vec);
				}
			}
			if (line[0] == 'v' && line[1] == 'n') {
				Vector vec;
				sscanf(line, "vn %lf %lf %lf\n", &vec.coords[0], &vec.coords[2], &vec.coords[1]); //girl
				normals.push_back(vec);
			}
			if (line[0] == 'v' && line[1] == 't') {
				Vector vec;
				sscanf(line, "vt %lf %lf\n", &vec.coords[0], &vec.coords[1]);
				uvs.push_back(vec);
			}
			if (line[0] == 'f') {
				TriangleIndices t;
				int i0, i1, i2, i3;
				int j0, j1, j2, j3;
				int k0, k1, k2, k3;
				int nn;

				char* consumedline = line + 1;
				int offset;
				t.faceGroup = curGroup;
				nn = sscanf(consumedline, "%u/%u/%u %u/%u/%u %u/%u/%u%n", &i0, &j0, &k0, &i1, &j1, &k1, &i2, &j2, &k2, &offset);
				if (nn == 9) {
					if (i0 < 0) t.vtxi = vertices.size() + i0; else	t.vtxi = i0 - 1;
					if (i1 < 0) t.vtxj = vertices.size() + i1; else	t.vtxj = i1 - 1;
					if (i2 < 0) t.vtxk = vertices.size() + i2; else	t.vtxk = i2 - 1;
					if (j0 < 0) t.uvi = uvs.size() + j0; else	t.uvi = j0 - 1;
					if (j1 < 0) t.uvj = uvs.size() + j1; else	t.uvj = j1 - 1;
					if (j2 < 0) t.uvk = uvs.size() + j2; else	t.uvk = j2 - 1;
					if (k0 < 0) t.ni = normals.size() + k0; else	t.ni = k0 - 1;
					if (k1 < 0) t.nj = normals.size() + k1; else	t.nj = k1 - 1;
					if (k2 < 0) t.nk = normals.size() + k2; else	t.nk = k2 - 1;

					indices.push_back(t);
				}
				else {
					nn = sscanf(consumedline, "%u/%u %u/%u %u/%u%n", &i0, &j0, &i1, &j1, &i2, &j2, &offset);
					if (nn == 6) {
						if (i0 < 0) t.vtxi = vertices.size() + i0; else	t.vtxi = i0 - 1;
						if (i1 < 0) t.vtxj = vertices.size() + i1; else	t.vtxj = i1 - 1;
						if (i2 < 0) t.vtxk = vertices.size() + i2; else	t.vtxk = i2 - 1;
						if (j0 < 0) t.uvi = uvs.size() + j0; else	t.uvi = j0 - 1;
						if (j1 < 0) t.uvj = uvs.size() + j1; else	t.uvj = j1 - 1;
						if (j2 < 0) t.uvk = uvs.size() + j2; else	t.uvk = j2 - 1;
						indices.push_back(t);
					}
					else {
						nn = sscanf(consumedline, "%u %u %u%n", &i0, &i1, &i2, &offset);
						if (nn == 3) {
							if (i0 < 0) t.vtxi = vertices.size() + i0; else	t.vtxi = i0 - 1;
							if (i1 < 0) t.vtxj = vertices.size() + i1; else	t.vtxj = i1 - 1;
							if (i2 < 0) t.vtxk = vertices.size() + i2; else	t.vtxk = i2 - 1;
							indices.push_back(t);
						}
						else {
							nn = sscanf(consumedline, "%u//%u %u//%u %u//%u%n", &i0, &k0, &i1, &k1, &i2, &k2, &offset);
							if (i0 < 0) t.vtxi = vertices.size() + i0; else	t.vtxi = i0 - 1;
							if (i1 < 0) t.vtxj = vertices.size() + i1; else	t.vtxj = i1 - 1;
							if (i2 < 0) t.vtxk = vertices.size() + i2; else	t.vtxk = i2 - 1;
							if (k0 < 0) t.ni = normals.size() + k0; else	t.ni = k0 - 1;
							if (k1 < 0) t.nj = normals.size() + k1; else	t.nj = k1 - 1;
							if (k2 < 0) t.nk = normals.size() + k2; else	t.nk = k2 - 1;
							indices.push_back(t);
						}
					}
				}


				consumedline = consumedline + offset;

				while (true) {
					if (consumedline[0] == '\n') break;
					if (consumedline[0] == '\0') break;
					nn = sscanf(consumedline, "%u/%u/%u%n", &i3, &j3, &k3, &offset);
					TriangleIndices t2;
					t2.faceGroup = curGroup;
					if (nn == 3) {
						if (i0 < 0) t2.vtxi = vertices.size() + i0; else	t2.vtxi = i0 - 1;
						if (i2 < 0) t2.vtxj = vertices.size() + i2; else	t2.vtxj = i2 - 1;
						if (i3 < 0) t2.vtxk = vertices.size() + i3; else	t2.vtxk = i3 - 1;
						if (j0 < 0) t2.uvi = uvs.size() + j0; else	t2.uvi = j0 - 1;
						if (j2 < 0) t2.uvj = uvs.size() + j2; else	t2.uvj = j2 - 1;
						if (j3 < 0) t2.uvk = uvs.size() + j3; else	t2.uvk = j3 - 1;
						if (k0 < 0) t2.ni = normals.size() + k0; else	t2.ni = k0 - 1;
						if (k2 < 0) t2.nj = normals.size() + k2; else	t2.nj = k2 - 1;
						if (k3 < 0) t2.nk = normals.size() + k3; else	t2.nk = k3 - 1;
						indices.push_back(t2);
						consumedline = consumedline + offset;
						i2 = i3;
						j2 = j3;
						k2 = k3;
					}
					else {
						nn = sscanf(consumedline, "%u/%u%n", &i3, &j3, &offset);
						if (nn == 2) {
							if (i0 < 0) t2.vtxi = vertices.size() + i0; else	t2.vtxi = i0 - 1;
							if (i2 < 0) t2.vtxj = vertices.size() + i2; else	t2.vtxj = i2 - 1;
							if (i3 < 0) t2.vtxk = vertices.size() + i3; else	t2.vtxk = i3 - 1;
							if (j0 < 0) t2.uvi = uvs.size() + j0; else	t2.uvi = j0 - 1;
							if (j2 < 0) t2.uvj = uvs.size() + j2; else	t2.uvj = j2 - 1;
							if (j3 < 0) t2.uvk = uvs.size() + j3; else	t2.uvk = j3 - 1;
							consumedline = consumedline + offset;
							i2 = i3;
							j2 = j3;
							indices.push_back(t2);
						}
						else {
							nn = sscanf(consumedline, "%u//%u%n", &i3, &k3, &offset);
							if (nn == 2) {
								if (i0 < 0) t2.vtxi = vertices.size() + i0; else	t2.vtxi = i0 - 1;
								if (i2 < 0) t2.vtxj = vertices.size() + i2; else	t2.vtxj = i2 - 1;
								if (i3 < 0) t2.vtxk = vertices.size() + i3; else	t2.vtxk = i3 - 1;
								if (k0 < 0) t2.ni = normals.size() + k0; else	t2.ni = k0 - 1;
								if (k2 < 0) t2.nj = normals.size() + k2; else	t2.nj = k2 - 1;
								if (k3 < 0) t2.nk = normals.size() + k3; else	t2.nk = k3 - 1;
								consumedline = consumedline + offset;
								i2 = i3;
								k2 = k3;
								indices.push_back(t2);
							}
							else {
								nn = sscanf(consumedline, "%u%n", &i3, &offset);
								if (nn == 1) {
									if (i0 < 0) t2.vtxi = vertices.size() + i0; else	t2.vtxi = i0 - 1;
									if (i2 < 0) t2.vtxj = vertices.size() + i2; else	t2.vtxj = i2 - 1;
									if (i3 < 0) t2.vtxk = vertices.size() + i3; else	t2.vtxk = i3 - 1;
									consumedline = consumedline + offset;
									i2 = i3;
									indices.push_back(t2);
								}
								else {
									consumedline = consumedline + 1;
								}
							}
						}
					}
				}

			}


		}
		fclose(f);
		//cout<<vertices.size()<<endl;
		//BBox b = build_bbox(0,indices.size());
		build_bvh(&bvh,0,indices.size());
		
	}

	BBox build_bbox(int i0, int i1){
              BBox res;
		res.bmax = vertices[indices[i0].vtxi];
		res.bmin = vertices[indices[i0].vtxi];
		int index =-1;
            for(int i=i0;i<i1;i++){
		for(int j=0;j<3;j++){
			if(j==0) index = indices[i].vtxi;
			else{if(j==1) index = indices[i].vtxj;
				else index = indices[i].vtxk;}
			for(int k=0; k<3;k++){			
				res.bmin.coords[k] = std::min(res.bmin[k],vertices[index][k]);
				res.bmax.coords[k] = std::max(res.bmax[k],vertices[index][k]);
			}
		}
	     }
	return res;
	}

	void build_bvh(BVH* node, int i0, int i1){
		node->bbox = build_bbox(i0,i1);
		node->i0 = i0;
		node->i1 = i1;
		node->fg = NULL;
		node->fd = NULL;
		Vector diag = node->bbox.bmax - node->bbox.bmin;
		int split_dim;
		if((diag[0] >diag[1]) && (diag[0] >diag[2])){
			split_dim =0;
		}else{
			if((diag[1] >diag[0]) && (diag[1] >diag[2])){
				split_dim = 1;	
			}else{
				split_dim = 2;			
			}
		}

		double split_val = node->bbox.bmin[split_dim] + diag[split_dim]*0.5;

		int pivot = i0 - 1;
		for(int i =i0;i<i1;i++){
			double centre_split_dim = (vertices[indices[i].vtxi].coords[split_dim] + vertices[indices[i].vtxj].coords[split_dim] + vertices[indices[i].vtxk].coords[split_dim])/3.;		
			if(centre_split_dim < split_val){
				pivot++;                              
				//cout<<indices[i].vtxi<<" : "<<indices[pivot].vtxi<<endl;
				std::swap(indices[i].vtxi,indices[pivot].vtxi);
				std::swap(indices[i].vtxj,indices[pivot].vtxj);
				std::swap(indices[i].vtxk,indices[pivot].vtxk);
				//cout<<indices[i].vtxi<<" : "<<indices[pivot].vtxi<<endl;

				std::swap(indices[i].ni,indices[pivot].ni);
				std::swap(indices[i].nj,indices[pivot].nj);
				std::swap(indices[i].nk,indices[pivot].nk);
			}
		}

	if(pivot <= i0 || pivot>=i1) 
			return;

	node->fg = new BVH();
	build_bvh(node->fg, i0, pivot);

	node->fd = new BVH();
	build_bvh(node->fd, pivot, i1);
	
	}


	void add_texture(const char* filename) {

		textures.resize(textures.size() + 1);
		w.resize(w.size() + 1);
		h.resize(h.size() + 1);

		FILE* f;
		f = fopen(filename, "rb");
		unsigned char info[54];
		fread(info, sizeof(unsigned char), 54, f); // read the 54-byte header

		w[w.size() - 1] = *(int*)&info[18]; // extract image height and width from header
		h[h.size() - 1] = *(int*)&info[22];

		int size = 3 * w[w.size() - 1] * h[h.size() - 1];
		textures[textures.size() - 1].resize(size); // allocate 3 bytes per pixel
		fread(&textures[textures.size() - 1][0], sizeof(unsigned char), size, f); // read the rest of the data at once
		fclose(f);

		for (int i = 0; i < size; i += 3) {
			std::swap(textures[textures.size() - 1][i], textures[textures.size() - 1][i + 2]);
		}
	}

	std::vector<TriangleIndices> indices;
	std::vector<Vector> vertices;
	std::vector<Vector> normals;
	std::vector<Vector> uvs; // Vector en 3D mais on n'utilise que 2 composantes
	std::vector<Vector> vertexcolors;

	std::vector<std::vector<unsigned char> > textures;
	std::vector<int> w, h;

};


class Scene{
public :
    Scene(){};

    void addSphere(const Sphere& s){
        objects.push_back(&s);
    }
    void addTriangle(const Triangle& s){
        objects.push_back(&s);
    }
    void addGeometry(const Geometry& s){
        objects.push_back(&s);
    }

   bool intersect(const Ray& r, Vector& P, Vector& N, int& indice, double& t) const {
        bool has_inter = false;
        t = std::numeric_limits<double>::max();
        
        for(int i=0; i<objects.size();i++){
           Vector Plocal, Nlocal;
           double tlocal;
           bool interlocal = objects[i]->intersect(r,Plocal, Nlocal, tlocal);
         
           
           if(interlocal){
                has_inter = true;
                //cout<<t<<" ; "<<tlocal<<endl;
                if(tlocal < t){//if(i ==6) cout<<"ok : "<<i<<endl;
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
        if(has_intersection){//cout<<objects.size()<<" ;ok  "<<indice_sphere<<endl;
            if(indice_sphere == 0){
                   return L->albedo * intensiteL / (4* MATH_PI * L->R * L->R);
            }

            if(objects[indice_sphere]->mirro){
                Vector direction_miroir = r.u - 2*dot(r.u,N)*N;
                Ray rray_miroir(P+0.001*N, direction_miroir);
                intensite_pixel = getColor(rray_miroir, numrebond-1);
            }

            else {
		  if(objects[indice_sphere]->transp){
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
                               Ray rayon_refracte(P - 0.001*NforTransp,direction_refracte);
                               intensite_pixel = getColor(rayon_refracte, numrebond-1);
                          }
		        
		    }
	            else{   //eclairage direct
                            /*Vector vlocal = (L-P);
                            vlocal.normalize();
			    Ray shadowRay(P+0.0001*N, vlocal);     
			    Vector Pprime, Nprime;
			    int indiceprime;
			    double tprime;
			    double d_lightlocal = (L-P).getNormSquare(); 
			    if(intersect(shadowRay,Pprime, Nprime, indiceprime,tprime) && tprime*tprime < d_lightlocal){
				 intensite_pixel = Vector(0.,0.,0.);
			    }
			    else{
				intensite_pixel = spheres[indice_sphere]->albedo / MATH_PI * intensiteL * std::max(0., dot(vlocal,N)) / d_lightlocal;
			     }*/
                             
                             
                             Vector axePO = (P - L->O);
                             axePO.normalize();
                             Vector dir_aleatoire = random_cos(axePO);
                             Vector point_aleatoire = dir_aleatoire * L->R + L->O;
                             Vector wi = (point_aleatoire - P);
                             wi.normalize();
                             double d_lightlocal = (point_aleatoire - P).getNormSquare();
                             //Vector Np = dir_aleatoire;

                           
			    Ray shadowRay(P+0.001*N, wi);     
			    Vector Pprime, Nprime;
			    int indiceprime;
			    double tprime;
                            bool has_inter = intersect(shadowRay,Pprime, Nprime, indiceprime,tprime);
                             
                              if(has_inter && tprime*tprime < d_lightlocal*0.99){
				 intensite_pixel = Vector(0.,0.,0.);
			    }else{//if(indice_sphere==6) cout<<"ok1 : "<<indice_sphere<<endl;
                             intensite_pixel = (intensiteL / (4*MATH_PI*d_lightlocal)*std::max(0.,dot(N,wi))*dot(dir_aleatoire,-1*wi) /dot(axePO, dir_aleatoire))* objects[indice_sphere]->albedo ;                        
                             }
                             //eclairage indirect      
                             Vector direction_aleatoir = random_cos(N);
                             Ray rray_aleatoir(P+0.001*N, direction_aleatoir);
                             intensite_pixel += getColor(rray_aleatoir, numrebond-1) * objects[indice_sphere]->albedo;

		       }
        	}
	}
        return intensite_pixel;
    }


    std::vector<const Object*> objects;
    Sphere* L ;
    double intensiteL;
};

int main()
{
    int W = 512;
    int H = 512;
    Scene s;
    Sphere s_lumiere(Vector(15, 70, -30),15.,Vector(1.,1.,1.));
    /*Sphere s1(Vector(0., 0., -55.), 10, Vector(1.,1.,1.));
    Sphere s2(Vector(-15., 0., -35.), 10, Vector(1.,1.,1.),false,true);
    Sphere s3(Vector(15., 0., -75.), 10, Vector(1.,1.,1.),true);*/
    Geometry g1("BeautifulGirl.obj",25.,Vector(0.,-10.,-20.),Vector(1.,1.,1.));
    Sphere ssol(Vector(0., -2000-20, 0.), 2000, Vector(1.,1.,1.));
    Sphere splafond(Vector(0., 2000+100, 0.), 2000, Vector(1.,1.,1.));
    Sphere smurgauche(Vector(-2000-50, 0., 0.), 2000, Vector(1.,0.,1.));
    Sphere smurdroit(Vector(2000+50, 0., 0.), 2000, Vector(0.,1.,0.));
    Sphere smurfond(Vector(0., 0., -2000-100), 2000, Vector(0.,1.,1.));
    
    //Triangle tri(Vector(-10,-10,-55),Vector(10,-10,-55),Vector(0,10,-55),Vector(1,0,0)); 

    s.addSphere(s_lumiere);    
    /*s.addSphere(s1);
    s.addSphere(s2);
    s.addSphere(s3);*/
    s.addGeometry(g1);
    s.addSphere(ssol);
    s.addSphere(splafond);
    s.addSphere(smurgauche);
    s.addSphere(smurdroit);
    s.addSphere(smurfond);
    //s.addTriangle(tri);
    

    s.L = &s_lumiere;
    s.intensiteL = 3000000000;

    double alphafov = 60 * MATH_PI /180;
    double d = W / (2 * tan(alphafov/2.));
    const int nbr_ray =10; 
    Vector position_camera(0.,0.,0.);
    double focus_distance = 35.;
    double aperture=0.5;
    std::vector<unsigned char> img(W*H * 3, 0);
#pragma omp parallel for  schedule(dynamic,1)
    for(int i = 0 ; i<H ; i++){
        
        for(int j = 0; j < W; j++){
                
                Vector I = Vector(0.,0.,0.);

		for(int k =0;k<nbr_ray;k++){                       
                        //Box Muller
                        double r1 = distrib(engine[omp_get_thread_num()]);
                        double r2 = distrib(engine[omp_get_thread_num()]);
                        double R = sqrt(-2 * log(r1));
                        double dx = R * cos(2 * MATH_PI * r2);
                        double dy = R * sin(2 * MATH_PI * r2);
                           
                        double dx_aperture = (distrib(engine[omp_get_thread_num()])-0.5)*aperture;
                        double dy_aperture = (distrib(engine[omp_get_thread_num()])-0.5)*aperture;
                        Vector u(j-W / 2. + 0.5 + dx, -i+H /2. + 0.5 + dy, -d);
                        u.normalize();
                        Vector destination = position_camera + focus_distance * u; 
                        Vector new_origin = position_camera + Vector(dx_aperture,dy_aperture,0);
                        Vector new_direction = (destination - new_origin);
                        new_direction.normalize();
                        Ray r(new_origin  , new_direction);
                 	I += s.getColor(r, 5)/nbr_ray;
                }

                img[(i*W + j)*3 + 0] = std::min(255., pow(I[0], 0.45));
                img[(i*W + j)*3 + 1] = std::min(255., pow(I[1], 0.45));
                img[(i*W + j)*3 + 2] = std::min(255., pow(I[2], 0.45));
           }
        }
    

    //img[(10 * W + 50)*3] = 255;  // pixel at (x, y) = (50, 10) is red

    stbi_write_png("image16.png", W, H,3, &img[0], 0);

    return 0;
}
