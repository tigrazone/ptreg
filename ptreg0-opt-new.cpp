#define _USE_MATH_DEFINES
#define erand48(dummy) (double(rand()) / RAND_MAX)




#include <math.h> // smallpt, a Path Tracer by Kevin Beason, 2008
#include <stdlib.h> // Make : g++ -O3 -fopenmp smallpt.cpp -o smallpt
#include <stdio.h>

#include <time.h>
clock_t t0, t1;


typedef struct Vec { // Usage: time ./explicit 16 && xv image.ppm
double x, y, z; // position, also color (r,g,b)
Vec(double x_=0, double y_=0, double z_=0){ x=x_; y=y_; z=z_; }
Vec operator+(const Vec &b) const { return Vec(x+b.x,y+b.y,z+b.z); }
Vec operator-(const Vec &b) const { return Vec(x-b.x,y-b.y,z-b.z); }
Vec operator*(double b) const { return Vec(x*b,y*b,z*b); }
Vec mult(const Vec &b) const { return Vec(x*b.x,y*b.y,z*b.z); }
Vec& norm(){ return *this = *this * (1/sqrt(x*x+y*y+z*z)); }
double dot(const Vec &b) const { return x*b.x+y*b.y+z*b.z; } // cross:
Vec operator%(Vec&b){return Vec(y*b.z-z*b.y,z*b.x-x*b.z,x*b.y-y*b.x);}
} const cVec;
struct Ray { Vec o, d; Ray(Vec o_, Vec d_) : o(o_), d(d_) {} };
enum Refl_t { DIFF, SPEC, REFR }; // material types, used in radiance()
struct Sphere {
double rad, rad2; // radius
Vec p, e, c; // position, emission, color
Refl_t refl; // reflection type (DIFFuse, SPECular, REFRactive)
Sphere(double rad_, Vec p_, Vec e_, Vec c_, Refl_t refl_):
rad(rad_), p(p_), e(e_), c(c_), refl(refl_) { rad2 = rad_*rad_; }
double intersect(const Ray &r) const { // returns distance, 0 if nohit
Vec op = p-r.o; // Solve t^2*d.d + 2*t*(o-p).d + (o-p).(o-p)-R^2 = 0
double t, eps=1e-4, b=op.dot(r.d), det=b*b-op.dot(op)+rad2;
if (det<0) return 0; else det=sqrt(det);
return (t=b-det)>eps ? t : ((t=b+det)>eps ? t : 0);
}
};


#define PI2 (M_PI + M_PI)


Sphere spheres[] = {//Scene: radius, position, emission, color, material
Sphere(1e5, Vec( 1e5+1,40.8,81.6), Vec(),Vec(.75,.25,.25),DIFF),//Left
Sphere(1e5, Vec(-1e5+99,40.8,81.6),Vec(),Vec(.25,.25,.75),DIFF),//Rght
Sphere(1e5, Vec(50,40.8, 1e5), Vec(),Vec(.75,.75,.75),DIFF),//Back
Sphere(1e5, Vec(50,40.8,-1e5+170), Vec(),Vec(), DIFF),//Frnt
Sphere(1e5, Vec(50, 1e5, 81.6), Vec(),Vec(.75,.75,.75),DIFF),//Botm
Sphere(1e5, Vec(50,-1e5+81.6,81.6),Vec(),Vec(.75,.75,.75),DIFF),//Top
Sphere(16.5,Vec(27,16.5,47), Vec(),Vec(1,1,1)*.999, SPEC),//Mirr
Sphere(16.5,Vec(73,16.5,78), Vec(),Vec(1,1,1)*.999, REFR),//Glas
//
Sphere(5e-3,Vec(50,81.6-36.5,81.6),Vec(4,4,4)*1e7, Vec(), DIFF),//Lite
};


/*
Sphere spheres[] = {//Scene: radius, position, emission, color, material
Sphere(1e5, Vec( 1e5+1,40.8,81.6), Vec(),Vec(.75,.25,.25),DIFF),//Left
Sphere(1e5, Vec(-1e5+99,40.8,81.6),Vec(),Vec(.25,.25,.75),DIFF),//Rght
Sphere(1e5, Vec(50,40.8, 1e5), Vec(),Vec(.75,.75,.75),DIFF),//Back
Sphere(1e5, Vec(50,40.8,-1e5+170), Vec(),Vec(), DIFF),//Frnt
Sphere(1e5, Vec(50, 1e5, 81.6), Vec(),Vec(.75,.75,.75),DIFF),//Botm
Sphere(1e5, Vec(50,-1e5+81.6,81.6),Vec(),Vec(.75,.75,.75),DIFF),//Top
Sphere(16.5,Vec(27,16.5,47), Vec(),Vec(1,1,1)*.999, SPEC),//Mirr
Sphere(16.5,Vec(73,16.5,78), Vec(),Vec(1,1,1)*.999, REFR),//Glas
Sphere(600, Vec(50,681.6-.27,81.6),Vec(12,12,12), Vec(), DIFF) //Lite
};
*/

/*
//from explicit.cpp
Sphere spheres[] = {//Scene: radius, position, emission, color, material
  Sphere(1e5, Vec( 1e5+1,40.8,81.6), Vec(),Vec(.75,.25,.25),DIFF),//Left
  Sphere(1e5, Vec(-1e5+99,40.8,81.6),Vec(),Vec(.25,.25,.75),DIFF),//Rght
  Sphere(1e5, Vec(50,40.8, 1e5),     Vec(),Vec(.75,.75,.75),DIFF),//Back
  Sphere(1e5, Vec(50,40.8,-1e5+170), Vec(),Vec(),           DIFF),//Frnt
  Sphere(1e5, Vec(50, 1e5, 81.6),    Vec(),Vec(.75,.75,.75),DIFF),//Botm
  Sphere(1e5, Vec(50,-1e5+81.6,81.6),Vec(),Vec(.75,.75,.75),DIFF),//Top
  Sphere(16.5,Vec(27,16.5,47),       Vec(),Vec(1,1,1)*.999, SPEC),//Mirr
  Sphere(16.5,Vec(73,16.5,78),       Vec(),Vec(1,1,1)*.999, REFR),//Glas
  Sphere(1.5, Vec(50,81.6-16.5,81.6),Vec(4,4,4)*100,  Vec(), DIFF),//Lite
};
*/

	
//	3/ 1* 2* 1* 2/ 2dot   => 2/ 1* 1* 1dot 1* 1/
//****** start
double molif_r = 1.; // Global mollification radius, shrinks per sample
double molif_r2 = 1.; // Global mollification radius, shrinks per sample

double mollify(Vec& l, const Vec rd, Vec& n, Vec& nl, double dist, int type){	

double cos_max=1./sqrt(1.+(molif_r2/(dist*dist) ));// Cone angle


Vec out;

if(type == REFR){
	
	bool bigger = n.dot(nl)>0;
	
	// Compute refraction vector
double nc=1,nt=1.5,nnt=bigger?nc/nt:nt/nc,ddn=rd.dot(nl),cos2t;

if((cos2t=1-nnt*nnt*(1-ddn*ddn))>0) // Refraction vector
{
out =(rd*nnt-n*(
	//(n.dot(nl)>0?1:-1)*(ddn*nnt+sqrt(cos2t))
		(bigger
			?
			(ddn*nnt+sqrt(cos2t)) : -(ddn*nnt+sqrt(cos2t))
		)
	  )
	).norm();
}
} else
{
	out = rd -(n+n)*n.dot(rd); // Reflection vector
}

/*
double solid_angle=PI2*(1.-cos_max); // Solid angle of the cone
return l.dot(out)>=cos_max ? (1./(solid_angle*l.dot(out))):0.; // Mollify
*/

//double solid_angle=PI2*(1.-cos_max); // Solid angle of the cone
return l.dot(out)>=cos_max ? (1./(PI2*(1.-cos_max)*l.dot(out))):0.; // Mollify
}
//****** end


int numSpheres = sizeof(spheres)/sizeof(Sphere);
inline double clamp(double x){ return x<0 ? 0 : x>1 ? 1 : x; }
inline int toInt(double x){ return int(pow(clamp(x),1/2.2)*255+.5); }
inline bool intersect(const Ray &r, double &t, int &id){
double n=sizeof(spheres)/sizeof(Sphere), d, inf=t=1e20;
for(int i=int(n);i--;) if((d=spheres[i].intersect(r))&&d<t){t=d;id=i;}
return t<inf;
}
Vec radiance(const Ray &r, int depth, unsigned short *Xi,int E=1){
double t; // distance to intersection
int id=0; // id of intersected object
if (!intersect(r, t, id)||depth>10) return Vec();// if miss / too deep
const Sphere &obj = spheres[id]; // the hit object
Vec x=r.o+r.d*t, n=(x-obj.p).norm(), nl=n.dot(r.d)<0?n:n*-1,f=obj.c,e;
double p = f.x>f.y && f.x>f.z ? f.x : f.y>f.z ? f.y : f.z; // max refl
if (++depth>5||!p) if (erand48(Xi)<p) f=f*(1/p); else return obj.e*E;

double tmpp;


for (int i=0; i<numSpheres; i++){ // Explicit connections
	const Sphere &s = spheres[i];
	if (s.e.x<=0 && s.e.y<=0 && s.e.z<=0) continue; // skip non-lights

	Vec sw=s.p-x,su=((fabs(sw.x)>.1?Vec(0,1):Vec(1))%sw).norm(),sv=sw%su;

	double cos_a_max = sqrt(1-s.rad2/(x-s.p).dot(x-s.p));

	double eps1 = erand48(Xi), eps2 = erand48(Xi);

	double cos_a = 1-eps1+eps1*cos_a_max;
	double sin_a = sqrt(1-cos_a*cos_a), phi = PI2*eps2;

	Vec l = su*cos(phi)*sin_a + sv*sin(phi)*sin_a + sw*cos_a;

	if (intersect(Ray(x,l.norm()), t, id) && id==i) // shadow ray
	{
		//***** start
		
		tmpp = 
		(
			(obj.refl==DIFF) ?
								(l.dot(nl)*M_1_PI) : 
								mollify(l,r.d,n,nl,t,obj.refl)
		);
		
		if(tmpp>0.0)
		e = e + f.mult(
		s.e*PI2*(1-cos_a_max)
		*
		tmpp
		
		);// BRDF
		//***** end
	}
}


if (obj.refl == DIFF){ // Ideal DIFFUSE reflection
double r1=PI2*erand48(Xi), r2=erand48(Xi), r2s=sqrt(r2);
Vec w=nl, u=((fabs(w.x)>.1?Vec(0,1):Vec(1))%w).norm(), v=w%u;
Vec d = (u*cos(r1)*r2s + v*sin(r1)*r2s + w*sqrt(1-r2)).norm();
return obj.e*E+e+f.mult(radiance(Ray(x,d),depth,Xi,0));
} else if (obj.refl == SPEC) // Ideal SPECULAR reflection
return obj.e+e+f.mult(radiance(Ray(x,r.d-n*2*n.dot(r.d)),depth,Xi));
Ray reflRay(x, r.d-n*2*n.dot(r.d)); // Ideal dielectric REFRACTION
bool into = n.dot(nl)>0; // Ray from outside going in?
double nc=1, nt=1.5, nnt=into?nc/nt:nt/nc, ddn=r.d.dot(nl), cos2t;
if ((cos2t=1-nnt*nnt*(1-ddn*ddn))<0) // Total internal reflection
return obj.e+e+f.mult(radiance(reflRay,depth,Xi));
Vec tdir = (r.d*nnt - n*((into?1:-1)*(ddn*nnt+sqrt(cos2t)))).norm();
double a=nt-nc, b=nt+nc, R0=a*a/(b*b), c = 1-(into?-ddn:tdir.dot(n));
double Re=R0+(1-R0)*c*c*c*c*c,Tr=1-Re,P=.25+.5*Re,RP=Re/P,TP=Tr/(1-P);
return obj.e+e+f.mult(depth>2 ? (erand48(Xi)<P ? // Russian roulette
radiance(reflRay,depth,Xi)*RP:radiance(Ray(x,tdir),depth,Xi)*TP) :
radiance(reflRay,depth,Xi)*Re+radiance(Ray(x,tdir),depth,Xi)*Tr);
}
int main(int argc, char *argv[]){
int w=1024,
	h=768,
	samps = argc==2 ? atoi(argv[1]) : 2048; // # samples
	
	t0 = clock();
	
Ray cam(Vec(50,52,295.6), Vec(0,-0.042612,-1).norm()); // cam pos, dir
Vec cx=Vec(w*.5135/h), cy=(cx%cam.d).norm()*.5135, r, *c=new Vec[w*h];

long ttl_samples = samps*h*w;
long sample_n = 0;

for (int s=0; s<samps; s++)
{

for(int y=0; y<h;y++){unsigned short Xi[3]={0,0,y*y*y};// Looping rows
fprintf(stderr,"\rRendering (%d spp) %5.2f%%",samps,
			(100.*
				float(sample_n)/float(ttl_samples)
			)
);
for (unsigned short x=0; x<w; x++) // Loop cols
{
//for (
int sy=0, 
i=(h-y-1)*w+x;

// sy<2; sy++) // 2x2 subpixel rows
//for (
int sx=-1;

//, r=Vec()
	//)
{ // 2x2 subpixel cols

r=Vec();

{

sx++;
if(sx==2)
{
	sx=0;
	sy++;
	if(sy==2) sy=0;
}
	
//***** start	
molif_r = 1.*pow(1.+s,-1./6); // Mollification shrinkage
molif_r2 = molif_r * molif_r; 
//***** end

double r1=2*erand48(Xi), dx=r1<1 ? sqrt(r1)-1: 1-sqrt(2-r1);
double r2=2*erand48(Xi), dy=r2<1 ? sqrt(r2)-1: 1-sqrt(2-r2);
Vec d = cx*( ( (sx+.5 + dx)/2 + x)/w - .5) +
cy*( ( (sy+.5 + dy)/2 + y)/h - .5) + cam.d;
r = r + radiance(Ray(cam.o+d*140,d.norm()),0,Xi);


sample_n++;
}
 // Camera rays are pushed ^^^^^ forward to start in interior
c[i] = c[i] + r;
}
}
}

} //for (int s=0; s<samps; s++)

t1 = clock();

printf("\nrender time %.2fs\n", float(t1-t0)/float(CLOCKS_PER_SEC));

FILE *f = fopen("image-ptreg0-opt-new1.ppm", "w"); // Write image to PPM file.
fprintf(f, "P3\n%d %d\n%d\n", w, h, 255);
for (int i=0; i<w*h; i++) {
	r = c[i]*(1./(double(samps)));
	c[i] = Vec(clamp(r.x),clamp(r.y),clamp(r.z));
fprintf(f,"%d %d %d ", toInt(c[i].x), toInt(c[i].y), toInt(c[i].z));
}
}