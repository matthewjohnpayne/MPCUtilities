// The idea here is to make a version of kepcart that works for any orbit type,
// circular, elliptical, parabolic, or hyperbolic.

// We can still use a general gm variable for the G*Mass term.
// However, the elements will now need to be something like
// q, e, incl, longnode, argperi, timeperi


#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define PI 3.14159265358979323846
#define ECL	(84381.448*(1./3600)*PI/180.) /*Obliquity of ecliptic at J2000*/

#define MAX_ITER 10000000

typedef struct {
  double x, y, z, xd, yd, zd;
} State;

typedef struct {
  double a, e, incl, longnode, argperi, meananom;
} Elements;

double machine_epsilon = 2e-15;

double sign(double x) {
    return (x > 0.0) - (x < 0.0);
}

double principal_value(double theta)
{
  theta -= 2.0*PI*floor(theta/(2.0*PI));
  return(theta);
}

// This version is older.  It takes a cartesian state and loads values into pointers to
// element variables.  It assumes an elliptical (or circular) orbit.
void keplerian(double gm, State state, 
	  double *a, double *e, double *incl, double *longnode, double *argperi, double *meananom)
{
  double rxv_x, rxv_y, rxv_z, hs, h, parameter;
  double r, vs, rdotv, rdot, ecostrueanom, esintrueanom, cosnode, sinnode;
  double rcosu, rsinu, u, trueanom, eccanom;
  double nx, ny, nz;

  double nxr_x, nxr_y, nxr_z;
  double nxrs, nxr;

  // This part is independent of the type of orbit.
  /* find direction of angular momentum vector */
  rxv_x = state.y * state.zd - state.z * state.yd;
  rxv_y = state.z * state.xd - state.x * state.zd;
  rxv_z = state.x * state.yd - state.y * state.xd;
  hs = rxv_x * rxv_x + rxv_y * rxv_y + rxv_z * rxv_z;
  h = sqrt(hs);

  r = sqrt(state.x * state.x + state.y * state.y + state.z * state.z);
  vs = state.xd * state.xd + state.yd * state.yd + state.zd * state.zd;
  /* v = sqrt(vs);  unnecessary */
  rdotv = state.x * state.xd + state.y * state.yd + state.z * state.zd;
  rdot = rdotv / r;

  // The parameter is independent of the orbit type, but its in interpretation
  // depends on the orbit type.
  parameter = hs / gm;

  *incl = acos(rxv_z / h);

  if(rxv_x!=0.0 || rxv_y!=0.0) {
    *longnode = atan2(rxv_x, -rxv_y);
  } else {
    *longnode = 0.0;
  }

  // Here a might be negative.
  double a_inv = (2.0/r - vs/gm);
  *a = 1.0/a_inv;

  ecostrueanom = parameter / r - 1.0;
  esintrueanom = rdot * h / gm;
  *e = sqrt(ecostrueanom * ecostrueanom + esintrueanom * esintrueanom);

  if(esintrueanom!=0.0 || ecostrueanom!=0.0) {
    trueanom = atan2(esintrueanom, ecostrueanom);
  } else {
    trueanom = 0.0;
  }

  cosnode = cos(*longnode);
  sinnode = sin(*longnode);

  nx = cosnode;
  ny = sinnode;
  nz = 0.0;

  /* find node vector cross r vector */
  nxr_x =   ny * state.z;  
  nxr_y = - nx * state.z;  
  nxr_z = nx * state.y - ny * state.x;
  nxrs = nxr_x*nxr_x + nxr_y*nxr_y + nxr_z*nxr_z;
  nxr = sqrt(nxrs);

  double ori = (nxr_x * rxv_x + nxr_y * rxv_y + nxr_z * rxv_z)/(nxr*h);
  //printf("%lf\n", ori);

  /* u is the argument of latitude */
  rcosu = state.x * cosnode + state.y * sinnode;
  //rsinu = (state.y * cosnode - state.x * sinnode)/cos(*incl);  // The term in the denominator might be a problem.
  //printf("%lf %lf %lf\n", rsinu, rxn, rsinu-rxn);
  rsinu = nxr*ori;

  if(rsinu!=0.0 || rcosu!=0.0) {
    u = atan2(rsinu, rcosu);
  } else {
    u = 0.0;
  }

  *argperi = principal_value(u - trueanom);
  eccanom = 2.0 * atan(sqrt((1.0 - *e)/(1.0 + *e)) * tan(trueanom/2.0));
  *meananom = eccanom - *e * sin(eccanom);

  return;
}

// This looks like an imcomplete attempt to take a cartesiant state, with an epoch
// and load values into pointers to a different set of orbital elements that can
// describe ellitical, parabolic, or hyperbolic orbits.
void universal_keplerian(double gm, State state, double tepoch,
			 double *q, double *e, double *incl, double *longnode, double *argperi, double *timeperi)
{
  double rxv_x, rxv_y, rxv_z, hs, h, parameter;
  double r, vs, rdotv, rdot, ecostrueanom, esintrueanom, cosnode, sinnode;
  double rcosu, rsinu, u, trueanom, eccanom;
  double nx, ny, nz;

  double nxr_x, nxr_y, nxr_z;
  double nxrs, nxr;

  // This part is independent of the type of orbit.
  /* find direction of angular momentum vector */
  rxv_x = state.y * state.zd - state.z * state.yd;
  rxv_y = state.z * state.xd - state.x * state.zd;
  rxv_z = state.x * state.yd - state.y * state.xd;
  hs = rxv_x * rxv_x + rxv_y * rxv_y + rxv_z * rxv_z;
  h = sqrt(hs);

  r = sqrt(state.x * state.x + state.y * state.y + state.z * state.z);
  vs = state.xd * state.xd + state.yd * state.yd + state.zd * state.zd;
  /* v = sqrt(vs);  unnecessary */
  rdotv = state.x * state.xd + state.y * state.yd + state.z * state.zd;
  rdot = rdotv / r;

  // The parameter is independent of the orbit type, but its in interpretation
  // depends on the orbit type.
  parameter = hs / gm;

  *incl = acos(rxv_z / h);

  if(rxv_x!=0.0 || rxv_y!=0.0) {
    *longnode = atan2(rxv_x, -rxv_y);
  } else {
    *longnode = 0.0;
  }

  // Here a might be negative.
  double a;
  double a_inv = (2.0/r - vs/gm);

  double temp = 1.0  +  parameter * (vs / gm  -  2.0 / r);
  if (temp < 0.0){
    *e = 0.0;
  }else{
    *e = sqrt(temp);
  }

  *q = parameter / (1.0 + *e);

  ecostrueanom = parameter / r - 1.0;
  esintrueanom = rdot * h / gm;
  *e = sqrt(ecostrueanom * ecostrueanom + esintrueanom * esintrueanom);

  if(esintrueanom!=0.0 || ecostrueanom!=0.0) {
    trueanom = atan2(esintrueanom, ecostrueanom);
  } else {
    trueanom = 0.0;
  }

  cosnode = cos(*longnode);
  sinnode = sin(*longnode);

  nx = cosnode;
  ny = sinnode;
  nz = 0.0;

  /* find node vector cross r vector */
  nxr_x =   ny * state.z;  
  nxr_y = - nx * state.z;  
  nxr_z = nx * state.y - ny * state.x;
  nxrs = nxr_x*nxr_x + nxr_y*nxr_y + nxr_z*nxr_z;
  nxr = sqrt(nxrs);

  double ori = (nxr_x * rxv_x + nxr_y * rxv_y + nxr_z * rxv_z)/(nxr*h);
  //printf("%lf\n", ori);

  /* u is the argument of latitude */
  rcosu = state.x * cosnode + state.y * sinnode;
  //rsinu = (state.y * cosnode - state.x * sinnode)/cos(*incl);  // The term in the denominator might be a problem.
  //printf("%lf %lf %lf\n", rsinu, rxn, rsinu-rxn);
  rsinu = nxr*ori;

  if(rsinu!=0.0 || rcosu!=0.0) {
    u = atan2(rsinu, rcosu);
  } else {
    u = 0.0;
  }

  *argperi = principal_value(u - trueanom);

  // At this point we have most of the elements:
  // a, e, incl, longnode, argperi
  // However, we still need timeperi

  printf("a_inv %.16lf\n", a_inv);  
  if(a_inv > 0.0){

    // This is where we left off
    // This section of the code is entered when we use e=1.0 because
    // roundoff error somehow yields an elliptic orbit for which
    // a_inv is small and a is very large.
    // Then we see that eccanom~0 because e=1
    // and meananom = 0.0
    // We end up with a problem getting timeperi because a and
    // the period are not well defined.

    // Chambers' X2EL routine has a way of approaching this
    // that seems more robust.
    
    a = 1.0 / a_inv;
    eccanom = 2.0 * atan(sqrt((1.0 - *e)/(1.0 + *e)) * tan(trueanom/2.0));
    double meananom = eccanom - *e * sin(eccanom);
    printf("%.16lf %.16lf %.16lf\n", meananom, eccanom, trueanom);
    *timeperi = tepoch - meananom * sqrt( a*a*a / gm);

    //*q = a*(1 - *e);
    
  }else if(a_inv < 0.0){

    a = 1.0 / a_inv;
    double Z = 1.0/(*e) * (1.0 - r * a_inv);
    // This is a hack.  Need a real solution.
    if(Z<1.0) Z=1.0;
    if(Z>-1.0) Z=-1.0;
    double F = log(Z + sqrt(Z*Z - 1.0))*sign(rdotv);
    *timeperi = tepoch - (*e * sinh(F) - F) * sqrt( - a*a*a / gm);
    //printf("foo F %lf Z %lf %lf %lf %lf %.16lf\n", F, Z, sign(rdotv), log(Z + sqrt(Z*Z - 1.0)), Z + sqrt(Z*Z - 1.0), Z*Z-1.0);
    //*q = a*(1 - *e);    
    
  }else{    // a_inv == 0.0

    printf("hs %lf\n", hs);
    double qp = *q;

    double tau = tan(trueanom/2.0);
    *timeperi = tepoch - (1.0/3.0 * tau * tau * tau + tau)* sqrt( 2.0 * qp * qp * qp / gm);

  }

  return;
}

void keplerians(int num, double gm, State *state, 
	  double *a, double *e, double *incl, double *longnode, double *argperi, double *meananom)
{
  int i;
  for(i=0; i<num; i++){
    keplerian(gm, state[i], &a[i], &e[i], &incl[i], &longnode[i], &argperi[i], &meananom[i]);
  }
}


void cartesian(double gm, 
	       double a, double e, double i, double longnode, double argperi, double meananom, 
	       State *state)
{
  double meanmotion, cosE, sinE, foo;
  double x, y, z, xd, yd, zd;
  double xp, yp, zp, xdp, ydp, zdp;
  double cosw, sinw, cosi, sini, cosnode, sinnode;
  double E0, E1, E2, den;
  int n;
  double principal_value(double theta);
  FILE *fout;

  /* first compute eccentric anomaly */
  /* try Steffensen's method */

  meananom = principal_value(meananom);
  n = 0;
  E0 = meananom; 
  do {
    E1 = meananom + e * sin(E0);
    E2 = meananom + e * sin(E1);

    den = E2 - 2.0*E1 + E0;
    if(fabs(den) > machine_epsilon) {
      E0 = E0 - (E1-E0)*(E1-E0)/den;
    }
    else {
      E0 = E2;
      E2 = E1;
    }
    n++;
  } while(fabs(E0-E2) > machine_epsilon && n<MAX_ITER);

  if(n>=MAX_ITER){ /* Steffensen's method failed.  Try direct iteration */
    n = 0;
    E1 = meananom; 
    do {
      E0 = E1;
      E1 = meananom + e * sin(E0);
      n++;
    } while(fabs(E0-E1) > machine_epsilon && n<MAX_ITER);

  }

  if(n>=MAX_ITER){
    fout = fopen("test_point.txt", "w");
    fprintf(fout, "%.16lf %.16lf %.16lf %.16lf %.16lf %.16lf\n", a, e, i, longnode, argperi, meananom);
    fprintf(fout, "%.16lf %.16lf %.16le %d\n", E0, E1, machine_epsilon, n);
    exit(-1);
  }

  cosE = cos(E0);
  sinE = sin(E0);

  /* compute unrotated positions and velocities */
  foo = sqrt(1.0 - e*e);
  meanmotion = sqrt(gm/(a*a*a));
  x = a * (cosE - e);
  y = foo * a * sinE;
  z = 0.0;
  xd = -a * meanmotion * sinE / (1.0 - e * cosE);
  yd = foo * a * meanmotion * cosE / (1.0 - e * cosE);
  zd = 0.0;

  /* rotate by argument of perihelion in orbit plane*/
  cosw = cos(argperi);
  sinw = sin(argperi);
  xp = x * cosw - y * sinw;
  yp = x * sinw + y * cosw;
  zp = z;
  xdp = xd * cosw - yd * sinw;
  ydp = xd * sinw + yd * cosw;
  zdp = zd;

  /* rotate by inclination about x axis */
  cosi = cos(i);
  sini = sin(i);
  x = xp;
  y = yp * cosi - zp * sini;
  z = yp * sini + zp * cosi;
  xd = xdp;
  yd = ydp * cosi - zdp * sini;
  zd = ydp * sini + zdp * cosi;

  /* rotate by longitude of node about z axis */
  cosnode = cos(longnode);
  sinnode = sin(longnode);
  state->x = x * cosnode - y * sinnode;
  state->y = x * sinnode + y * cosnode;
  state->z = z;
  state->xd = xd * cosnode - yd * sinnode;
  state->yd = xd * sinnode + yd * cosnode;
  state->zd = zd;

  return;
}

void cartesians(int num, double gm, 
		double *a, double *e, double *incl, double *longnode, double *argperi, double *meananom, 
		State *state)
{
  int i;
  for(i=0; i<num; i++){
    cartesian(gm, a[i], e[i], incl[i], longnode[i], argperi[i], meananom[i], &state[i]);
  }

}

void cartesian_vectors(int num, double gm, 
		       double *a, double *e, double *incl, double *longnode, double *argperi, double *meananom,
		       double *positions,
		       double *velocities)		       
{
  State state;
  int i;
  for(i=0; i<num; i++){
    cartesian(gm, a[i], e[i], incl[i], longnode[i], argperi[i], meananom[i], &state);
    positions[i*3+0]=state.x;
    positions[i*3+1]=state.y;
    positions[i*3+2]=state.z;        
    velocities[i*3+0]=state.xd;
    velocities[i*3+1]=state.yd;
    velocities[i*3+2]=state.zd;    
  }
}

void cartesian_elements(int num, double gm, 
			double *elements,
			double *positions,
			double *velocities)		       
{
  State state;
  int i;
  double a, e, incl, longnode, argperi, meananom;
  for(i=0; i<num; i++){

    a         = elements[6*i+0];
    e         = elements[6*i+1];
    incl      = elements[6*i+2];
    longnode  = elements[6*i+3];
    argperi   = elements[6*i+4];
    meananom  = elements[6*i+5];

    cartesian(gm, a, e, incl, longnode, argperi, meananom, &state);
    positions[i*3+0]=state.x;
    positions[i*3+1]=state.y;
    positions[i*3+2]=state.z;        
    velocities[i*3+0]=state.xd;
    velocities[i*3+1]=state.yd;
    velocities[i*3+2]=state.zd;    
  }
}

/* And transform x,y,z from ecliptic to eq */
void xyz_ec_to_eq(double x_ec, double y_ec, double z_ec,
		  double *x_eq, double *y_eq, double *z_eq)
{
  double se,ce;

  se = sin(-ECL);
  ce = cos(ECL);

  *x_eq = x_ec;
  *y_eq = ce*y_ec + se*z_ec;
  *z_eq = -se*y_ec + ce*z_ec;

  return;
}

#define KEPLER_FLAG 1
#define NEWTON_MAX 6
#define LAGCON_MAX 15
#define PI 3.14159265358979323846
#define EPS 1e-13

/* kepler step in universal variables from Danby (1988) p. 178 */

int universal_stepper(double gm, double dt, State *s0, State *s)
{
  double r0, v0s, u, alpha;
  double zeta;
  double r;
  double f, fdot, g, gdot;
  double ss;
  double g0, g1, g2, g3, g4, g5;
  int flag, kepler();

  flag = 0;
  r0 = sqrt(s0->x*s0->x + s0->y*s0->y + s0->z*s0->z);
  v0s = s0->xd*s0->xd + s0->yd*s0->yd + s0->zd*s0->zd;
  u = s0->x*s0->xd + s0->y*s0->yd + s0->z*s0->zd;
  alpha = 2.0*gm/r0 - v0s;
  zeta = gm - alpha*r0;

  /* solve kepler's equation */
  flag = kepler(gm, dt, r0, alpha, u, zeta, &r, &ss, &g0, &g1, &g2, &g3, &g4, &g5);

  f = 1.0 - (gm/r0)*g2;
  g = dt - gm*g3;
  fdot = - (gm/(r*r0))*g1;
  /*gdot = 1.0 - (gm/r)*g2;*/
  gdot = (1.0 + g*fdot)/f;

  s->x = f*s0->x + g*s0->xd;
  s->y = f*s0->y + g*s0->yd;
  s->z = f*s0->z + g*s0->zd;
  s->xd = fdot*s0->x + gdot*s0->xd;
  s->yd = fdot*s0->y + gdot*s0->yd;
  s->zd = fdot*s0->z + gdot*s0->zd;

  return(flag);
}

int kepler(double gm, double dt, double r0, double alpha, double u, double zeta,
	   double *rx, double *s, double *g0, double *g1, double *g2, double *g3, double *g4, double *g5)
{

  double x; 
  //double xx, yy, xx1, yy1, omx, h;
  //double k0x, k0y, k1x, k1y, k2x, k2y, k3y;
  double c0, c1, c2, c3, c4, c5;
  //double a, en, ch, sh, e, ec, es, dm;
  double f, fp, fpp, fppp;
  double ss, sst, dss;
  double sgn();
  void cfun(), stumpff();
  double ln;
  int nc;

  /* solve kepler's equation */
 
  // if dt=0, we're finished.
  
  // If dt>0 then ss>0 and vice versa.
  // Let's find a bracket for ss.

  double ssprev = 0.0;
  
  if(dt==0.0){
    ss = 0.0;
    ssprev = 0.0;
  }else{

    ss = 0.0;
    x = ss*ss*alpha;
    stumpff(&c0, &c1, &c2, &c3, x);
    c1 *= ss; c2 *= ss*ss; c3 *= ss*ss*ss;
    f = r0*c1 + u*c2 + gm*c3 - dt;
    double fprev;
    
    ss = 0.0;
    do{
      ssprev = ss;
      fprev = f;
      ss += dt/r0;
      x = ss*ss*alpha;
      stumpff(&c0, &c1, &c2, &c3, x);
      c1 *= ss; c2 *= ss*ss; c3 *= ss*ss*ss;
      f = r0*c1 + u*c2 + gm*c3 - dt;
    }while(f<0.0);

    do{
      ssprev = ss;
      fprev = f;
      ss *= 0.5;
      x = ss*ss*alpha;
      stumpff(&c0, &c1, &c2, &c3, x);
      c1 *= ss; c2 *= ss*ss; c3 *= ss*ss*ss;
      f = r0*c1 + u*c2 + gm*c3 - dt;
    }while(f>0.0);
    //printf("%lf %lf %lf %lf\n", ssprev, fprev, ss, f);

  }

  double ssp = ssprev;
  double ssn = ss;
  // we have a bracket now.  Let's bisect.

  nc = 0;
  while(fabs(ssp-ssn)>1e-4){
    double ss = (ssp + ssn)/2;
    x = ss*ss*alpha;
    stumpff(&c0, &c1, &c2, &c3, x);
    c1 *= ss; c2 *= ss*ss; c3 *= ss*ss*ss;
    f = r0*c1 + u*c2 + gm*c3 - dt;
    if(f<0.0){
      ssn = ss;
    }else{
      ssp = ss;
    }
    nc++;
  }

  //printf("%d ", nc);
  
  sst = ss;
  nc = 0;

  /* Newton Iteration */
  do{
    x = ss*ss*alpha;
    //cfun(x, &c0, &c1, &c2, &c3, &c4, &c5);
    stumpff(&c0, &c1, &c2, &c3, x);
    c1 *= ss; c2 *= ss*ss; c3 *= ss*ss*ss;
    
    f    = r0*c1 + u*c2 + gm*c3 - dt;
    fp   = r0*c0 + u*c1 + gm*c2;
    fpp  = zeta*c1 + u*c0;
    fppp = zeta*c0 - u*alpha*c1;
    dss = -f/fp;
    dss = -f/(fp + dss*fpp/2.0);
    dss = -f/(fp + dss*fpp/2.0 + dss*dss*fppp/6.0);
    ss += dss;
    nc++;

  }while(fabs(dss) > EPS && nc<NEWTON_MAX);

  //printf("%d ", nc);

  // Perhaps use bisection if Newton's method fails?
  
  /* Laguerre-Conway iteration if Newton fails */
  if(fabs(dss) >= EPS && nc == NEWTON_MAX){
    /*printf("Newton failed\n");*/
    ss = sst;
    ln = 5.0;
    nc = 0;
    do{
      x = ss*ss*alpha;
      //cfun(x, &c0, &c1, &c2, &c3, &c4, &c5);
      stumpff(&c0, &c1, &c2, &c3, x);
      c1 *= ss; c2 *= ss*ss; c3 *= ss*ss*ss;
    
      f   = r0*c1 + u*c2 + gm*c3 - dt;
      fp  = r0*c0 + u*c1 + gm*c2;
      fpp = zeta*c1 + u*c0;
      dss = -ln*f/(fp + sgn(fp)*
		sqrt(fabs((ln-1.0)*(ln-1.0)*fp*fp - (ln-1.0)*ln*f*fpp)));
      ss += dss;
      nc++;
    }while(fabs(dss) > EPS && nc<LAGCON_MAX);

    if(fabs(dss) >= EPS && nc == LAGCON_MAX){
      return(KEPLER_FLAG);
    }
  }
  x = ss*ss*alpha;
  cfun(x, &c0, &c1, &c2, &c3, &c4, &c5);

  c1 *= ss; c2 *= ss*ss; c3 *= ss*ss*ss; c4 *= ss*ss*ss*ss; c5 *= ss*ss*ss*ss*ss;
  *g0 = c0; *g1 = c1; *g2 = c2; *g3 = c3; *g4 = c4; *g5 = c5;
  *rx = fp;
  *s = ss;
  //printf("%lf\n", ss);
  return(0);
}

void stumpff(c0, c1, c2, c3, x)
     double *c0, *c1, *c2, *c3, x;
{
    
  int n;
  double xm;
  double d0, d1, d2, d3;

  n = 0; xm = 0.1;
  while(fabs(x)>xm){
    n++;
    x /= 4;
  }


  d2=(1-x*(1-x*(1-x*(1-x*(1-x*(1-x/182.0)/132.0)/90.0)/56.0)/30.0)/12.0)/2.0;
  d3=(1-x*(1-x*(1-x*(1-x*(1-x*(1-x/210.0)/156.0)/110.0)/72.0)/42.0)/20.0)/6.0;

  d1=1.0-x*d3;
  d0=1.0-x*d2;

  while(n>0){
    n--;
    d3=(d2+d0*d3)/4.0;
    d2=d1*d1/2.0;
    d1=d0*d1;
    d0=2.0*d0*d0-1.0;
  }
  
  *c0 = d0;
  *c1 = d1;
  *c2 = d2;
  *c3 = d3;
  return;
}
    
void cfun(double z, double *c0, double *c1, double *c2, double *c3, double *c4, double *c5)
{

/* Stumpf c-functions by argument four-folding, Mikkola */

  double C6=1.0/6.0,C132=1.0/132.0,C56=1.0/56.0;
  double C30=1.0/30.0,C24=1.0/24.0,C156=1.0/156.0;
  double C90=1.0/90.0,C110=1.0/110.0,C16=1.0/16.0,C8=1.0/8.0;
  double C72=1.0/72.0,C42=1.0/42.0,C120=1.0/120.0,U=1.0;

  double h;
  int i, k;
  h = z;
  k = 0;
  while(fabs(h)>=0.1){
    h=0.25*h;
    k++;
  }
  
  *c4 = (U-h*(U-h*(U-h*C90/(U+h*C132))*C56)*C30)*C24;
  *c5 = (U-h*(U-h*(U-h*C110/(U+h*C156))*C72)*C42)*C120;

  for(i=0;i<k;i++){
    *c3 = C6 - h*(*c5);
    *c2 = 0.5 - h*(*c4);
    *c5 = (*c5 + *c4 + (*c2)*(*c3))*C16;
    *c4 = (*c3)*(2.0 - h*(*c3))*C8;
    h=4.0*h;
  }
   
  *c3 = C6 - z*(*c5);
  *c2 = 0.5 - z*(*c4);
  *c1 = 1.0 - z*(*c3);
  *c0 = 1.0 - z*(*c2);
  
  return;
}

double sgn(double x)
{
  if(x>0.0)
    return(1.0);
  else
    return(-1.0);
}


void universal_cartesian(double gm, double q, double e, double i, double longnode, double argperi, double timeperi, double time, State *state){

  // Here's the plan:
  // 1. Determine the position and velocity vector for the object at pericenter in its orbital plane at timeperi
  // 2. Rotate those vectors by the longnode, incl, and argperi angles to orient them in space.
  // 3. Use universal function to advance those vectors to the position and velocity to time.

  double h = sqrt(gm*q*(1.0+e));
  double x, y, z, xd, yd, zd;
  double xp, yp, zp, xdp, ydp, zdp;

  int universal_kepler(double gm, double dt, State *s0, State *s);

  x = q;
  y = 0.0;
  z = 0.0;
  xd = 0.0;
  yd = h/q;
  zd = 0.0;

  /* rotate by argument of perihelion in orbit plane*/
  double cosw = cos(argperi);
  double sinw = sin(argperi);
  xp = x * cosw - y * sinw;
  yp = x * sinw + y * cosw;
  zp = z;
  xdp = xd * cosw - yd * sinw;
  ydp = xd * sinw + yd * cosw;
  zdp = zd;

  /* rotate by inclination about x axis */
  double cosi = cos(i);
  double sini = sin(i);
  x = xp;
  y = yp * cosi - zp * sini;
  z = yp * sini + zp * cosi;
  xd = xdp;
  yd = ydp * cosi - zdp * sini;
  zd = ydp * sini + zdp * cosi;

  State s0;

  /* rotate by longitude of node about z axis */
  double cosnode = cos(longnode);
  double sinnode = sin(longnode);
  s0.x = x * cosnode - y * sinnode;
  s0.y = x * sinnode + y * cosnode;
  s0.z = z;
  s0.xd = xd * cosnode - yd * sinnode;
  s0.yd = xd * sinnode + yd * cosnode;
  s0.zd = zd;

  int flag = universal_stepper(gm, time-timeperi, &s0, state);
  //printf("%d\n", flag);
  flag = 0;

  return;
  
}





void xv2el(double gm, State state, 
	  double *q_, double *a_, double *e_, double *i_, double *longnode, double *argperi, double *meananom, double *TimePeri)
{
  double r,v2;
  double hx,hy,hz, h2, h;
  double rdotv, energy; 
  double fac, inc, CapOmega, u;
  double a, e, q, face, cw, sw, w, Big_E, Big_F, Big_M, tmpf;
  double omega, TP;
  int ORBIT_TYPE; 

  // Find direction of angular momentum vector
  r  = sqrt(state.x * state.x + state.y * state.y + state.z * state.z);
  v2 = state.xd * state.xd + state.yd * state.yd + state.zd * state.zd; 
  hx = state.y * state.zd - state.z * state.yd;
  hy = state.z * state.xd - state.x * state.zd;
  hz = state.x * state.yd - state.y * state.xd;
  h2 = hx*hx + hy*hy + hz*hz;
  h = sqrt(h2);

  // r.v & Energy
  rdotv  = state.x * state.xd + state.y * state.yd + state.z * state.zd;
  energy = 0.5*v2 - gm/r;

  // Inclination
  // --- 'fac' is just being used as a temporary variable  
  fac    = hz/h;
  if( fac < -1.0) {
    inc = PI;
  }
  else if( fac < 1.0 ){ 
    inc = acos(fac);
  }
  else{
    inc = 0.0;
  }

  // CapOmega = longitude of ascending node
  // u = ... 
  fac = sqrt(hx*hx + hy*hy);        // Danby 6.15.1 would suggest that this = /h * sin(i) ? / 
  if (fac < TINY){
    CapOmega = 0.0; 
    u = atan2(state.y,state.x);
    if ( hz < 0.0) u = -u;
  }
  else{
    CapOmega = atan2(hx,-hy);
    u        = atan2( state.z / sin(inc) , state.x*cos(CapOmega) + state.y*sin(CapOmega) );
  }
  if ( CapOmega < 0.0) CapOmega = CapOmega + 2.0*PI;
  if ( u < 0.0)           u        = u + 2.0*PI;


  // ORBIT_TYPE
  if (abs(energy*r/gm) < sqrt(TINY)){
    ORBIT_TYPE = 1;                                 // 1 == PARABOLA 
  }
  else{
    a = -0.5*gm/energy;
    if ( a < 0.0 ){
      fac = -h2/(gm*a);
      if (fac > TINY){
	ORBIT_TYPE = 2;                             // 2 == HYPERBOLA
      }
      else{
	ORBIT_TYPE = 1;                             // 1 == PARABOLA
      }
    }
    else{
	ORBIT_TYPE = 0;                             // 0 == ELLIPSE
    }
  }

  // Evaluate cases
  // - Find Big_E == Eccentric Anom, Big_M = Mean Anom

  // ELLIPSE
  if (ORBIT_TYPE == 0){
    fac = 1.0 - h2/(gm*a);
    if (fac > TINY){
      e = sqrt(fac);
      q = a*(1.0-e); // Not used, added by MJP 
      face = (a - r)/(a*e);
      if (face < -1.0){
	Big_E = PI;
      }
      else if (face < 1.0){
	Big_E = acos(face);
      }
      else {      
	Big_E = 0.0;
      }
      if (rdotv < 0.0 ) Big_E = 2.0*PI - Big_E;
      fac = 1.0 - e*cos(Big_E);                   // 1 - e.cos(E) = r/a [Danby 6.15.8]
      cw = (cos(Big_E) - e)/fac;                  // 
      sw = sqrt(1.0 - e*e)*sin(Big_E)/fac;
      w = atan2(sw, cw);
      if (w < 0.0) w = w + 2.0*PI;
    }
    else{
      e     = 0.0; 
      q     = a;
      Big_E = u;
      w     = u;
    }
    Big_M = Big_E - e*sin(Big_E);
  }

  // PARABOLA
  if (ORBIT_TYPE == 1){
    a = 0.0; // Not used, added by MJP, also not "correct":  should be infinity
    q = 0.5*h2/gm;
    e = 1.0;
    fac = 2.0*q/r - 1.0;
    if (fac < -1.0){
      w = PI;
    }
    else if (fac < 1.0){
      w = acos(fac);
    }
    else{
      w = 0.0;
    }
    if (rdotv < 0.0) w = 2.0*PI - w;
    tmpf  = tan(0.5*w);
    Big_M = tmpf*(1.0 + tmpf*tmpf/3.0);
  }

  // HYPERBOLA
  if (ORBIT_TYPE == 2){
    e = sqrt(1.0 + fac);                       // 'fac' remains from setting the ORBIT_TYPE: fac = -h2/(gm*a) 
    q = a*(1-e);                               // Not used, added by MJP 
    tmpf = (a - r)/(a*e);
    if (tmpf < 1.0) tmpf = 1.0;
    Big_F = log(tmpf + sqrt(tmpf*tmpf - 1.0));
    if (rdotv < 0.0) Big_F = -Big_F;
    fac = e*cosh(Big_F) - 1.0;
    cw = (e - cosh(Big_F))/fac;
    sw = sqrt(e*e - 1.0)*sinh(Big_F)/fac;
    w = atan2(sw, cw);
    if (w < 0.0) w = w + 2.0*PI;
    Big_M = e*sinh(Big_F) - Big_F;
  }

  // omega = argument of pericenter
  omega = u - w; 
  if (omega < 0.0) omega = omega + 2.0*PI;

  // time of pericenter
  TP = 1.0;
  /*
  if (ORBIT_TYPE == 0){ // ELLIPSE
    
  }
  if (ORBIT_TYPE == 1){ // PARABOLA

  }
  if (ORBIT_TYPE == 2){ // HYPERBOLA
    double tau = tan(trueanom/2.0);
    TP = tepoch - (1.0/3.0 * tau * tau * tau + tau)* sqrt( 2.0 * q * q * q / gm);
  }
  */

  // Populate output 
  // --- Doing it here because I'm trying to leave a clean comparison to the fortran code 
  *q_        = q;
  *a_        = a;
  *e_        = e;
  *i_        = inc;
  *longnode  = CapOmega;
  *argperi   = omega; 
  *meananom  = Big_M;
  *TimePeri  = TP;


  return;
}

