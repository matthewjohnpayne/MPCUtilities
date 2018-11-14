/* There must be a better version of this around somewhere in the Stark directory */

#include <string.h>
#include <stdio.h>
#include <math.h>

#define KEPLER_FLAG 1

typedef struct {
  double x, y, z, xd, yd, zd;
} State;

extern double machine_epsilon;

#define NEWTON_MAX 6
#define LAGCON_MAX 15
#define PI 3.14159265358979323846
#define EPS 1e-13

/* kepler step in universal variables from Danby (1988) p. 178 */

int universal_stepper(double gm, double dt, State *s0, State *s)
{
  double r0, v0s, u, alpha;
  double zeta, zetapr;
  double r, x0dx, x0dv, v0dx, v0dv;
  double r0pr, rpr, upr, alphapr;
  double dx, dy, dz, dxd, dyd, dzd;
  double f, fdot, g, gdot;
  double fpr, fdotpr, gpr, gdotpr;
  double ss, sspr;
  double g0, g1, g2, g3, g4, g5;
  double g1a, g2a, g3a, g1pr, g2pr, g3pr;
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

int kepler(gm, dt, r0, alpha, u, zeta, rx, s, g0, g1, g2, g3, g4, g5)
     double gm, dt, r0, alpha, u, zeta, *rx, *s, *g0, *g1, *g2, *g3, *g4, *g5;
{

  double x, y, z; 
  double xx, yy, xx1, yy1, omx, h;
  double k0x, k0y, k1x, k1y, k2x, k2y, k3y;
  double c0, c1, c2, c3, c4, c5;
  double a, en, ch, sh, e, ec, es, dm, sigma;
  double f, fp, fpp, fppp, fdt;
  double ss, sst, dss;
  double sgn();
  void cfun(), stumpff();
  double ln;
  int nc;

  /* solve kepler's equation */
 
  /* determine initial guesses */
  if(fabs(dt/r0) > 0.2){
    if(alpha<=0.0){
      /* hyperbolic motion */
      a = gm/alpha;
      en = sqrt(-gm/(a*a*a));
      ch = 1.0 - r0/a;
      sh = u/sqrt(-a*gm);
      e = sqrt(ch*ch - sh*sh);
      dm = en*dt;
      if(dm<0.0){
	ss = -log((-2.0*dm + 1.8*e)/(ch - sh))/sqrt(-alpha);
      }else{
	ss = log((2.0*dm + 1.8*e)/(ch + sh))/sqrt(-alpha);
      }
    }else{
      /* elliptic motion */
      a = gm/alpha;
      en = sqrt(gm/(a*a*a));
      ec = 1.0 - r0/a;
      es = u/(en*a*a);
      e = sqrt(ec*ec + es*es);
      /*dt -= trunc(en*dt/(2.0*PI))*2.0*PI/en;*/
      dt -= ((int)(en*dt/(2.0*PI)))*2.0*PI/en;
      y = en*dt - es;

      /* RK4 step for initial guess */

      xx = ec;
      yy = es;
      h = en*dt;
      omx = h/(1.0 - xx);
      k0x = - yy * omx;
      k0y = xx * omx;
      xx1 = xx + k0x/2.0;
      yy1 = yy + k0y/2.0;
      omx = h/(1.0 - xx1);
      k1x = - yy1 * omx;
      k1y = xx1 * omx;
      xx1 = xx + k1x/2.0;
      yy1 = yy + k1y/2.0;
      omx = h/(1.0 - xx1);
      k2x = - yy1 * omx;
      k2y = xx1 * omx;
      xx1 = xx + k2x;
      yy1 = yy + k2y;
      omx = h/(1.0 - xx1);
      k3y = xx1 * omx;
      yy += (k0y + 2.0*(k1y + k2y) + k3y)/6.0;
      
      x = y + yy;

      /* Danby Guess */
      /*
      sigma = sgn((es*cos(y) + ec*sin(y)));
      x = y + sigma*0.85*e;
      printf("DAN: %lf %lf\n", y, x/sqrt(alpha)); 
      */
      ss = x/sqrt(alpha);
      
    }
      
  }else{
    /* close to parabolic */
    ss = dt/r0 - (dt*dt*u)/(2.0*r0*r0*r0);
  }

  sst = ss;
  nc = 0;

  /* Newton Iteration */
  do{
    x = ss*ss*alpha;
    cfun(x, &c0, &c1, &c2, &c3, &c4, &c5);
    /*stumpff(&c0, &c1, &c2, &c3, x);*/
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

    /*    
    fdt = f/dt; 
    }while(fabs(fdt*fdt) > EPS*EPS && nc<NEWTON_MAX);
    */

  // Perhaps use bisection if Newton's method fails?
  
  /* Laguerre-Conway iteration if Newton fails */
  if(fabs(dss) >= EPS && nc == NEWTON_MAX){
    /*printf("Newton failed\n");*/
    ss = sst;
    ln = 5.0;
    nc = 0;
    do{
      x = ss*ss*alpha;
      cfun(x, &c0, &c1, &c2, &c3, &c4, &c5);
      /*stumpff(&c0, &c1, &c2, &c3, x);*/
      c1 *= ss; c2 *= ss*ss; c3 *= ss*ss*ss;
    
      f   = r0*c1 + u*c2 + gm*c3 - dt;
      fp  = r0*c0 + u*c1 + gm*c2;
      fpp = zeta*c1 + u*c0;
      dss = -ln*f/(fp + sgn(fp)*
		sqrt(fabs((ln-1.0)*(ln-1.0)*fp*fp - (ln-1.0)*ln*f*fpp)));
      ss += dss;
      nc++;
    }while(fabs(dss) > EPS && nc<LAGCON_MAX);

      /*
      fdt = f/dt; 
      }while(fabs(fdt*fdt) > EPS*EPS && nc<LAGCON_MAX);
      */

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

