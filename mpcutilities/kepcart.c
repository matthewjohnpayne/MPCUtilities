#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/*
# --------------------------------------------------------------
# Oct 2018
# Payne
#
# conversion funcs for keplerian <-> cartesian coordinates
#
# straight port/copy/development of Holman's routines
# - Note some name-changes & duplication to make the purpose of 
#   the routines more obvious
#
# I suspect that some of these routines are inaccurate / die for 
# cases of hyperbolic orbits: INVESTIGATE & FIX
#  
# --------------------------------------------------------------
*/


// --- Given time & energy we should move these constants ------
// --- and classes to a more central location ------------------
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
// --------------------------------------------------------------




void cart2kep(double gm, State state, double *a, double *e, double *incl, double *longnode, double *argperi, double *meananom)
{
    /*
     Converts a "state" of cartesian coordinates to keplerian coordinates
     Keplerians elements are: (a, e, incl, longnode, argperi, meananom)
     
     Parameters
     ----------
     GM          : float,
        Gravitational Constant
     cartState   : "CartState" Object-type as defined above.
        Assumes HELIOCENTRIC ECLIPTIC CARTESIAN initial conditions
     
     Evaluates
     -------
     a, e, incl, longnode, argperi, meananom : pointers to doubles
     
     Examples
     --------
     >>> ...
    */

    double rxv_x, rxv_y, rxv_z, hs, h, parameter;
    double r, vs, rdotv, rdot, ecostrueanom, esintrueanom, cosnode, sinnode;
    double rcosu, rsinu, u, trueanom, eccanom;

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
    parameter = hs / gm;

    *incl = acos(rxv_z / h);

    if(rxv_x!=0.0 || rxv_y!=0.0) {
    *longnode = atan2(rxv_x, -rxv_y);
    } else {
    *longnode = 0.0;
    }

    *a = 1.0 / (2.0 / r - vs / gm);

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

    /* u is the argument of latitude */
    rcosu = state.x * cosnode + state.y * sinnode;
    rsinu = (state.y * cosnode - state.x * sinnode)/cos(*incl);

    if(rsinu!=0.0 || rcosu!=0.0) {
    u = atan2(rsinu, rcosu);
    } else {
    u = 0.0;
    }

    *argperi = u - trueanom;

    eccanom = 2.0 * atan(sqrt((1.0 - *e)/(1.0 + *e)) * tan(trueanom/2.0));
    *meananom = eccanom - *e * sin(eccanom);

    return;
}


void keplerian(double gm, State state, double *a, double *e, double *incl, double *longnode, double *argperi, double *meananom)
{
    /*
     Exactly the same as "cart2kep" (see above)
     Provided for backward compatibility with Holman's legacy code
     */
    cart2kep( gm,  state,  a,  e,  incl, longnode,  argperi,  meananom);
}


void cart2kep_array(int num, double gm, State *state,
                    double *a, double *e, double *incl, double *longnode, double *argperi, double *meananom)
{
    /*
     Converts arrays of cartesian coordinates (stateArrays) to arrays of keplerian coordinates
     Keplerians elements are: (a, e, incl, longnode, argperi, meananom)
     
     *** MAKES MULTIPLE CALLS TO "cart2kep" ABOVE ***
     
     Parameters
     ----------
     num         : int
     Length of the input stateArray
     GM          : float,
     Gravitational Constant
     cartState   : "CartStateArray" Object-type as defined in MPCFormat.
     Assumes HELIOCENTRIC ECLIPTIC CARTESIAN initial conditions
     
     Evaluates
     -------
     a, e, incl, longnode, argperi, meananom : arrays
     
     Examples
     --------
     >>> ...
     
     */
    
    int i;
    for(i=0; i<num; i++){
        cart2kep(gm, state[i], &a[i], &e[i], &incl[i], &longnode[i], &argperi[i], &meananom[i]);
    }
}

void keplerians(int num, double gm, State *state,
                double *a, double *e, double *incl, double *longnode, double *argperi, double *meananom)
{
    /*
     Exactly the same as "cart2kep_array" (see above)
     Provided for backward compatibility with Holman's legacy code
     */
    cart2kep_array( num,  gm,  state,  a,  e,  incl,  longnode,  argperi,  meananom);
}



double principal_value(double theta)
{
    /*
     Forces theta
     
     Parameters
     ----------
     theta         : double
     radians
     
     Evaluates
     -------
     theta          : double
     
     Examples
     --------
     >>> ...
     
     */
    theta -= 2.0*PI*floor(theta/(2.0*PI));
    return(theta);
}

void kep2cartState(double gm, double a, double e, double i, double longnode, double argperi, double meananom, State *state)
{
    /*
     Converts keplerian coordinates to a cartesian state
     Keplerians elements are: (a, e, incl, longnode, argperi, meananom)
     Assumes HELIOCENTRIC ECLIPTIC KEPLERIAN initial conditions
     
     
     Parameters
     ----------
     GM          :   float
     Gravity
     a           :   float
     Semi-major axis
     e           :   float
     Eccentricity
     incl        :   float
     Inclination
     longnode    :   float
     Longitude of ascending node
     argperi     :   float
     Argument of pericenter
     meananom    :   float
     Mean anomaly
     
     
     Evaluates
     -------
     cartState   :   "CartState" Object-type as defined above.
     
     Examples
     --------
     >>> ...
     
     */
    
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



void cartesian(double gm, double a, double e, double i, double longnode, double argperi, double meananom, State *state)
{
    /*
     Exactly the same as "kep2cartState" (see above)
     Provided for backward compatibility with Holman's legacy code
     */
    kep2cartState( gm, a,  e,  i,  longnode,  argperi,  meananom,  state);
}



/*
 void cartesian_prev(double gm,
 double a, double e, double incl, double longnode, double argperi, double meananom,
 State *state)
 {
 double meanmotion, cosE, sinE, foo;
 double x, y, z, xd, yd, zd;
 double xp, yp, zp, xdp, ydp, zdp;
 double cosw, sinw, cosi, sini, cosnode, sinnode;
 double E0, E1, E2, den;
 
 // first compute eccentric anomaly
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
 } while(fabs(E0-E2) > machine_epsilon);
 
 cosE = cos(E0);
 sinE = sin(E0);
 
 // compute unrotated positions and velocities
 foo = sqrt(1.0 - e*e);
 meanmotion = sqrt(gm/(a*a*a));
 x = a * (cosE - e);
 y = foo * a * sinE;
 z = 0.0;
 xd = -a * meanmotion * sinE / (1.0 - e * cosE);
 yd = foo * a * meanmotion * cosE / (1.0 - e * cosE);
 zd = 0.0;
 
 // rotate by argument of perihelion in orbit plane
 cosw = cos(argperi);
 sinw = sin(argperi);
 xp = x * cosw - y * sinw;
 yp = x * sinw + y * cosw;
 zp = z;
 xdp = xd * cosw - yd * sinw;
 ydp = xd * sinw + yd * cosw;
 zdp = zd;
 
 // rotate by inclination about x axis
 cosi = cos(incl);
 sini = sin(incl);
 x = xp;
 y = yp * cosi - zp * sini;
 z = yp * sini + zp * cosi;
 xd = xdp;
 yd = ydp * cosi - zdp * sini;
 zd = ydp * sini + zdp * cosi;
 
 // rotate by longitude of node about z axis
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
 */


void kep2cartStateArray(int num, double gm,
                        double *a, double *e, double *incl, double *longnode, double *argperi, double *meananom,
                        State *state)
{
    /*
     Converts arrays of keplerian coordinates to a cartesian state array
     Keplerians elements are: (a, e, incl, longnode, argperi, meananom)
     Assumes HELIOCENTRIC ECLIPTIC KEPLERIAN initial conditions
     
     *** MAKES MULTIPLE CALLS TO "kep2cartState" ABOVE ***
     
     Parameters
     ----------
     GM          :   float
     Gravity
     a           :   array
     Semi-major axis
     e           :   array
     Eccentricity
     incl        :   array
     Inclination
     longnode    :   array
     Longitude of ascending node
     argperi     :   array
     Argument of pericenter
     meananom    :   array
     Mean anomaly
     
     
     Returns
     -------
     cartStateArray      :   "CartState" Object-type as defined above.
     
     Examples
     --------
     >>> ...
     
     */
    
    int i;
    for(i=0; i<num; i++){
        kep2cartState(gm, a[i], e[i], incl[i], longnode[i], argperi[i], meananom[i], &state[i]);
    }
}

void cartesians(int num, double gm,
                double *a, double *e, double *incl, double *longnode, double *argperi, double *meananom,
                State *state)
{
    /*
     Exactly the same as "kep2cartStateArray" (see above)
     Provided for backward compatibility with Holman's legacy code
     */
    kep2cartStateArray( num,  gm, a,  e,  incl,  longnode,  argperi,  meananom, state);
}







void kep2cartPV(int num, double gm,
                double *a, double *e, double *incl, double *longnode, double *argperi, double *meananom,
                double *positions,
                double *velocities)
{
    /*
     Converts arrays of keplerian coordinates to arrays of cartesian positions and velocities
     Keplerians elements are: (a, e, incl, longnode, argperi, meananom)
     Assumes HELIOCENTRIC ECLIPTIC KEPLERIAN initial conditions
     
     *** MAKES MULTIPLE CALLS TO "kep2cartState" ABOVE ***
     
     Parameters
     ----------
     num         :   int
     Length of arrays of keplerian coordinates (a, e, ...)
     GM          :   float
     Gravity
     a           :   array
     Semi-major axis
     e           :   array
     Eccentricity
     incl        :   array
     Inclination
     longnode    :   array
     Longitude of ascending node
     argperi     :   array
     Argument of pericenter
     meananom    :   array
     Mean anomaly
     
     
     Returns
     -------
     pos_arr, vel_arr    :
     ** DESCRIBE THESE ***
     
     
     Examples
     --------
     >>> ...
     
     */
    
    State state;
    int i;
    for(i=0; i<num; i++){
        kep2cartState(gm, a[i], e[i], incl[i], longnode[i], argperi[i], meananom[i], &state);
        
        positions[i*3+0]=state.x;
        positions[i*3+1]=state.y;
        positions[i*3+2]=state.z;
        
        velocities[i*3+0]=state.xd;
        velocities[i*3+1]=state.yd;
        velocities[i*3+2]=state.zd;
    }
}

void cartesian_vectors(int num, double gm,
                       double *a, double *e, double *incl, double *longnode, double *argperi, double *meananom,
                       double *positions,
                       double *velocities)
{
    /*
     Exactly the same as "kep2cartPV" (see above)
     Provided for backward compatibility with Holman's legacy code
     */
    kep2cartPV( num,  gm,
                a,  e,  incl,  longnode,  argperi, meananom,
                positions,
                velocities);
}






void kepState2cartPV(int num, double gm,
                     Elements *elements,
                     double *positions,
                     double *velocities)
{
    /*
     Converts arrays of keplerian element-objects to arrays of cartesian positions and velocities
     Assumes HELIOCENTRIC ECLIPTIC KEPLERIAN initial conditions
     
     *** MAKES MULTIPLE CALLS TO "kep2cartState" ABOVE ***
     
     Parameters
     ----------
     GM              :   float
     Gravity
     elementsArray   :   "elementsArray" Object-type as defined in MPCFormat.
     
     
     Returns
     -------
     pos_arr, vel_arr    :
     ** DESCRIBE THESE ***
     
     
     Examples
     --------
     >>> ...
     
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
        
        kep2cartState(gm, a, e, incl, longnode, argperi, meananom, &state);
        
        positions[i*3+0]=state.x;
        positions[i*3+1]=state.y;
        positions[i*3+2]=state.z;
        velocities[i*3+0]=state.xd;
        velocities[i*3+1]=state.yd;
        velocities[i*3+2]=state.zd;
    }
     */
    State state;
    int i;
    for(i=0; i<num; i++){
        kep2cartState(gm, elements[i].a, elements[i].e, elements[i].incl, elements[i].longnode, elements[i].argperi, elements[i].meananom, &state);
        
        positions[i*3+0]=state.x;
        positions[i*3+1]=state.y;
        positions[i*3+2]=state.z;
        
        velocities[i*3+0]=state.xd;
        velocities[i*3+1]=state.yd;
        velocities[i*3+2]=state.zd;
    }
}

void cartesian_elements(int num, double gm,
                     Elements *elements,
                     double *positions,
                     double *velocities)
{
    /*
     Exactly the same as "kepState2cartPV" (see above)
     Provided for backward compatibility with Holman's legacy code
     */
    kepState2cartPV( num,  gm,
                     elements,
                     positions,
                     velocities);
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



