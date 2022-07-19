/* Header file for X-Ray structure factor
 * Author		: Yin Zhang
 * Version		: 0.99
 * Date			: March 2nd, 2018 */
#include	<stdio.h>
#include	<stdlib.h>
#include	<math.h>
#include	<string.h>
#include	<malloc.h>
#include	<stdarg.h>

typedef struct{
	double	x, y, z;		// coords of this atom
	double	f;				// atom form factor
}atom;

static	double	Lx, Ly, Lz;

/* Read cfgfile*/
atom * read_cfg(char * cfgfile, int *N, double *qdx, double *qdy,double *f2sum);

/*caculcate the S(q)*/
double	Sq(atom *p, int Natom, double qx, double qy, double qz, double f2sum);
