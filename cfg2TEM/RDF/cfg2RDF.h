/* Header file for cfg2RDF.c
 * Author	: Yin Zhang
 * Version	: 1.0
 * Date		: Dec 25th 2016 */

static	int		N0 = 10;
static	double	dr = 5;		
static	double	H[4][4] = {0};
static	int		Natom = 0;
static	char	cfgfile[50], mass[50], element[50], output[50];
static	double	xlo, xhi, ylo, yhi;
static	double	** gray;
static	int		rmax_dr;
static	int		RDF_atom[10000];
static	int		N_RDF = 0;

typedef struct ATOM{
	int		ID;
	double	x, y, z;
}atom;

int		read_input();

/* Read CFG file */
atom *	Read_CFG();

/* Write JPEG image */
void	WriteJPEG();
