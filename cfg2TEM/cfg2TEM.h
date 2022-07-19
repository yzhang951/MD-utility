/* Header file for cfg2TEM.c
 * Author	: Yin Zhang
 * Version	: 1.0
 * Date		: Dec 25th 2016 */

static	int		MAX_RGB = 250;
static	int		w = 5;		// Standard deviation for Gaussian distribution
static	int		ROW = 1920;
static	int		COL = 1080;
static	double	H[4][4] = {0};
static	int		Natom = 0;
static	char	cfgfile[50], mass[50], element[50], output[50];
static	double	xlo, xhi, ylo, yhi, zlo, zhi;
static	double	** gray;
static	int		pro_dir;	// projection direction, 1=x, 2=y, 3=z;
static	double	b;

typedef struct ATOM{
	int		ID;
	double	x, y, z;
}atom;

int		read_input();

/* Read CFG file */
atom *	Read_CFG();

/* Write JPEG image */
void	WriteJPEG();
