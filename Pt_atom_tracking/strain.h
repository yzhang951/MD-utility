/* Header file for particle tracking
 * Author	:	Yin	Zhang
 * Version	:	1.3
 * Date		:	Oct 18th 2016 */

#define			MAX_image	20
#define			MAX_area	10

static	int		w = 5;			// w is a user defined parameter in eq (1)
static	int		w2= 3;
static	int 	lamda = 1;		// lamda is the half width of Gaussian noise
static	double 	r = 0.7;		// filter background
static	int		MAX = 900000;	// Maximum number of particle
static	double	L = 10;		// Maximum distance a point is allowed 
static	int		ROW = 1920;
static	int		COL = 1080;
static	int 	offset = 3;
static	double	mul = 16;	// factor that dinstinguish center and edge
static	char	file[10][50];
static	int		r1, r2, r3, R1, R2, R3, G1, G2, G3, B1, B2, B3, line_width;
static	double	cutoff;
static	int		Nimage;
static	int 	Narea;
static	int		Num_offset = 0;
static	int		x_o;
static	int 	y_o;
static	double	xlo[MAX_area], xhi[MAX_area], ylo[MAX_area], yhi[MAX_area];
static	double	XLO, XHI, YLO, YHI;

typedef struct{
	int		ID;				// ID number of this particle
	int		pair_prev;		// corresponding particle in previous photo
	int		pair_next;		// corresponding particle in next photo
							// pair = 0 means dummy pair particle
	unsigned char	R, G, B;
	double	A;				// Intensity
	double	x, y;			// x = row coordinate, y = col coordinate
	int		row, col;
	double	s1, s2;			// reduced coordinate for (x,y)
	double	m0, m2;
	bool	particle;

	int		neigh_list[50];	// start from 0
	int		neigh_num;
	double	distance_list[50];
	double	M[3][3];
	double	M_prime[3][3];
	double	detM;
	bool	t_l;
	int		z;
	double	z_;
}point;

typedef struct{
	double	a, b, c;
	/* equation of :  ax + by + c = 0, a^2+b^2 = 1 */

	double	xy_, x_, y_, x2_, y2_;
	double	xmin, xmax;
	int		id;
	int		num_p;
}LINE;

static 	LINE	p_line[1000];
static	int		n_line = 1;

int sgn(int	x){
	if(x>=0) return 1;

	if(x<0) return -1;
}

int	listcpy(point * p1, point * p2, int n){
	int i, j, k;
	for(i=1;i<=n;i++){
		p1[i].ID = p2[i].ID;
		p1[i].pair_prev = p2[i].pair_prev;
		p1[i].pair_next = p2[i].pair_next;
		p1[i].A = p2[i].A;
		p1[i].x = p2[i].x;
		p1[i].y = p2[i].y;
		p1[i].row = p2[i].row;
		p1[i].col = p2[i].col;
		p1[i].s1= p2[i].s1;
		p1[i].s2= p2[i].s2;
		p1[i].m0= p2[i].m0;
		p1[i].m2= p2[i].m2;
		p1[i].particle = p2[i].particle;
		p1[i].neigh_num = p2[i].neigh_num;
		p1[i].detM = p2[i].detM;
		for(j=0;j<p2[i].neigh_num;j++)
			p1[i].neigh_list[j] = p2[i].neigh_list[j];
		for(j=1;j<3;j++)
			for(k=1;k<3;k++)
				p1[i].M[j][k] = p2[i].M[j][k];
	}
	return 0;
}

int	bubble_sort(point *p){
	int		i, n, tmp_i;
	double	tmp_d;
	bool	flag = false;

	n = p->neigh_num;
	
	do{
		flag = false;
		for(i=1;i<n;i++){
			if(p->distance_list[i-1] > p->distance_list[i]){
				tmp_d = p->distance_list[i-1];
				p->distance_list[i-1] = p->distance_list[i];
				p->distance_list[i] = tmp_d;
				tmp_i = p->neigh_list[i-1];
				p->neigh_list[i-1] = p->neigh_list[i];
				p->neigh_list[i] = tmp_i;
				flag = true;
			}	
		}
	}while(flag);
	
	return 0;
}

int	build_neigh_for_new_atom(point *list, int i){
	// new atoms should be placed on the last position of the list
	double distance;
	int	j;

	list[i].M[1][1] = 0;
	list[i].M[1][2] = 0;
	list[i].M[2][2] = 0;
	list[i].M[2][1] = 0;

	for(j=1;j<i;j++){


	   	distance = (list[i].x - list[j].x) * (list[i].x - list[j].x);
       
  		distance+= (list[i].y - list[j].y) * (list[i].y - list[j].y);

		distance = sqrt(distance);

		if(distance < cutoff){
			list[i].neigh_list[list[i].neigh_num] = j;
			
			list[j].neigh_list[list[j].neigh_num] = i;
																               
		  	list[i].distance_list[list[i].neigh_num] = distance;
		
			list[j].distance_list[list[j].neigh_num] = distance;

			list[i].neigh_num = list[i].neigh_num + 1;

			list[j].neigh_num = list[j].neigh_num + 1;

			bubble_sort(&list[i]);

			bubble_sort(&list[j]);



            list[i].M[1][1] += (list[i].x - list[j].x) * (list[i].x - list[j].x);
			list[i].M[1][2] += (list[i].x - list[j].x) * (list[i].y - list[j].y);
			
			list[i].M[2][1] = list[i].M[1][2];								           
			list[i].M[2][2] += (list[i].y - list[j].y) * (list[i].y - list[j].y);

            list[j].M[1][1] += (list[i].x - list[j].x) * (list[i].x - list[j].x);
			list[j].M[1][2] += (list[i].x - list[j].x) * (list[i].y - list[j].y);
			
			list[j].M[2][1] = list[i].M[1][2];								           
			list[j].M[2][2] += (list[i].y - list[j].y) * (list[i].y - list[j].y);


		}
	}


	return 0;
}


int	read_input();

void abort_(const char * s, ...);

int	strain(double M0[3][3], double M1[3][3], double * exx, double * eyy, double * exy);

double ** LoadPNG (char * filename, unsigned long * row, unsigned long * col);
/* Load PNG file to grayscale matrix */
point ** LoadJPEG(char * filename, unsigned long * row, unsigned long * col);
/* Load JPEG file to grayscale matrix */

void	WriteJPEG(char * filename, point ** p, unsigned long row, unsigned long col);

int		Draw_circle(point **p, unsigned long row, unsigned long col, double s1, double s2, int Raidus, unsigned char R, unsigned char G, unsigned char B, int line_width);

int		Draw_line(point **p, unsigned long row, unsigned long col, double s1_1, double s1_2, double s2_1, double s2_2, unsigned char R, unsigned char G, unsigned char B, int line_width);
/* Draw a circle in image, s1 s2 are center's reduced coordinate, R is radius */

int Image_restoration(point ** p, unsigned long row, unsigned long col);
/* Image restoration part in 
 * < Feature point tracking and trajectory analysis for video imaging in cell 
 * biology >, equation 1-6													*/

bool	estimate_point(point ** p, unsigned long row, unsigned long col, long x, long y);

bool	refine_point(point ** p, unsigned long row, unsigned long col, long x, long y);

point *	sort_point(point ** p, unsigned long row, unsigned long col, int * n);
/* Return a list of particle point, start from 1, n = number of particle */


double	cost(point * list1, int i, point * list2, int j);
/* Return the cost function for topological matrix G */


int	** link_initialization(point * list1, int n1, point * list2, int n2);
/* First step for linking particles, 2.3.1 in
 *  <http://dx.doi.org/10.1016/j.jsb.2005.06.002 >	*/

double	link_optimization(point * list1, int n1, point * list2, int n2, int **G);
/* Optimize the total cost function, 2.3.2	*/

int	offset_translation(point * list1, int n1, unsigned long row2, unsigned long col2, point * list2, int n2, int **G);

int	find_lost_atoms(point ** p[MAX_image], point ** plist,int *n1, int *n2, int *n3, int NUM_list2);
/* using three images to find lost atoms, n = number of list2*/

point ** final_result(point * list1, int n1, point * list2, int n2, int **G1, point * list3, int n3, int **G2, int row, int col);
/* Draw the final result in to image */

int	mark_z(point *p, int n);

int	WriteCFG(point **p, int N, double X, double Y, double Z, double mass, char *ele);
