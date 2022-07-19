/* Read JPEG image 
 * Author	:	Yin Zhang
 * Version	:	1.3
 * Date		:	Oct 18th	2016 */

#include <stdio.h>
#include <jpeglib.h>
#include <jerror.h>
#include <string.h>
#include <stdarg.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "strain.h"

#define	PNG_DEBUG	3
#include <png.h>

int	read_input(){
	char	line[255];
	int		count = 0;
	int		i = 0;
	int		j = 0;

	while(fgets(line, 255, stdin)){
//		printf("%s", line);

		if(line[0] == '#')	continue;

		if(count == 0)
			sscanf(line, "%d%d%lf%d%lf%d", &w, &w2, &cutoff, &lamda, &r, &offset);
		if(count == 1)
			sscanf(line,"%d%lf%d", &MAX, &L, &mul);
		if(count == 2)
			sscanf(line,"%d%d", &ROW, &COL);
		if(count == 3){
			if(i==0){
				sscanf(line,"%d", &Nimage);
				Nimage = Nimage + 1;
				if(Nimage >= MAX_image)
					abort_("Number of image exceed limit!\n");
				i++;
			}
			else if(i<Nimage){
				sscanf(line,"%s", file[i]);
				i++;
			}
			else{
				sscanf(line,"%s", file[0]);
				count++;
			}
			count--;
		}
		if(count == 4)
			sscanf(line,"%d%d%d%d%d%d%d%d%d%d%d%d%d", &r1, &R1, &G1, &B1, &r2, &R2, &G2, &B2, &r3, &R3, &G3, &B3, &line_width);

		if(count == 5){
			if(j==0){
				sscanf(line,"%d", &Narea);

				if(Narea >= MAX_area)
					abort_("Number of area exceed limit!\n");
				j++;
			}
			else if(j<=Narea){

				sscanf(line,"%lf%lf%lf%lf", &xlo[j-1], &xhi[j-1], &ylo[j-1], &yhi[j-1]);

				j++;
			}
			else{

				count++;
			}
			count--;
		}

		if(count == 6){
			sscanf(line,"%d%d%d", &Num_offset, x_o, y_o);
		}

		count ++;

	}

	return 0;
}

void	abort_(const char * s, ...)
{
	va_list	args;
	va_start(args, s);
	vfprintf(stderr, s, args);
	fprintf(stderr, "\n");
	va_end(args);
	abort();
}

/****************************************************************************
 * 				Calculate the strain using deformation gradient				*
 ****************************************************************************/
int	strain(double M0[3][3], double M1[3][3], double * exx, double * eyy, double *exy){

	double	Minv[3][3], F[3][3] = {0}, C[3][3] = {0}, det;
	int		i, j, k;

	det = M0[1][1] * M0[2][2] - M0[1][2] * M0[2][1];

	Minv[1][1] = M0[2][2] / det;
	Minv[1][2] = -M0[1][2] /det;
	Minv[2][1] = Minv[1][2];
	Minv[2][2] = M0[1][1] / det;
	
	for(i=1;i<3;i++)
		for(j=1;j<3;j++)
			for(k=1;k<3;k++)
				F[i][j] += Minv[i][k] * M1[k][j];

	/* Left Cauchy Tensor, B = F * F' 			*/
	for(i=1;i<3;i++)
		for(j=1;j<3;j++)
			for(k=1;k<3;k++)
				C[i][j] += F[i][k]*F[j][k];

	*exx = (C[1][1] - 1) / 2;
	*eyy = (C[2][2] - 1) / 2;
	*exy = C[1][2] / 2;

	return 0;
}

/****************************************************************************
 * 							Loading	PNG	file								*/

double **	LoadPNG(char * filename, unsigned long * row, unsigned long * col)
/****************************************************************************/
{
	unsigned long	x, y, i, j;
	double	**		gray;

	png_byte	color_type;
	png_byte	bit_depth;
	png_structp	png_ptr;
	png_infop	info_ptr;
	int			number_of_passes;
	png_bytep *	row_pointers;
	
	unsigned char	header[8];		// 8 is the maximum size that can be checked
	/* open file and test for it being a png */

	FILE	*fp = fopen(filename, "rb");
	if (!fp)
		abort_("[LoadPNG] File %s could not be opened for reading", filename);
	fread(header, 1, 8, fp);
	if (png_sig_cmp(header, 0, 8))
		abort_("[LoadPNG] File %s is not recognized as a PNG file", filename);


	/* initialize stuff */
	png_ptr = png_create_read_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
	
	if (!png_ptr)
		abort_("[LoadPNG] png_create_read_struct failed");

	info_ptr = png_create_info_struct (png_ptr);
	if (!info_ptr)
		abort_("[LoadPNG] png_create_info_struct failed");

	if (setjmp(png_jmpbuf(png_ptr)))
		abort_("[LoadPNG] Error during init_io");

	png_init_io (png_ptr, fp);
	png_set_sig_bytes(png_ptr, 8);

	png_read_info(png_ptr, info_ptr);

	x	=	png_get_image_width(png_ptr, info_ptr);
	y	=	png_get_image_height(png_ptr, info_ptr);
	*(row) = y;
	*(col) = x;

	color_type = png_get_color_type(png_ptr, info_ptr);
	bit_depth  = png_get_bit_depth(png_ptr, info_ptr);

	number_of_passes = png_set_interlace_handling(png_ptr);
	png_read_update_info(png_ptr, info_ptr);
	
	/* read file */
	if (setjmp(png_jmpbuf(png_ptr)))
		abort_("[LoadPNG] Error during read_image");


	row_pointers = (png_bytep*) malloc(sizeof(png_bytep) * y);
	gray	= (double **) malloc(sizeof(point *) * y);
	for (i=0;i<y;i++){
		row_pointers[i] = (png_byte*)malloc(png_get_rowbytes(png_ptr,info_ptr));
		gray[i] = (double *) malloc(sizeof(point) * x);
	}

	png_read_image(png_ptr, row_pointers);
	
	printf("bit depth: %d \n",bit_depth);
	printf("color type: %d\n",png_get_color_type(png_ptr, info_ptr));

//	abort();
/*
	if (png_get_color_type(png_ptr, info_ptr) == PNG_COLOR_TYPE_RGB)
		abort_("[process_file] input file must be RGB type");

	if (png_get_color_type(png_ptr, info_ptr) == PNG_COLOR_TYPE_RGBA)
		abort_("[process_file] color_type of input file must be PNG_COLOR_TYPE_RGBA (%d) (is %d)",	PNG_COLOR_TYPE_RGBA, png_get_color_type(png_ptr, info_ptr));
*/
	
	for (i=0;i<y;i++){
		png_byte* row = row_pointers[i];
		for (j=0;j<x;j++){
			png_byte* ptr = &(row[j*3]);

			
//			printf("Pixel at position [ %d - %d] has RGB values: %d - %d - %d\n", i, j, ptr[0], ptr[1], ptr[2]);
		
			
			gray[i][j] = 0.2989 * ptr[0] + 0.5870 * ptr[1] + 0.1140 * ptr[2];

//			printf("Pixel at position [ %d - %d] has intensity values: %lf\n", i, j, gray[i][j]);

		}
	}


	fclose(fp);
	return	gray;
}

/****************************************************************************
 * 							Loading JPEG file								*/

point ** LoadJPEG(char * filename, unsigned long * row, unsigned long * col)
/*****************************************************************************/
{
	unsigned long	x, y, i, j;
	unsigned int	data_size;				// length of the file
	int				channels;				// 3=> RGB		4 => RGBA
	unsigned int	type;
	unsigned char *	rowptr[1];				// pointer to an array
	unsigned char *	jdata;					// data for the image

	unsigned int	red, green, blue;
	struct jpeg_decompress_struct	info;	// for our jpeg info
	struct jpeg_error_mgr			err;	// the error handler
	point	** gray;

	FILE * file = fopen (filename, "rb");	// open the file

	info.err = jpeg_std_error (& err);
	jpeg_create_decompress (& info);		// fills info structure

	// if the jpeg file doesn't load
	if(!file){
		fprintf (stderr, "Error reading JPEG file %s!", filename);
		return	0;
	}

	jpeg_stdio_src(&info, file);
	jpeg_read_header(&info, TRUE);			// read jpeg file header

	jpeg_start_decompress(&info);			// decompress the file

	// set width and height
	x = info.output_width;
	y = info.output_height;
	channels = info.num_components;
	
	printf("%d channels image!\n", channels);

	data_size	= x * y * channels;

	//---------------------------------------------
	// read scanlines one at a time & put bytes
	// 		in jdata[] array. Assumes an RGB image
	//---------------------------------------------
	gray  = (point **) malloc (sizeof(point *) * y);
	
	for (i=0;i<y;i++)
		gray[i] = (point *) malloc (sizeof(point) * x);

	jdata = (unsigned char *) malloc (sizeof(unsigned char) * data_size);

	while (info.output_scanline < info.output_height)	// loop
	{
		// Enable jpeg_read_scanlines() to fill our jdata array
		rowptr[0] = (unsigned char *) jdata +
				channels * info.output_width * info.output_scanline;

		i = info.output_scanline;

		jpeg_read_scanlines(&info, rowptr, 1);

		if(channels == 3)
			for (j=0;j<info.output_width;j++){

				gray[i][j].R	= *(jdata + 3 * (x * i + j));
	
				gray[i][j].G	= *(jdata + 3 * (x * i + j) + 1);
	
				gray[i][j].B	= *(jdata + 3 * (x * i + j) + 2);

				gray[i][j].A = 0.2126*gray[i][j].R + 0.7152*gray[i][j].G + 0.0722*gray[i][j].B;
										
			}
		else if(channels == 1)
			for (j=0;j<info.output_width;j++){
				gray[i][j].A = *(jdata + x * i + j);

				gray[i][j].R = gray[i][j].A;

				gray[i][j].G = gray[i][j].A;

				gray[i][j].B = gray[i][j].A;


			}

	}
	//----------------------------------------------
	
	jpeg_finish_decompress(&info);		// finish decompressing

	*row = y;
	*col = x;

	return	gray;					 
}

/*****************************************************************************
 * 						Write Image to JPEG									 *
 *****************************************************************************/
void	WriteJPEG(char * filename, point ** p, unsigned long row, unsigned long col){
	struct	jpeg_compress_struct	cinfo;
	struct	jpeg_error_mgr			jerr;
	unsigned char *row_pointer[1];
	unsigned char *image_buffer;
	unsigned long	i, j, row_stride;
	FILE	*fp;

	image_buffer = (unsigned char *) malloc (sizeof(unsigned char)*3*row*col);

	for(i=0;i<row;i++)
		for(j=0;j<col;j++){
			if(p[i][j].particle){
				image_buffer[3*(i*col + j)] = 0;
				image_buffer[3*(i*col + j) + 1] = 0;
				image_buffer[3*(i*col + j) + 2] = 255;
			}
			else{
				image_buffer[3*(i*col + j)] = p[i][j].R;
				image_buffer[3*(i*col + j) + 1] = p[i][j].G;
				image_buffer[3*(i*col + j) + 2] = p[i][j].B;

			}
		}


	fp = fopen(filename, "w+");

	jpeg_create_compress(&cinfo);

	jpeg_stdio_dest(&cinfo, fp);

	cinfo.image_width = col;
	cinfo.image_height= row;
	cinfo.input_components = 3;
	cinfo.in_color_space = JCS_RGB;
	cinfo.err = jpeg_std_error(&jerr);

	jpeg_set_defaults(&cinfo);

	jpeg_set_quality(&cinfo, 100, true);

	jpeg_start_compress(&cinfo, true);
	
	row_stride = col * 3;

	while(cinfo.next_scanline < cinfo.image_height){
		row_pointer[0] = & image_buffer[cinfo.next_scanline * row_stride];
		(void) jpeg_write_scanlines(&cinfo, row_pointer, 1);
	}

	jpeg_finish_compress(&cinfo);

	fclose(fp);

	jpeg_destroy_compress(&cinfo);
}

/*****************************************************************************
 * 							Draw a Circle									 *
 *****************************************************************************/

int	Draw_circle(point **p, unsigned long row, unsigned long col, double s1, double s2, int Radius, unsigned char R, unsigned char G, unsigned char B, int line_width){
	int		x, y;					// Coordinate of center
	int		i, j;

	x = int(row * s1);
	y = int(col * s2);

	if(s1<0 || s2<0)
		return 0;
	if(s1>=1 || s2>=1)
		return 0;

	for(i=-Radius-line_width;i<Radius+line_width;i++)
		for(j=-Radius-line_width;j<Radius+line_width;j++){
			if(x+i < 0 || y+j < 0)	continue;
			if(x+i >= row || y+j >= col) continue;

			if(abs(sqrt(i*i+j*j) - Radius) < line_width / 2){
				p[x+i][y+j].R	=	R;
				p[x+i][y+j].G	=	G;
				p[x+i][y+j].B	=	B;	
			}
		}

	return  0;
}

/*****************************************************************************
 * 						Draw a Line											 *
 *****************************************************************************/

int	Draw_line(point ** p, unsigned long row, unsigned long col, double s1_1, double s1_2, double s2_1, double s2_2, unsigned char R, unsigned char G, unsigned char B, int line_width){

	int		i, j;
	int		x1, y1, x2, y2;
	int		xmin = row, xmax = 0, ymin = col, ymax = 0;
	double	err;
	
	if(s1_1 < 0 || s1_2 < 0)
		return 0;
	if(s2_1 < 0 || s2_2 < 0)
		return 0;
	if(s1_1 >= 1 || s1_2 >= 1)
		return 0;
	if(s2_1 >= 1 || s2_2 >= 1)
		return 0;


	x1 = int(s1_1 * row);
	y1 = int(s1_2 * col);
	x2 = int(s2_1 * row);
	y2 = int(s2_2 * col);

	if(x1-line_width < xmin)
		xmin = x1-line_width;
	if(x2-line_width < xmin)
		xmin = x2-line_width;
	if(x1+line_width > xmax)
		xmax = x1+line_width;
	if(x2+line_width > xmax)
		xmax = x2+line_width;
	if(y1-line_width < ymin)
		ymin = y1-line_width;
	if(y2-line_width < ymin)
		ymin = y2-line_width;
	if(y1+line_width > ymax)
		ymax = y1+line_width;
	if(y2+line_width > ymax)
		ymax = y2+line_width;

	for(i=xmin;i<=xmax;i++){
		for(j=ymin;j<=ymax;j++){
				err = (x2-x1)*(j-y1) + (y2-y1)*(x1-i);
				err = err / sqrt((x2-x1)*(x2-x1) + (y2-y1)*(y2-y1));	

				if((i-x1)*(x2-x1)+(j-y1)*(y2-y1)<0 || (i-x2)*(x2-x1)+(j-y2)*(y2-y1)>0)
					continue;

				if(fabs(err) < line_width/2){
					p[i][j].R = R;
					p[i][j].G = G;
					p[i][j].B = B;
				}
		}
	}
	
	return 0;
}


/*****************************************************************************
 * 						Image Restoration									 *
 *****************************************************************************/

int Image_restoration(point ** p, unsigned long row, unsigned long col){
	long	x, y, i, j;
	long	xpi, ypj;
	double	Imax = 0, Imin = 255, sum;
	double	B = 0, K0 = 0, tmp = 0;
	double	Kij;
	double	**m0;

	/********* 	Using equation 4-6 *************/
	for (i=-w;i<=w;i++){
		B = B + exp(-i*i/(4*lamda*lamda));
		tmp = tmp + exp(-i*i/(2*lamda*lamda));
	}

	B = B * B;
	tmp = tmp * tmp;
	K0 = tmp / B - B / (2*w+1) / (2*w+1);

	for (x=0;x<row;x++){
		for (y=0;y<col;y++){
			sum	= 0;

			for (i=-w;i<=w;i++){
				for (j=-w;j<=w;j++){
					xpi = x - i;
					ypj = y - j;

					if (xpi<0)
						xpi = 0;
					if (xpi>=row)
						xpi = row - 1;

					if (ypj<0)
						ypj = 0;
					if (ypj>=col)
						ypj = col - 1;
					
					Kij = (exp(-(i*i + j*j)/(4*lamda*lamda)) / B - 1 / ((2*w+1)*(2*w+1))) / K0;

				/*	Equation (6) */	
					sum = sum + p[xpi][ypj].A * Kij; 

				}

			}

			if (sum < 0) sum = 0;

			p[x][y].A = sum;
		}
	}

	/*********  Normalize Intensity ************/

	for (x=0;x<row;x++){
		for (y=0;y<col;y++){
			if(p[x][y].A > Imax) Imax = p[x][y].A;
			if(p[x][y].A < Imin) Imin = p[x][y].A;
		}
	}

	for (x=0;x<row;x++)
		for (y=0;y<col;y++)
			p[x][y].A = (p[x][y].A - Imin) / (Imax - Imin);

	return 0;
}
/*****************************************************************************/
	/* check if point A[x][y] is a particle, return 1 if yes, 0 if no 
	 * if A[x][y] is a local minimum around nearst w atoms and A[x][y] > r,
	 * then it is a particle estimate center location 
	 * Also, calculate the 0 order as equation (8) */

bool estimate_point(point ** p, unsigned long row, unsigned long col, long x, long y){
	long	i, j;
	long	xpi, ypj;
	int		ans = 0;
	double	min = 1000;
	double	m0 = 0, m2 = 0;
	
	p[x][y].particle = true;

	for(i=-w2;i<=w2;i++){
		for(j=-w2;j<=w2;j++){
			xpi = x + i;
			ypj = y + j;

			if (xpi<0)
				xpi = 0;
			if (xpi>=row)
				xpi = row - 1;

			if (ypj<0)
				ypj = 0;
			if (ypj>=col)
				ypj = col - 1;

			if (i*i + j*j <= w2*w2){
				m0 = m0 + p[xpi][ypj].A;
				m2 = m2 + p[xpi][ypj].A * (i*i + j*j);
			}

			if (p[x][y].A < p[xpi][ypj].A)
				p[x][y].particle = false ;	
			
			if (p[xpi][ypj].A < min)
				min = p[xpi][ypj].A;
		}
	}
	
	p[x][y].m0 = m0;

	p[x][y].m2 = m2 / m0;

	if (p[x][y].A < r)
		p[x][y].particle = false;

	if (x<w2||y<w2||(row-x)<w2||(col-w2)<w2)
		p[x][y].particle = false;

	if (p[x][y].A - min < r)
		p[x][y].particle = false;

	return p[x][y].particle;
}

/*****************************************************************************
 *						Refine the location of point						 *
 *****************************************************************************/
/* Calculate the refined location according to equation (7)	*/
bool refine_point(point ** p, unsigned long row, unsigned long col, long x, long y)
{
	int		i, j, xpi, ypj;
	double ex = 0, ey = 0;
	
	if (p[x][y].particle == false)
		return false;

	do{
		p[x][y].particle = false;

		if (ex > 0.5)	x++;
		if (ex < -0.5)	x--;
		if (ey > 0.5)	y++;
		if (ey < -0.5)	y--;
		
		ex = 0;
		ey = 0;

		for(i=-w2;i<=w2;i++){
			for(j=-w2;j<=w2;j++){
				if (i*i + j*j > w2*w2)
					continue;

				xpi = x + i;
				ypj = y + j;
	
				if (xpi<0)
					xpi = 0;
				if (xpi>=row)
					xpi = row - 1;

				if (ypj<0)
					ypj = 0;
				if (ypj>=col)
					ypj = col - 1;
	
					ex = ex + i*p[xpi][ypj].A / p[x][y].m0;
					ey = ey + j*p[xpi][ypj].A / p[x][y].m0;
			}
		}
	}while(fabs(ex) > 0.5 || fabs(ey) > 0.5);

	p[x][y].x = x + ex;			// y = row coordinate
	p[x][y].y = y + ey;			// x = col coordinate
	p[x][y].particle = true;

	return	p[x][y].particle;
}

/*****************************************************************************/
/*	Sort point	*/

point * sort_point(point ** p, unsigned long row, unsigned long col, int *N){
	int		i, j, k, n = 0;
	int		I, J;
	double	ave0 = 0, ave2 = 0;		// Calculate the average of m0 m2
	double	sigma0 = 0, sigma2 = 0;	// Calculate the stander deviation of m0 m2
	double	S, dm0, dm2;
	double	distance;
	point	* list;

	list = (point *) malloc(sizeof(point) * MAX);

	for(i=0;i<row;i++){
		for(j=0;j<col;j++){
			if(p[i][j].particle){
				p[i][j].ID = n+1;
				memcpy(&list[n+1], &p[i][j], sizeof(p[i][j]));
				list[n+1].x = p[i][j].x;
				list[n+1].y = p[i][j].y;
				list[n+1].row = row;
				list[n+1].col = col;
				list[n+1].s1 = list[n+1].x / row;
				list[n+1].s2 = list[n+1].y / col;
				list[n+1].neigh_num = 0;
				list[n+1].z_ = 14;
				list[n+1].t_l = false;
				
				for(I=1;I<3;I++)
					for(J=1;J<3;J++)
						list[n+1].M[I][J] = 0;

				/* list[i] starts from 1, list[0] = dummy particle */
				n++;

				ave0 += p[i][j].m0;
				sigma0 += p[i][j].m0 * p[i][j].m0;

				ave2 += p[i][j].m2;
				sigma2 += p[i][j].m2 * p[i][j].m2;
			}
		}
	}

	if(n > MAX)
		abort_("Particle number %d exceed MAXIMUM allowed number!\n", n);

/* Using equation
 * sigma = sqrt (E(X^2) - E(X))

	ave0 = ave0 / n;
	ave2 = ave2 / n;
	sigma0 = sigma0 / n;
	sigma2 = sigma2 / n;
	
	sigma0 = sqrt(sigma0 - ave0 * ave0);
	sigma2 = sqrt(sigma2 - ave2 * ave2);
*/

/*	Build neighbor list */
	for(i=1;i<=n;i++){
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
			}
		}
	}

	for(i=1;i<=n;i++){
//		list[i].neigh_num = 28;

		for(k=0;k<list[i].neigh_num;k++){
			j = list[i].neigh_list[k];

			list[i].M[1][1] += (list[i].x - list[j].x) * (list[i].x - list[j].x);

			list[i].M[1][2] += (list[i].x - list[j].x) * (list[i].y - list[j].y);

			list[i].M[2][1] = list[i].M[1][2];

			list[i].M[2][2] += (list[i].y - list[j].y) * (list[i].y - list[j].y);
		}

	}
	*N = n;

	return list;
}
/*****************************************************************************
 * 					Return the Value of Cost Function						 *
 *****************************************************************************/
double	cost(point * list1, int i, point * list2, int j){
/* G must satisify topology constrain */
	double	phi = 0, inf = 1e200;
	double	weight = 1, s1, s2;

	if(i == 0){
		s1 = list2[j].s1;
		s2 = list2[j].s2;
//		weight += exp((-(s1-0.5)*(s1-0.5) - (s2-0.5)*(s2-0.5)) * mul);
		return	weight*L*L;
	}

	if(j == 0){
		s1 = list1[i].s1;
		s2 = list1[i].s2;
//		weight += exp((-(s1-0.5)*(s1-0.5) - (s2-0.5)*(s2-0.5)) * mul);
		return	weight*L*L;
	}

	phi+= (list1[i].x-list2[j].x) * (list1[i].x-list2[j].x);

	phi+= (list1[i].y-list2[j].y) * (list1[i].y-list2[j].y);

//	phi+= (list1[i].m0-list2[j].m0)*(list1[i].m0-list2[j].m0);

//	phi+= (list1[i].m2-list2[j].m2)*(list1[i].m2-list2[j].m2);	
/*
	if(phi > L*L)
		return inf;
*/
	return phi;
}



/*****************************************************************************
 * 		Linking Initialization, find nearst neighbor for each atom			 *
 *****************************************************************************/

int	** link_initialization(point * list1, int n1, point * list2, int n2){
/*	Return a (n1+1)*(n2+1) matrix, G, which is defined as equation (12)	*/
	int	**G;			// topology matrix
	int	i, j, best;
	double	min, phi_ij;
	double	inf = 1e150;
	
	G = (int **) malloc (sizeof(int *)*(n1+1));

	for(i=0;i<n1+1;i++)
		G[i] = (int *) malloc (sizeof(int)*(n2+1));

	list1[0].pair_next = 0;
	list2[0].pair_prev = 0;

	G[0][0] = 0;

	for(i=1;i<n1;i++){
		list1[i].pair_next = 0;
		G[i][0] = 1;
	}

	for(j=1;j<=n2;j++){
		list2[j].pair_prev = 0;
		G[0][j] = 1;
	}

	for(i=1;i<=n1;i++)
		for(j=1;j<=n2;j++)
			G[i][j] = 0;

	for(i=1;i<=n1;i++){
		min = L*L;

		best = 0;

		for(j=1;j<=n2;j++){

			if(list2[j].pair_prev)
				continue;

			phi_ij = cost(list1, i, list2, j);
			
			if(phi_ij < min){
				best = j;
				min = phi_ij;

			}else if(phi_ij > L*L){
				G[i][j] = -1;
			/* G[i][j] = -1 means gij will never be considered for this pair*/
			}
				
		}

		if(G[i][best] == -1){
			list1[i].pair_next = 0;
			G[i][0] = 1;
		}

		if(list2[best].pair_prev == 0){
			list1[i].pair_next = best;
			list2[best].pair_prev = i;
	/* Cancel the linking to dummy particle and create linking between i and j*/

			G[i][0] = 0;
			G[0][best] = 0;
			G[i][best] = 1;
		}
	}

	return G;
}

/*****************************************************************************
 * 						Total Cost Function Minimization					 *
 *****************************************************************************/

double	link_optimization(point * list1, int n1, point * list2, int n2, int ** G){
	int		i, j, k, l;
	int		m, t;
	int		better, tmp;
	double	z, min;				// z is the "reduced cost" using eq. (15);
	int		flag = 0;			// if all zij >=0, flag = 0, else flag = 1;
	double	phi_total = 0;
	int		n = 0;

	do{
		flag = 0;

		for(i=0;i<=n1;i++){
			min = 1e200;
			if(i)
				l = list1[i].pair_next;
			else
				l = 0;
			better = l;

			for(m=0;m<=list1[i].neigh_num;m++){
				k = list1[i].neigh_list[m];

				j = list1[k].pair_next;

				if(i==0&&j==0) continue;
				if(G[i][j] < 0) continue;
			
				if(i==0){
					/* an appearing particle ?*/	
					z = cost(list1, 0, list2, j) - \
						cost(list1, k, list2, j) + \
						cost(list1, k, list2, 0);
				
				}else if(j==0){
					/* a disappearing particle */
					z = cost(list1, i, list2, 0) - \
						cost(list1, i, list2, l) + \
						cost(list1, 0, list2, l);
				
				}else{	
					z = cost(list1, i, list2, j) - \
						cost(list1, i, list2, l) - \
						cost(list1, k, list2, j) + \
						cost(list1, k, list2, l);
				}
				
				if(z < min){
					better = j;
					min = z;
				}

			}

			if(min < -0.00001){
				j = better;
				l = list1[i].pair_next;
				k = list2[j].pair_prev;

				list1[i].pair_next = j;
				list2[j].pair_prev = i;
				list1[k].pair_next = l;
				list2[l].pair_prev = k;

				G[i][j] = 1;
				G[k][l] = 1;
				G[i][l] = 0;
				G[k][j] = 0;

				flag ++;
				
				if(n%10==0)
					printf("#%05d iteration!\t\tdphi = %g\n", n, min);

				n++;

				break;
			}
		}

	}while(flag&&n<999);
	
	G[0][0] = 0;

	for(i=1;i<=n1;i++)
		for(j=1;j<=n2;j++){
			if(G[i][j] < 0)
				G[i][j] = 0;
			else if(G[i][j] == 1){
				phi_total += cost(list1, i, list2, j);	
			//	printf("%d\t%d\t%lf\n",i,j,cost(list1, i, list2, j));
			}
		}

	printf("\t\tLinking Finished!\n");

	for(i=1;i<=n1;i++){
		k = list1[i].pair_next;

		for(m=1;m<3;m++)
			for(t=1;t<3;t++){
				list1[i].M_prime[m][t] = 0;
				list1[i].M[m][t] = 0;
			}

		for(m=0;m<list1[i].neigh_num;m++){
			j = list1[i].neigh_list[m];
			l = list1[j].pair_next;

			if(l==0) continue;
			
			list1[i].M[1][1] += (list1[i].x - list1[j].x) * (list1[i].x - list1[j].x);
			list1[i].M[1][2] += (list1[i].x - list1[j].x) * (list1[i].y - list1[j].y);
			list1[i].M[2][1] = list1[i].M[1][2];
			list1[i].M[2][2] += (list1[i].y - list1[j].y) * (list1[i].y - list1[j].y);



			list1[i].M_prime[1][1] += (list1[i].x - list1[j].x) * (list2[k].x - list2[l].x);
			list1[i].M_prime[1][2] += (list1[i].x - list1[j].x) * (list2[k].y - list2[l].y);
			list1[i].M_prime[2][1] += (list1[i].y - list1[j].y) * (list2[k].x - list2[l].x);
			list1[i].M_prime[2][2] += (list1[i].y - list1[j].y) * (list2[k].y - list2[l].y);

		}
	}

	return phi_total;
}

/******************************************************************************
 * 							Translation Offset								  *
 ******************************************************************************/

int	offset_translation(point * list1, int n1, unsigned long row2, unsigned long col2, point * list2, int n2, int **G){
	point	* list1_tmp, * list1_best, * list2_tmp, *list2_best;
	double	phi_best, phi1;
	int		**G_1;
	int		i, j, k, l;
	int		ii, jj;
	double	x1, x2, y1, y2;

	x1 = 0;
	x2 = 0;
	y1 = 0;
	y2 = 0;

	G_1 = (int **) malloc (sizeof(int *)*(n1+1));

	for(i=0;i<n1+1;i++)
		G_1[i] = (int *) malloc (sizeof(int)*(n2+1));

	for(i=1;i<n1+1;i++){
		x1 += list1[i].x;
		y1 += list1[i].y;
	}
	
	x1 = x1 / n1;
	y1 = y1 / n1;

	for(i=1;i<n2+1;i++){
		x2 += list2[i].x;
		y2 += list2[i].y;
	}

	x2 = x2 / n2;
	y2 = y2 / n2;

	printf("Center of mass: %lf %lf, %lf %lf\n", x1, y1, x2, y2);
/*
	for(i=1;i<n2+1;i++){
		list2[i].x = list2[i].x + (x1 - x2);
		list2[i].s1 = list2[i].x / row1;
		list2[i].y = list2[i].y + (y1 - y2);
		list2[i].s2 = list2[i].y / col1;
	}


	link_optimization(list1, n1, list2, n2, G);

	return 0;
*/
	list1_tmp = (point *) malloc (sizeof(point) * (n1+1));
	list2_tmp = (point *) malloc (sizeof(point) * (n2+1));
	list1_best = (point *) malloc (sizeof(point) *(n1+1));
	list2_best = (point *) malloc (sizeof(point) * (n2+1));
	
	listcpy(list1_tmp, list1, n1);
	listcpy(list2_tmp, list2, n2);
	listcpy(list1_best, list1, n1);
	listcpy(list2_best, list2, n2);

	phi_best = link_optimization(list1, n1, list2, n2, G);

	for(i=-offset;i<=offset;i++){
		for(j=-offset;j<=offset;j++){
			if(i==0 && j==0) continue;
			listcpy(list2_tmp, list2, n2);

			for(k=1;k<=n1;k++){
				list2_tmp[k].x = list2[k].x + i;
				list2_tmp[k].y = list2[k].y + j;
				list2_tmp[k].s1 = list2[k].s1 + double(i) / row2;
				list2_tmp[k].s2 = list2[k].s2 + double(j) / col2; 
			}

			G_1	 = link_initialization(list1_tmp, n1, list2_tmp, n2);
			phi1 = link_optimization(list1_tmp, n1, list2_tmp, n2, G_1);
			if(phi_best > phi1){

				listcpy(list2_best, list2_tmp, n1);
				listcpy(list1_best, list1_tmp, n2);
	/*			
				for(ii=0;ii<=n1;ii++)
					for(jj=0;jj<=n2;jj++)
						G[ii][jj] = G_1[ii][jj];
*/
				phi_best = phi1;
			}

			printf("%g\t%g\n", phi_best, phi1);

		}
	}

	listcpy(list1, list1_best, n1);
	listcpy(list2, list2_best, n2);

	return 0;
}

/*****************************************************************************
 * 			finding lost atoms using three images							 *
 *****************************************************************************/

int	find_lost_atoms(point ** p[MAX_image], point ** plist, int *n1, int *n2, int *n3, int NUM_list2){
	/* three scenarios
	 * 1.lost atoms in list1
	 * 2.lost atoms in list2
	 * 3.lost atoms in list3
	 * other scenarios are ignored	*/
	int		i, j, k;
	int		ii, jj, kk;
	int		I, J, K;
	int		x, y;
	bool	flag;
	double	min, dis;
	int		num_atom = 0;		// number of lost atoms
	point	*list1, *list2, *list3;

	list1 = plist[NUM_list2 - 1];
	list2 = plist[NUM_list2];
	list3 = plist[NUM_list2 + 1];

	//	scenario 1 & 2
	for(k=1;k<=*n3;k++){
		flag = false;	
		for(ii=0;ii<Narea;ii++){
			if(list3[k].s2 < xlo[ii])	flag = true;
			if(list3[k].s2 > xhi[ii])	flag = true;
			if(list3[k].s1 < ylo[ii])	flag = true;
			if(list3[k].s1 > yhi[ii])	flag = true;
		}
		
		if(flag) continue;

		j = list3[k].pair_prev;

		if(j == 0){
			// scenario 2
			flag = true;
			ii	= 0;	
			i = 0;
			min = L*L;

			for(K=0;K<list3[k].neigh_num;K++){
				kk = list3[k].neigh_list[K];
				jj = list3[kk].pair_prev;
				if(jj == 0) continue;
				ii = list2[jj].pair_prev;
				if(ii == 0) continue;
							for(J=0;J<list1[ii].neigh_num;J++){
					I = list1[ii].neigh_list[J];
					if(list1[I].pair_next == 0){
						dis = (list1[I].x - list3[k].x) * (list1[I].x - list3[k].x) + (list1[I].y - list3[k].y) * (list1[I].y - list3[k].y);
						if(dis < min){
							min = dis;
							i	= I;
						}
					}
				}

			}

			if(i){
				num_atom++;
				*n2 = *n2 + 1;
				/*using the nearest four atoms*/
				list2[*n2].x = (list1[i].x + list3[k].x) / 2;
				list2[*n2].y = (list1[i].y + list3[k].y) / 2;
				list2[*n2].neigh_num = 0;
				build_neigh_for_new_atom(list2, *n2);
				list2[*n2].x = 0;
				list2[*n2].y = 0;

				for(I=0;I<5;I++){

					list2[*n2].x += list2[list2[*n2].neigh_list[I]].x/5;

					list2[*n2].y += list2[list2[*n2].neigh_list[I]].y/5;
				}

				list2[*n2].s1 = list2[*n2].x / list2[1].row;
				list2[*n2].s2 = list2[*n2].y / list2[1].col;

				if(list2[*n2].distance_list[0] < w2){
					*n2 = *n2 - 1;
					num_atom --;
				}

			}
		}
		else if (list2[j].pair_prev == 0){
			// scenario 1
				num_atom++;
				*n1 = *n1 + 1;
				/*using the nearest four atoms*/
				list1[*n1].x = (list2[j].x + list3[k].x) / 2;
				list1[*n1].y = (list2[j].y + list3[k].y) / 2;
				list1[*n1].neigh_num = 0;
				build_neigh_for_new_atom(list1, *n1);
				list1[*n1].x = 0;
				list1[*n1].y = 0;

				for(I=0;I<5;I++){

					list1[*n1].x += list1[list1[*n1].neigh_list[I]].x/5;

					list1[*n1].y += list1[list1[*n1].neigh_list[I]].y/5;
				}

				list1[*n1].s1 = list1[*n1].x / list1[1].row;
				list1[*n1].s2 = list1[*n1].y / list1[1].col;


				if(list1[*n1].distance_list[0] < w2){
					*n1 = *n1 - 1;
					num_atom --;
				}

		}
	}	

		
	for(j=1;j<=*n2;j++){
		flag = false;	
		for(ii=0;ii<Narea;ii++){
			if(list2[j].s2 < xlo[ii])	flag = true;
			if(list2[j].s2 > xhi[ii])	flag = true;
			if(list2[j].s1 < ylo[ii])	flag = true;
			if(list2[j].s1 > yhi[ii])	flag = true;
		}
		
		if(flag) continue;
		i = list2[j].pair_prev;
		if(i){
			if(list2[j].pair_next == 0){
				// scenario 3

				num_atom++;
				*n3 = *n3 + 1;
				/*using the nearest four atoms*/
				list3[*n3].x = (list1[i].x + list2[j].x) / 2;
				list3[*n3].y = (list1[i].y + list2[j].y) / 2;
				list3[*n3].neigh_num = 0;
				build_neigh_for_new_atom(list3, *n3);
				list3[*n3].x = 0;
				list3[*n3].y = 0;

				for(I=0;I<5;I++){

					list3[*n3].x += list3[list3[*n3].neigh_list[I]].x/5;

					list3[*n3].y += list3[list3[*n3].neigh_list[I]].y/5;
				}

				list1[*n3].s1 = list1[*n3].x / list3[1].row;
				list1[*n3].s2 = list1[*n3].y / list3[1].col;

				if(list3[*n3].distance_list[0] < w2){
					*n3 = *n3 - 1;
					num_atom --;
				}

			}
			
		}
	}


	return num_atom;
}


/*****************************************************************************
 *						final_result	print result into image 			 *
 *****************************************************************************/

point ** final_result(point * list1, int n1, point * list2, int n2, int **G_1, point * list3, int n3, int **G_2, int row, int col){
	int		i, j, k;
	point	** p;
	double	exx, eyy, exy;
	double	vmise;
	FILE	*fstrain;

	fstrain = fopen("strain","w");

	p = (point **) malloc(sizeof(point *) * row);

	for(i=0;i<row;i++){
		p[i] = (point *) malloc(sizeof(point) * col);
		for(j=0;j<col;j++){
			p[i][j].particle = 0;
			p[i][j].R  = 255;
			p[i][j].G  = 255;
			p[i][j].B  = 255;
		}
	}
	
	for(i=1;i<=n1;i++){
		Draw_circle(p, row, col, list1[i].s1, list1[i].s2, r1, R1, G1, B1, line_width);

		j = list1[i].pair_next;

//		printf("%d\t%d\t%lf\t%lf\n",i,j,list1[i].s1, list2[j].s1);

		if(j!=0){
			Draw_line(p, row, col, list1[i].s1, list1[i].s2, list2[j].s1, list2[j].s2, R1, G1, B1, line_width);

			strain(list1[i].M, list1[i].M_prime, &exx, &eyy, &exy);

			vmise = sqrt(exy*exy + (eyy*eyy+exx*exx+(exx-eyy)*(exx-eyy))/6);

			if(fabs(100*vmise) < 1)
				fprintf(fstrain, "%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%d\n", list1[i].x, list1[i].y, 100*exx, 100*eyy, 100*exy, 100*vmise, list1[i].neigh_num);
		}
	}

	for(i=1;i<=n2;i++){	
		Draw_circle(p, row, col, list2[i].s1, list2[i].s2, r2, R2, G2, B2, line_width);	

		j = list2[i].pair_next;

		if(j!=0){
			Draw_line(p, row, col, list2[i].s1, list2[i].s2, list3[j].s1, list3[j].s2, 5, 5, 5, line_width);
		}

	}

	for(i=1;i<=n3;i++){
		Draw_circle(p, row, col, list3[i].s1, list3[i].s2, r3, R3, G3, B3, line_width);
	}

	return	p;
}


int WriteCFG(point **p, int N, double X, double Y, double Z, double mass, char *ele){
	    int     i, j, k, count = 0;
	    FILE    *fp1, *fp2, *fp3;
			
	    fp1 = fopen("01.cfg","w");
		
		fp2 = fopen("02.cfg","w");

//		fp3 = fopen("03.cfg","w");
			
		fprintf(fp1,"Number of particles = %d\n\n",N);
		fprintf(fp1,"A = 1.0 Angstrom (basic length-scale)\n\n");
		fprintf(fp1,"H0(1,1) = %.3lf A\nH0(1,2) = 0 A\nH0(1,3) = 0 A\n\n", X);
		fprintf(fp1,"H0(2,1) = 0 A\nH0(2,2) = %.3lf A\nH0(2,3) = 0 A\n\n", Y);
		fprintf(fp1,"H0(3,1) = 0 A\nH0(2,3) = 0 A\nH0(3,3) = %.3lf A\n\n", Z);
		fprintf(fp1,"Transform(1,1) = 1\nTransform(1,2) = 0\nTransform(1,3) = 0\nTransform(2,1) = 0\nTransform(2,2) = 1\nTransform(2,3) = 0\nTransform(3,1) = 0\nTransform(3,2) = 0\nTransform(3,3) = 1\n\n");
		fprintf(fp1,"eta(1,1) = 0\neta(1,2) = 0\neta(1,3) = 0\neta(2,2) = 0\neta(2,3) = 0\neta(3,3) = 0\n\n\n\n");
		fprintf(fp1,"R = 1.0 [ns^-1]\n\n");


		fprintf(fp2,"Number of particles = %d\n\n",N);
		fprintf(fp2,"A = 1.0 Angstrom (basic length-scale)\n\n");
		fprintf(fp2,"H0(1,1) = %.3lf A\nH0(1,2) = 0 A\nH0(1,3) = 0 A\n\n", X);
		fprintf(fp2,"H0(2,1) = 0 A\nH0(2,2) = %.3lf A\nH0(2,3) = 0 A\n\n", Y);
		fprintf(fp2,"H0(3,1) = 0 A\nH0(2,3) = 0 A\nH0(3,3) = %.3lf A\n\n", Z);
		fprintf(fp2,"Transform(1,1) = 1\nTransform(1,2) = 0\nTransform(1,3) = 0\nTransform(2,1) = 0\nTransform(2,2) = 1\nTransform(2,3) = 0\nTransform(3,1) = 0\nTransform(3,2) = 0\nTransform(3,3) = 1\n\n");
		fprintf(fp2,"eta(1,1) = 0\neta(1,2) = 0\neta(1,3) = 0\neta(2,2) = 0\neta(2,3) = 0\neta(3,3) = 0\n\n\n\n");
		fprintf(fp2,"R = 1.0 [ns^-1]\n\n");

/*
		fprintf(fp3,"Number of particles = %d\n\n",2020);
		fprintf(fp3,"A = 1.0 Angstrom (basic length-scale)\n\n");
		fprintf(fp3,"H0(1,1) = %.3lf A\nH0(1,2) = 0 A\nH0(1,3) = 0 A\n\n", X);
		fprintf(fp3,"H0(2,1) = 0 A\nH0(2,2) = %.3lf A\nH0(2,3) = 0 A\n\n", Y);
		fprintf(fp3,"H0(3,1) = 0 A\nH0(2,3) = 0 A\nH0(3,3) = %.3lf A\n\n", Z);
		fprintf(fp3,"Transform(1,1) = 1\nTransform(1,2) = 0\nTransform(1,3) = 0\nTransform(2,1) = 0\nTransform(2,2) = 1\nTransform(2,3) = 0\nTransform(3,1) = 0\nTransform(3,2) = 0\nTransform(3,3) = 1\n\n");
		fprintf(fp3,"eta(1,1) = 0\neta(1,2) = 0\neta(1,3) = 0\neta(2,2) = 0\neta(2,3) = 0\neta(3,3) = 0\n\n\n\n");
		fprintf(fp3,"R = 1.0 [ns^-1]\n\n");
*/
		for(i=1;i<=N;i++){
//			if(p[1][i].pair_next)
//				if(p[2][p[1][i].pair_next].pair_next){
//					j = p[1][i].pair_next;
//					k = p[2][j].pair_next;
					
					
						
			fprintf(fp1,"%.2lf %s %.7lf %.7lf %.7lf 0 0 0\n", mass, ele, p[1][i].s2, p[1][i].s1, p[1][i].z_/Z+0.5);
			
			count++;
					

//					fprintf(fp2,"%.2lf %s %.7lf %.7lf %.7lf 0 0 0\n", mass, ele, p[2][j].s2, p[2][j].s1, 0.25);
//					fprintf(fp3,"%.2lf %s %.7lf %.7lf %.7lf 0 0 0\n", mass, ele, p[3][k].s2, p[3][k].s1, 0.25);
//					fprintf(fp1,"%.2lf %s %.7lf %.7lf %.7lf 0 0 0\n", mass, ele, p[1][i].s2, p[1][i].s1, 0.75);
//					fprintf(fp2,"%.2lf %s %.7lf %.7lf %.7lf 0 0 0\n", mass, ele, p[2][j].s2, p[2][j].s1, 0.75);
//					fprintf(fp3,"%.2lf %s %.7lf %.7lf %.7lf 0 0 0\n", mass, ele, p[3][k].s2, p[3][k].s1, 0.75);
	
			//	}


		}
		printf("%d\n",count);
		return 0;
}

int	mark_z(point *p, int n){
	int		i, j;
	int		m, k, max;
	double	xy_, x_, y_, x2_, y2_;
	double	A = 0, B = 0, r = 0, delta;
	double	sum0 = 0, sum1 = 0;

	if(p[n].neigh_num > 8) max = 8;
	else max = p[n].neigh_num;

	for(i=0;i<max;i++){	
		m = p[n].neigh_list[i];
		r = (p[n].x-p[m].x)*(p[n].x-p[m].x);
		r = r + (p[n].y-p[m].y)*(p[n].y-p[m].y);
		sum0 = sum0 + sqrt(r+(p[n].z_-p[m].z_)*(p[n].z_-p[m].z_));
		sum1 = sum1 + sqrt(r+(p[n].z_+p[m].z_)*(p[n].z_+p[m].z_));
	}

	if(sum0>sum1)
		return 1;
	else
		return -1;
}

/*****************************************************************************
 * 						Main Program										 *
 *****************************************************************************/

int	main(){
	point	** p[MAX_image];
	point	** output;
	point	* list[MAX_image];
	char	ele[5] = "Pt";

	/* 0 order intensity momentum, 2 order intensity momentum */
	int		n[MAX_image];				// number of particle in IMAGE;
	unsigned long	row[MAX_image], col[MAX_image], i, j;
	int		k;
	int		**G[MAX_image];				// topology matrix for list1 and list2
	int		lost, iii = 0;
	bool	flag = false;
	int		tmp, n1 = 1;
	
	read_input();
	
	for(k=1;k<Nimage;k++){
		printf("Reading file %s\n", file[k]);
		p[k] = LoadJPEG(file[k], &row[k], &col[k]);
		Image_restoration(p[k], row[k], col[k]);
		for(i=0;i<row[k];i++){
			for(j=0;j<col[k];j++){
				estimate_point(p[k], row[k], col[k], i, j);
				refine_point(p[k], row[k], col[k], i, j);
			}
		}

		list[k] = sort_point(p[k], row[k], col[k], &n[k]);
	}

//	while(flag==false){
		flag = true;
		for(i=1;i<=n[1];i++){	
			tmp = mark_z(list[1],i);
			list[1][i].z_ = tmp *list[1][i].z_;
			if(tmp<0) flag = false;
		}
//	}


	WriteCFG(list, n[1], col[1], row[1], 3.92*5.414, 195, ele);

	return 0;
}
