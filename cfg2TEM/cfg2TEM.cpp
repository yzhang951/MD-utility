/* Convert CFG file to TEM style projection image
 * Author		: Yin Zhang
 * Version		: 1.0
 * Data			: Dec 15th 2016 */

#include <stdio.h>
#include <jpeglib.h>
#include <jerror.h>
#include <string.h>
#include <stdarg.h>
#include <stdlib.h>
#include <ctime>
#include <math.h>

#include "cfg2TEM.h"

int		read_input(){
	char		line[255];
	int			count = 0;
	int			i = 0, j = 0;

	while(fgets(line, 255, stdin)){
		if(line[0] == '#') continue;

		printf("%s", line);

		if(count == 0)
			sscanf(line, "%s%s%s%d", cfgfile, mass, element, &pro_dir);
		if(count == 1)
			sscanf(line, "%lf%lf%lf%lf%lf%lf", &xlo, &xhi, &ylo, &yhi, &zlo, &zhi);
		if(count == 2)
			sscanf(line, "%s%d%d", output, &ROW, &COL);
		if(count == 3)
			sscanf(line, "%d%d%lf", &MAX_RGB, &w, &b);
		count++;
	}
}

/*******************************************************************************
*			Read atom info from CFG file									   ********************************************************************************/

atom *	Read_CFG(){
	char 	command[50];
	FILE 	*fp, *fo;
	int		i, j, k, flag;
	double	s1, s2, s3, tmp;
	atom 	*Cu;

	/* Read DATA! */
	for(i=0;i<50;i++) command[i]='\0';
	strcpy (command,"mul ");
	strcat (command, cfgfile);
	strcat (command, " 1 1 1 ");
	strcat (command, cfgfile);
	printf ("%s\n",command);
	system (command);

    /* cut -f 5 -d" "5" the fifth words  
    READ DATA!*/
    for(i=0;i<50;i++) command[i]='\0';
    strcpy (command,"grep \"Number of particles = \" ");
    strcat (command,cfgfile);
    strcat (command," | cut -f 5 -d\" \" > log1");
    system (command);
    fp = fopen("log1","r");
    fscanf(fp,"%d",&Natom);
    fclose(fp);
    system ("rm log1");
    Cu=(atom*)malloc(sizeof(atom)*Natom);
    
    for(i=0;i<50;i++) command[i]='\0';
    strcpy (command,"grep \"H0(1,1) =\" ");
    strcat (command,cfgfile);
    strcat (command," | cut -f 3 -d\" \" > log1");
    system (command);
    fp = fopen("log1","r");
    fscanf(fp,"%lf",&H[1][1]);
    fclose(fp);
    system ("rm log1");

    
    for(i=0;i<50;i++) command[i]='\0';
    strcpy (command,"grep \"H0(1,2) =\" ");
    strcat (command,cfgfile);
    strcat (command," | cut -f 3 -d\" \" > log1");
    system (command);
    fp = fopen("log1","r");
    fscanf(fp,"%lf",&H[1][2]);
    fclose(fp);
    system ("rm log1");
    
    for(i=0;i<50;i++) command[i]='\0';
    strcpy (command,"grep \"H0(1,3) =\" ");
    strcat (command,cfgfile);
    strcat (command," | cut -f 3 -d\" \" > log1");
    system (command);
    fp = fopen("log1","r");
    fscanf(fp,"%lf",&H[1][3]);
    fclose(fp);
    system ("rm log1");

    for(i=0;i<50;i++) command[i]='\0';
    strcpy (command,"grep \"H0(2,1) =\" ");
    strcat (command,cfgfile);
    strcat (command," | cut -f 3 -d\" \" > log1");
    system (command);
    fp = fopen("log1","r");
    fscanf(fp,"%lf",&H[2][1]);
    fclose(fp);
    system ("rm log1");

    for(i=0;i<50;i++) command[i]='\0';
    strcpy (command,"grep \"H0(2,2) =\" ");
    strcat (command,cfgfile);
    strcat (command," | cut -f 3 -d\" \" > log1");
    system (command);
    fp = fopen("log1","r");
    fscanf(fp,"%lf",&H[2][2]);
    fclose(fp);
    system ("rm log1");
    
    for(i=0;i<50;i++) command[i]='\0';
    strcpy (command,"grep \"H0(2,3) =\" ");
    strcat (command,cfgfile);
    strcat (command," | cut -f 3 -d\" \" > log1");
    system (command);
    fp = fopen("log1","r");
    fscanf(fp,"%lf",&H[2][3]);
    fclose(fp);
    system ("rm log1");
    
    for(i=0;i<50;i++) command[i]='\0';
    strcpy (command,"grep \"H0(3,1) =\" ");
    strcat (command,cfgfile);
    strcat (command," | cut -f 3 -d\" \" > log1");
    system (command);
    fp = fopen("log1","r");
    fscanf(fp,"%lf",&H[3][1]);
    fclose(fp);
    system ("rm log1");

    for(i=0;i<50;i++) command[i]='\0';
    strcpy (command,"grep \"H0(3,2) =\" ");
    strcat (command,cfgfile);
    strcat (command," | cut -f 3 -d\" \" > log1");
    system (command);
    fp = fopen("log1","r");
    fscanf(fp,"%lf",&H[3][2]);
    fclose(fp);
    system ("rm log1");
    
    for(i=0;i<50;i++) command[i]='\0';
    strcpy (command,"grep \"H0(3,3) =\" ");
    strcat (command,cfgfile);
    strcat (command," | cut -f 3 -d\" \" > log1");
    system (command);
    fp = fopen("log1","r");
    fscanf(fp,"%lf",&H[3][3]);
    fclose(fp);
    system ("rm log1");
    
    for(i=0;i<50;i++) command[i]='\0';
    sprintf(command,"grep  \"%s %s\" ",mass,element);
    strcat (command,cfgfile);
    strcat (command,"  | cut -f 3-5 -d \" \">  testdata");
    system (command);
    fp = fopen("testdata","r");
//  fscanf(fp,"%s",&temp);

	if(pro_dir == 1){
		tmp = H[1][1];
		H[1][1] = H[2][2];
		H[2][2] = H[3][3];
		H[3][3] = tmp;		
	}
	else if(pro_dir == 2){
		tmp = H[1][1];
		H[1][1] = H[3][3];
		H[3][3] = H[2][2];
		H[2][2] = tmp;
	}
	else if(pro_dir == 4){
		H[1][1] = H[1][1]/sqrt(2);
		H[3][3] = H[3][3]/sqrt(2);
	}

    for(i=0;i<Natom;i++){
		if(pro_dir == 1)
			flag=fscanf(fp,"%lf%lf%lf",&s3,&s2,&s1);
		else if (pro_dir == 2)
			flag=fscanf(fp,"%lf%lf%lf",&s1,&s3,&s2);
		else if (pro_dir > 2)
			flag=fscanf(fp,"%lf%lf%lf",&s1,&s2,&s3);

		if(flag<3){
			Natom=i;
			break;
		}
        (Cu+i)->x = s1; 
        (Cu+i)->y = s2;
        (Cu+i)->z = s3;
		if(pro_dir == 4){	// xz as 110 direction
			(Cu+i)->x = (s1-s3+0.5);
			(Cu+i)->y = s2;
			(Cu+i)->z = (s1+s3-0.5);
		}

   	}
    fclose(fp);

	return Cu; 
}

/*****************************************************************************
 * 			Convert Grayscale to Hot-to-Cold								 *
 *****************************************************************************/
void	hot2cold(unsigned char *RGB, double x){
	if(x<64){
		RGB[0] = 0;
		RGB[1] = 255*64/x;
		RGB[2] = 255;
	}else if(x<128){
		RGB[0] = 0;
		RGB[1] = 255;
		RGB[2] = 255-255*64/(x-64);
	}else if(x<192){
		RGB[0] = 255*64/(x-128);
		RGB[1] = 255;
		RGB[2] = 0;
	}else{
		RGB[0] = 255;
		RGB[1] = 255-255*64/(x-192);
		RGB[2] = 0;
	}
}

/*****************************************************************************
 * 						Write Image to JPEG									 *
 *****************************************************************************/
void	WriteJPEG(){
	struct	jpeg_compress_struct	cinfo;
	struct	jpeg_error_mgr			jerr;
	unsigned char *row_pointer[1];
	unsigned char *image_buffer;
	unsigned long	i, j, row_stride;
	FILE	*fp;

	image_buffer = (unsigned char *) malloc (sizeof(unsigned char)*3*ROW*COL);

	for(i=0;i<ROW;i++)
		for(j=0;j<COL;j++){
//			image_buffer[3*(i*COL + j)] = 0;//(int) gray[i][j];
			image_buffer[3*(i*COL + j)] = (int) gray[i][j];
			image_buffer[3*(i*COL + j) + 1] = (int) gray[i][j];
//			image_buffer[3*(i*COL + j) + 2] = 0;//(int) gray[i][j];
			image_buffer[3*(i*COL + j) + 2] = (int) gray[i][j];
/* Hot-to-cold color map*/			
//			hot2cold(image_buffer+3*(i*COL + j), gray[i][j]);
		}


	fp = fopen(output, "w+");

	jpeg_create_compress(&cinfo);

	jpeg_stdio_dest(&cinfo, fp);

	cinfo.image_width = COL;
	cinfo.image_height= ROW;
	cinfo.input_components = 3;
	cinfo.in_color_space = JCS_RGB;
	cinfo.err = jpeg_std_error(&jerr);

	jpeg_set_defaults(&cinfo);

	jpeg_set_quality(&cinfo, 100, true);

	jpeg_start_compress(&cinfo, true);
	
	row_stride = COL * 3;

	while(cinfo.next_scanline < cinfo.image_height){
		row_pointer[0] = & image_buffer[cinfo.next_scanline * row_stride];
		(void) jpeg_write_scanlines(&cinfo, row_pointer, 1);
	}

	jpeg_finish_compress(&cinfo);

	fclose(fp);

	jpeg_destroy_compress(&cinfo);
}





/*					MAIN		PROGRAM					*/
int	main(){
	int		i, j;
	int		m;
	int		ilo, ihi, jlo, jhi;
	double	s1, s2, r, max, af;
	int		background;
	atom	*Cu;
	FILE	*fp;

	srand(time(NULL));
	read_input();

	gray = (double **) malloc (sizeof(double *)*ROW);

	for(i=0;i<ROW;i++){
		gray[i] =  (double *) malloc (sizeof(double)*COL);
/*		for(j=0;j<COL;j++){

			background = rand()%5;  
			gray[i][j] = 58 + background;
		}*/
	}

	Cu = Read_CFG();
	max = MAX_RGB;

	for(m=0;m<Natom;m++){
		if(((Cu+m)->x - xlo)*((Cu+m)->x - xhi) >= 0) continue;
		if(((Cu+m)->y - ylo)*((Cu+m)->y - yhi) >= 0) continue;
		if(((Cu+m)->z - zlo)*((Cu+m)->z - zhi) >= 0) continue;

		s2 = ((Cu+m)->x - xlo) / (xhi - xlo) * COL;
		s1 = ((Cu+m)->y - ylo) / (yhi - ylo) * ROW;

		ilo = s1 - 6*w;
		ihi = s1 + 6*w;
		jlo = s2 - 6*w;
		jhi = s2 + 6*w;

		if(ilo < 0) ilo = 0;
		if(jlo < 0) jlo = 0;
		if(ihi >= ROW) ihi = ROW;
		if(jhi >= COL) jhi = COL;

		for(i=ilo;i<ihi;i++){
			for(j=jlo;j<jhi;j++){
				r = (i-s1)*(i-s1) + (j-s2)*(j-s2);
				gray[i][j] += MAX_RGB * exp(-r/w/w/2);
			}
		}
	}

	for(i=0;i<ROW;i++)
		for(j=0;j<COL;j++){
			if(gray[i][j]>max) max = gray[i][j];
		}
	printf("max/MAX_RGB is %lf\n",max/MAX_RGB);

	for(i=0;i<ROW;i++)
		for(j=0;j<COL;j++)
			gray[i][j]=gray[i][j]*MAX_RGB/max;

	fp = fopen("atom_fraction.dat","w");

	for(j=0;j<COL;j++){
		af = 0;
		for(i=0;i<ROW;i++)
			af = af + gray[i][j]/ROW;
		fprintf(fp,"%d %lf\n", j, af);
	}

	fclose(fp);

	WriteJPEG();

	return 0;
}
