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

#include "cfg2RDF.h"

int		read_input(){
	char		line[255];
	int			count = 0;
	int			i = 0, j = 0;

	while(fgets(line, 255, stdin)){
		if(line[0] == '#') continue;

		printf("%s", line);

		if(count == 0)
			sscanf(line, "%s%s%s", cfgfile, mass, element);
		if(count == 1)
			sscanf(line, "%lf%lf%lf%lf", &xlo, &xhi, &ylo, &yhi);
		if(count == 2)
			sscanf(line, "%s", output);
		if(count == 3)
			sscanf(line, "%d%lf%d", &N0, &dr, &rmax_dr);
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
    for(i=0;i<Natom;i++){
		flag=fscanf(fp,"%lf%lf%lf",&s1,&s2,&s3);

		if(flag<3){
			Natom=i;
			break;
		}
        (Cu+i)->x = s1; 
        (Cu+i)->y = s2;
        (Cu+i)->z = s3;
		if(s1>xlo&&s1<xhi){
			if(s2>xlo&&s2<xhi){
				if(s3>xlo&&s3<xhi){
					RDF_atom[N_RDF] = i;
					N_RDF = N_RDF + 1;
				}
			}
		}
   	}
    fclose(fp);

	return Cu; 
}



/*					MAIN		PROGRAM					*/
int	main(){
	int		i, j;
	int		m,n;
	int		ilo, ihi, jlo, jhi;
	double	x1,y1,z1,x2,y2,z2,r, PI, rho;
	double	*rdf;
	atom	*Cu;
	FILE	*fp;

	PI = atan(1.0)*4.0;

	srand(time(NULL));
	read_input();
	Cu = Read_CFG();
	rho = Natom/H[1][1]/H[2][2]/H[3][3];
	rdf = (double*)malloc(sizeof(double)*rmax_dr);
	for(i=0;i<rmax_dr;i++)
		rdf[i] = 0.0;


	for(m=0;m<N0;m++){
		j = RDF_atom[rand()%N_RDF];
		x1 = (Cu+j)->x*H[1][1];
		y1 = (Cu+j)->y*H[2][2];
		z1 = (Cu+j)->z*H[3][3];

		for(i=0;i<Natom;i++){
			if(i==j) continue;
			x2 = (Cu+i)->x*H[1][1];
			y2 = (Cu+i)->y*H[2][2];
			z2 = (Cu+i)->z*H[3][3];
			r = (x2-x1)*(x2-x1) + (y2-y1)*(y2-y1) + (z2-z1)*(z2-z1);
			r = sqrt(r);
			n = int(r/dr);
			if(n<rmax_dr){
				rdf[n]=rdf[n]+1.0/(4*PI*(n+1)*(n+1)*dr*dr*rho*dr)/N0;
			}
		}
	}

	fp = fopen(output,"w+");
	for(i=0;i<rmax_dr;i++){
		fprintf(fp,"%lf\t%lf\n",i*dr,rdf[i]);
	}

	fclose(fp);

	return 0;
}

