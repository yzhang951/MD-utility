/* Struture factor
 * Author		:	Yin Zhang
 * Version		:	0.99
 * Date			:	Mar 2nd 2018		*/

#include "structure_factor.h"

atom * read_cfg(char * cfgfile, int *N, double *qdx, double *qdy,double *f2sum){
	FILE 	*fp;
	int		i, j, Natom;
	double	H[4][4], s1, s2, s3;
	double	f_100[3];
	atom 	*p;
	char	command[100], element[3][3], tmp[3], mass[3][10];
	*f2sum 	= 0;
	strcpy(mass[0],"55.845");
	strcpy(mass[1],"195");
	strcpy(mass[2],"58.71");
//	atom form factor at 100keV
	f_100[0] = 26;
	f_100[1] = 24;
	f_100[2] = 28;

	strcpy(element[0],"Fe");
	strcpy(element[1],"Pt");
	strcpy(element[2],"Ni");

	strcpy(command,"mul ");
	strcat(command,cfgfile);
	strcat(command," 1 1 1 ");
	strcat(command,cfgfile);
	system(command);
	
    for(i=0;i<50;i++) 
		command[i]='\0';
	    
	strcpy (command,"grep \"Number of particles = \" ");
	strcat (command,cfgfile);
	strcat (command," | cut -f 5 -d\" \" > log1");
	system (command);
	fp = fopen("log1","r");
	fscanf(fp,"%d",&Natom);
	printf("There are %d atoms in this cfg file\n", Natom);
	fclose(fp);
	system ("rm log1");
	p=(atom*)malloc(sizeof(atom)*Natom);
	*N = Natom;

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

	Lx = H[1][1];
	Ly = H[2][2];
	Lz = H[3][3];

	*qdx = 2*(4*atan(1.0))/H[1][1];
	*qdy = 2*(4*atan(1.0))/H[2][2];
   	
	for(j=0;j<3;j++){
    	for(i=0;i<100;i++) command[i]='\0';
    	sprintf(command,"grep  \"%s %s\" ",mass[j],element[j]);
    	strcat (command,cfgfile);
    	if(j==0)
			strcat (command,"  | cut -f 2-5 -d \" \">  testdata");
		else
			strcat (command,"  | cut -f 2-5 -d \" \">>  testdata");

    	system (command);
	}

	fp = fopen("testdata","r");

	for(i=0;i<Natom;i++){
		fscanf(fp,"%s%lf%lf%lf",tmp,&s1,&s2,&s3);
		(p+i)->x = s1*H[1][1]+s2*H[2][1]+s3*H[3][1];
		(p+i)->y = s1*H[1][2]+s2*H[2][2]+s3*H[3][2];
		(p+i)->z = s1*H[1][3]+s2*H[2][3]+s3*H[3][3];
		(p+i)->f = 0;
		for(j=0;j<3;j++){
			if(strcmp(element[j],tmp)==0){
				(p+i)->f = f_100[j];
				*f2sum = *f2sum + f_100[j]*f_100[j];
			}
		}	
	}
	return p;
}

double	Sq(atom *p, int Natom, double qx, double qy, double qz, double f2sum){
	int		i, j;
	double	Re = 0, Im = 0, S;
	double	fi, fj, x, y, z, p1, p2;
	int		n = 1;
	double  y_min, y_max;
	y_min = 27;
	y_max = 96;
	
	for(i=0;i<Natom;i++){
		fi 	= (p+i)->f;
		x	= (p+i)->x;
		y	= (p+i)->y;
		z	= (p+i)->z;
//		if(y>y_max||y<y_min) continue;
		Re 	= Re + fi*cos(qx*x+qy*y+qz*z);	
		Im	= Im - fi*sin(qx*x+qy*y+qz*z);
	}
	S = Re*Re + Im*Im;
	S = S/f2sum;
	return S;	
}

int main(int argc, char **argv){
	atom	*p;
	FILE	*fp;
	FILE	*result;
	int		i, j, N, k, n;
	char	cfgfile[100], output[100];
	double 	F, qx=0.0, qy=0.0;
	double	dqx, dqy, dq = 0.0001;
	double	Sdq;
	double	k_guess[4][2];
	double	m1_x, m1_y, m0;
	double	k_result[18][4][2];

	k_guess[0][0] =  0.0188; 
  	k_guess[0][1] =  0.7712;
	k_guess[1][0] = -0.5112;
 	k_guess[1][1] =  0.3738;
	k_guess[2][0] = -0.2388;
   	k_guess[2][1] =  0.7263;	
	k_guess[3][0] = -0.6112;
 	k_guess[3][1] =  0.1739;

	/* k-vector of left grain and right grain, 0,1 left grain, 2,3 right grain*/

	for(k=1;k<19;k++){
		qx = 0;
		sprintf(cfgfile,"../CFG/%02d.cfg",k);
		p = read_cfg(cfgfile,&N,&dqx,&dqy,&F);
		printf("\n");

		for(n=0;n<4;n++){
			sprintf(cfgfile,"k_space_refined/k%02d_%d.dat",k, n+1);	
			fp = fopen(cfgfile,"w");


			m1_x = 0;
			m1_y = 0;
			m0 = 0;

			for(i=-500;i<=500;i++){
				for(j=-500;j<=500;j++){		
					qx = k_guess[n][0]+i*dq;
					qy = k_guess[n][1]+j*dq;
					Sdq = Sq(p,N,qx,qy,0,F);
					fprintf(fp,"%lf ",Sdq);

					m1_x = m1_x + Sdq*qx;
					m1_y = m1_y + Sdq*qy;
					m0 = m0 + Sdq;
					printf("\rIn progress %3.1lf %%",(i+500)/10.0);
				}
				fprintf(fp,"\n");
			}

			k_result[k-1][n][0] = m1_x/m0;
			k_result[k-1][n][1] = m1_y/m0;

			fclose(fp);
			printf("\nFinished k-vecotr #%d\n%lf\t%lf\n",n+1,m1_x/m0,m1_y/m0);
		}
	}


	printf("\n");

	result = fopen("k_vector_xy.dat","w");
	for(k=0;k<18;k++){
		for(n=0;n<4;n++){
			fprintf(result,"%% #%d vector of photo #%d\n",n+1,k+1);
			fprintf(result,"%lf\t%lf\n", k_result[k][n][0], k_result[k][n][1]);
		}
	}
	fclose(result);

	return 0;
}