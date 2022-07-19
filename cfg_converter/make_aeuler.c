/*  Yin Zhang 16th, Aug
 *  generate random aeuler angle */
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

int main(){
	int 	i,n;
	double 	ph_min,ph_max,th_min,th_max,om_min,om_max;
	double	r,ph,th,om;
	char	output[50];
	FILE	*fp;

/*	Input parameter */

	printf("This is a program that generates euler angle for CPFE.\n");
	printf("Please input the name of euler file:\n");
	scanf("%s",output);
	printf("Please input the number of element:\n");
	scanf("%d",&n);
	printf("Please input the range of PHI,  as 'min max':\n");
	scanf("%lf%lf",&ph_min,&ph_max);
	printf("Please input the range of THETA,as 'min max':\n");
	scanf("%lf%lf",&th_min,&th_max);
	printf("Please input the range of OMEGA,as 'min max':\n");
	scanf("%lf%lf",&om_min,&om_max);

/*	Random Number	*/
	fp = fopen(output,"w");
	fprintf(fp,"KFLAG : 1 --- SAME CRYSTAL OR SET OF CRYSTALS AT ALL INTEGRATION POINTS\n        2 --- DIFFERENT CRYSTALS AT INTEGRATION POINTS\n2\nNCRYS\n1\nEULER ANGLES FOR EACH CRYSTAL AT AN INTEGRATION POINT\nPHI THETA OMEGA\n");
	for(i=0;i<n;i++){
		srand(time(NULL)+i);
		r = (rand()%10000)/10000.0;
		ph = ph_min + (ph_max - ph_min)*r;	
		srand(time(NULL)+i*i);
		r = (rand()%10000)/10000.0;
		th = th_min + (th_max - th_min)*r;	
		srand(time(NULL)+i*i*i);
		r = (rand()%10000)/10000.0;
		om = om_min + (om_max - om_min)*r;

		fprintf(fp,"%9g%9g%9g\n",ph,th,om);
	}
	fclose(fp);
	printf("Complete!\n");
	return 0;
}
