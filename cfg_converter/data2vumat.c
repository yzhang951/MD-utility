/* Yin Zhang 26th Aug, 2016 */
#include <stdio.h>
#include <stdlib.h>

int main(){
	FILE	*fp_in,*fp_out;
	int		i,j,k,NOEL=130407;
	float	c11,c12,c44,gdoto,am,s0,ah0,ssat,rhard,ql;
	char	in[10]="para_inp",out[10]="test",temp[100];
	fp_in = fopen(in,"r");
	fp_out= fopen(out,"w");
	
	fprintf(fp_out,"    REAL DIMENSION(%d) :: YIN = (/ &\n",NOEL*10);

	fgets(temp,100,fp_in);
	fgets(temp,100,fp_in);
	for(i=0;i<NOEL;i++){
		fscanf(fp_in,"%f%f%f%f%f%f%f%f%f%f",&c11,&c12,&c44,&gdoto,&am,&s0,&ah0,&ssat,&rhard,&ql);
		fprintf(fp_out,"\t%g, %g, %g, %g, %g, &\n",c11,c12,c44,gdoto,am);
		fprintf(fp_out,"\t%g, %g, %g, %g, %g",s0,ah0,ssat,rhard,ql);
		if(i<NOEL-1) fprintf(fp_out,", &\n");
		else	fprintf(fp_out,"/)\n");
	}
	return 	0;
}
