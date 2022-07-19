#include <stdio.h>
#include <math.h>
#include <string.h> 
#include <malloc.h>

int Skip_line(FILE *fid, int n)
{
    int i;
    for (i=0;i<n;i++)
	while(fgetc(fid)!='\n')
        continue;
    return 0;
}

int main(int argc, char **argv){
    char cfgfile[20],mass[10],element[5],output[20],command[50],temp[5],line[80];
    FILE *fp,*fo;
    int i,j,k=0,frame,Natom,ID_layer_1[500],n1=0,ID_layer_2[500],n2=0;
    double H[4][4],s1,s2,s3;
    struct atom {
           int id;
           double x,y,z;
    };
    struct atom *Cu;
    strcpy (cfgfile,argv[1]);
    strcpy (mass,argv[2]);
    strcpy (element,argv[3]);
    
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
    Cu=(struct atom*)malloc(sizeof(struct atom)*Natom);
    
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
                         fscanf(fp,"%lf%lf%lf",&s1,&s2,&s3);
                         (Cu+i)->x = s1 * H[1][1]  +  s2 * H[2][1]  +  s3 * H[3][1];
                         (Cu+i)->y = s1 * H[1][2]  +  s2 * H[2][2]  +  s3 * H[3][2];
                         (Cu+i)->z = s1 * H[1][3]  +  s2 * H[2][3]  +  s3 * H[3][3];
    	}
    fclose(fp);
    
/*  i = 0;
    fp = fopen("bccID1", "r");
    while(fgets(line, 80, fp) !=NULL){
	sscanf(line, "%d", &ID_layer_1[i]);
	i++;
	n1++;
    }
    fclose(fp);

    i = 0;
    fp = fopen("bccID2", "r");
    while(fgets(line, 80, fp) !=NULL){
	sscanf(line, "%d", &ID_layer_2[i]);
	i++;
	n2++;
    }
    fclose(fp);*/
    sprintf(output,"cfg.xyz");
    fo = fopen(output,"w");
/*  
    fprintf(fo, "Dimond Silicon\n");
    fprintf(fo, "\n");
  
    fprintf(fo, "      %d atoms\n", Natom);
    fprintf(fo, "      0 bonds\n");
    fprintf(fo, "      0 angles\n");
    fprintf(fo, "      0 dihedrals\n");
    fprintf(fo, "      0 impropers\n");
    fprintf(fo, "\n");
  
    fprintf(fo, "      %d atom types\n",1);
    fprintf(fo, "      0 bond types\n");
    fprintf(fo, "      0 angle types\n");
    fprintf(fo, "      0 dihedral types\n");
    fprintf(fo, "      0 improper types\n");
    fprintf(fo, "\n");
  
    fprintf(fo, "      0 %g xlo xhi\n",H[1][1]);
    fprintf(fo, "      0 %g ylo yhi\n",H[2][2]);
    fprintf(fo, "      0 %g zlo zhi\n",H[3][3]);
    fprintf(fo, "\n");
    fprintf(fo, "Atoms\n");
    fprintf(fo, "\n");
*/
/*  for(j=0;j<n1;j++){
	i = ID_layer_1[j]-1;
	fprintf(fo,"      %d   %d     %g %g %g   0  0  0\n",j+1,1,(Cu+i)->x,(Cu+i)->y,(Cu+i)->z);
	}
    for(j=0;j<n2;j++){
	i = ID_layer_2[j]-1;
	fprintf(fo,"      %d   %d     %g %g %g   0  0  0\n",n1+j+1,1,(Cu+i)->x,(Cu+i)->y,(Cu+i)->z);
	}*/
    for(i=0;i<Natom;i++) {
		frame = 0;
/*		for(j=0;j<n1;j++)
			if(ID_layer_1[j] == i+1){
				k++;frame++;
				continue;
				}
		for(j=0;j<n2;j++)
			if(ID_layer_2[j] == i+1){
				k++;frame++;
				continue;
				}*/
			if(frame !=0) continue;
                         	fprintf(fo,"%s %g %g %g\n",element,(Cu+i)->x,(Cu+i)->y,(Cu+i)->z);
                         }
    fclose(fo);
    
    return 0;
}
