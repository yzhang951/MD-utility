clearvars;

filename0 = 'perfect.cfg';


% disl = [slip_norm slip_dir line_dir x y z] of dislocation
nu = 0.3;
% Unit vector in XYZ simulation box coord
Q(1,:) = [1 0 0];
Q(2,:) = [0 1 0];
Q(3,:) = [0 0 1];
inv_Q = inv(Q);
ex = inv_Q(1,:);
ey = inv_Q(2,:);
ez = inv_Q(3,:);

a0 = 3.147;
b_full = a0*sqrt(3)/2;
b_partial = b_full
d = 10;% stacking falut width

n_123 = (ex + ey)/sqrt(2);
b_1 = b_partial*(  ex  -  ey  +  ez)/sqrt(3);
b_2 = b_partial*(  ex  -  ey  -  ez)/sqrt(3);

n_456 = (ex - ey)/sqrt(2);
b_4 = b_partial*(  ex  +  ey  +  ez)/sqrt(3);
b_5 = b_partial*(  ex  +  ey  -  ez)/sqrt(3);



%system(command);
    
cfgfile = filename0;    
command = ['mul ',cfgfile,' 1 1 1 ',cfgfile];  
system(command);
   
command = ['grep "Number of particles = " ',cfgfile,' | cut -f 5 -d" " > log1'];
system (command);
  
load log1; 
Natom = log1; 
system ('rm log1');


   
command = ['grep "H0(1,1) =" ',cfgfile,' | cut -f 3 -d" " > log1'];  
system (command);   
load log1;
H(1,1) = log1; 
system ('rm log1');

    
command = ['grep "H0(1,2) =" ',cfgfile,' | cut -f 3 -d" " > log1'];
    
system (command);   
load log1;    
H(1,2) = log1;     
system ('rm log1');
    
    
command = ['grep "H0(1,3) =" ',cfgfile,' | cut -f 3 -d" " > log1'];
system (command);
load log1;
H(1,3) = log1; 
system ('rm log1');
    
command = ['grep "H0(2,1) =" ',cfgfile,' | cut -f 3 -d" " > log1'];
system (command);
load log1;
H(2,1) = log1; 
system ('rm log1');
    
command = ['grep "H0(2,2) =" ',cfgfile,' | cut -f 3 -d" " > log1'];
system (command);
load log1;
H(2,2) = log1; 
system ('rm log1');

command = ['grep "H0(2,3) =" ',cfgfile,' | cut -f 3 -d" " > log1'];
system (command);
load log1;
H(2,3) = log1; 
system ('rm log1');

command = ['grep "H0(3,1) =" ',cfgfile,' | cut -f 3 -d" " > log1'];
system (command);
load log1;
H(3,1) = log1; 
system ('rm log1');

command = ['grep "H0(3,2) =" ',cfgfile,' | cut -f 3 -d" " > log1'];
system (command);
load log1;
H(3,2) = log1; 
system ('rm log1');
    
command = ['grep "H0(3,3) =" ',cfgfile,' | cut -f 3 -d" " > log1'];
system (command);
load log1;
H(3,3) = log1; 
system ('rm log1');
    
    % break; %

command = ['grep "95.94 Mo" ',cfgfile,'  | cut -f 3-5 -d" " >  testdata']; %
system (command); 

load testdata;
s1 = testdata(:,1);
s2 = testdata(:,2);
s3 = testdata(:,3);
    %system ('rm testdata');
    
x  =  s1 * H(1,1)  +  s2 * H(2,1)  +  s3 * H(3,1);
y  =  s1 * H(1,2)  +  s2 * H(2,2)  +  s3 * H(3,2);
z  =  s1 * H(1,3)  +  s2 * H(2,3)  +  s3 * H(3,3);
id = 1:length(x);
xyz = [id' x y z];

output = ['Mo_disl.in'];
fp = fopen(output,'w');    
fprintf(fp,'Mo\n\n\t%d atoms\n',length(x));    
fprintf(fp,'\t0 bonds\n\t0 angles\n\t0 dihedrals\n\t0 impropers\n\n');    
fprintf(fp,'\t1 atom types\n\t0 bond types\n\t0 angle types\n\t0 dihedral types\n\t0 improper types\n\n');    
% fprintf(fp,'\t0 %g xlo xhi\n\t-100 %g ylo yhi\n\t0 %g zlo zhi\n\n',H(1,1),H(2,2)+100,H(3,3));
fprintf(fp,'\t-100 %g xlo xhi\n\t-100 %g ylo yhi\n\t-100 %g zlo zhi\n\n',H(1,1)+100,H(2,2)+100,H(3,3)+100);
fprintf(fp,'Atoms\n\n');

Lx = H(1,1);
Ly = H(2,2);
disl = [ 
%        n_123   b_1   b_1   75 70 H(3,3)/2;         
%        n_456  -b_5   b_5   75 75 H(3,3)/2;
        n_123   b_1   0 0 1   75 70 H(3,3)/2;         
        n_456   b_5   0 0 1  75 75 H(3,3)/2;
         ];


for i=1:length(x)
    uxyz = zeros(1,3);
    for j=1:length(disl(:,1))
% new coord in Norm-Line coord
        xyz_disl = disl(j,10:12);
        slip_norm = disl(j,1:3);
        Burgers = disl(j,4:6);

        line_dir = disl(j,7:9)/norm(disl(j,7:9));
        line_perp = cross(slip_norm, line_dir);

        b_screw = Burgers*line_dir';
        b_edge = Burgers*line_perp';
        
        xyz_bnl = [x(i) y(i) z(i)] - xyz_disl;
        x0 = xyz_bnl*line_perp';
        y0 = xyz_bnl*slip_norm';
        r  = sqrt(x0*x0+y0*y0);
        z0 = xyz_bnl*line_dir';
        theta = atan2(y0,x0) + pi;

        % Edge componet and screw componet in Norm-Line coord
        u = [ b_edge/2/pi*(theta + sin(2*theta)/4/(1-nu))
             -b_edge/2/pi*((1-2*nu)*log(r)/2/(1-2*nu) + cos(2*theta)/4/(1-nu))
              b_screw*theta/2/pi]';
        
        % Convert to simulation cell
        uxyz = uxyz + u*[line_perp' slip_norm' line_dir']';
    end

       
    fprintf(fp,'      %d   %d     %g %g %g   0  0  0\n',i, 1,x(i)+uxyz(1),y(i)+uxyz(2),z(i)+uxyz(3));    
    
end

fclose(fp);
