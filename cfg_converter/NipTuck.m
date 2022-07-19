% NipTuck.m %
% -- 09/30/2015 --  Yin Zhang%
% input cfgfile and choose atoms to form a circular dislocation line %
clc
clear all

cfgfile = ['relax.cfg'];

% cut -f 5 -d" "5" the fifth words %
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
% system ('rm testdata'); %
%command = ['grep -A 39632 "Si" pillar.cfg | tail -n +2 > testdata'];
% command = ['grep "Cu ." ',cfgfile,'  | cut -f 3-5 -d" " >  testdata']; %
 command = ['grep "63.546 Cu" ',cfgfile,'  | cut -f 3-5 -d" " >  testdata']; %
system (command); 

load testdata;
s1 = testdata(:,1);
s2 = testdata(:,2);
s3 = testdata(:,3);
%system ('rm testdata');

x  =  s1 * H(1,1)  +  s2 * H(2,1)  +  s3 * H(3,1);
y  =  s1 * H(1,2)  +  s2 * H(2,2)  +  s3 * H(3,2);
z  =  s1 * H(1,3)  +  s2 * H(2,3)  +  s3 * H(3,3);

xyz = [x y z];


center1 = 294;% layer 1
center2 = 595;% layer 2
r = 10000; % radius of dislocation line

h1 = [x( center1 + 1 ) y( center1 + 1 ) z( center1 + 1 )]; 
h2 = [x( center2 + 1 ) y( center2 + 1 ) z( center2 + 1 )];

%the center of the circular line
x1 = h1(1); y1 = h1(2); z1 = h1(3);
x2 = h2(1); y2 = h2(2); z2 = h2(3);

n = [-1.0887    1.3472         0];%the normal vector of slip plane

m = 1;

for i = 1:Natom
    i;
    
    dot_product = n(1)*(x(i) - x1) + n(2)*(y(i) - y1) + n(3)*(z(i) - z1);%In the plane?
    
    IN = z(i)-90 ;%Inside the circle?
    
    circle = (x(i)-x1)^2 + (y(i)-y1)^2 + (z(i)-z1)^2 - r^2;
    
    if abs(dot_product) < 0.17&& circle < 0% && IN > 0 
        layer(m) = i-1;
        m = m + 1;
    end   

end

fo = fopen(sprintf('%s','bccID1'),'w');

for i = 1:m-1
    fprintf( fo, '%g\n', layer(i)+1 );
end;

fclose(fo);

m = 1;

for i = 1:Natom
    i;
    
    dot_product = n(1)*(x(i) - x2) + n(2)*(y(i) - y2) + n(3)*(z(i) - z2);%In the plane?
    
    IN = z(i)-90 ;%Inside the circle?
    
    circle = (x(i)-x2)^2 + (y(i)-y2)^2 + (z(i)-z2)^2 - r^2;
    
    if abs(dot_product) < 0.17&& circle < 0% && IN > 0 
        layer(m) = i-1;
        m = m + 1;
    end   

end
    
    
fo = fopen(sprintf('%s','bccID2'),'w');

for i = 1:m-1
    fprintf( fo, '%g\n', layer(i)+1 );
end;

fclose(fo);

command = ['gcc NipTuck.c -o NipTuck'];
system(command);

command = ['./NipTuck relax.cfg'];
system(command);

%command = ['mpiexec -n 8 ./lmp_g++ < in.crack'];
%system(command);
