clearvars;
% Seed of random number generator
rng(1245213);
% Element Num
ele_num = 5;
% Cell number in xyz directions
nx = 30;
ny = 3;
nz = 3;
% Wavefactor k
k_norm = 1:0.01:10;
kn = [1 1 0];
kn = kn/norm(kn);
% Fourier transformation for each element
Re = zeros(ele_num+1,length(k_norm));
Im = zeros(ele_num+1,length(k_norm));

% FCC lattice, xyz= [x y z Element_type]
% 1:Fe, 2:Ni, 3:Co, 4:Cr, 5:Mn
% Lattice constant, a = 1.0
a0 = 2.86;
xyz = zeros(nx*ny*nz*4,4);
i_atom = 0;
for i=1:nx
    for j=1:ny
        for k=1:nz
%          [0,0,0] atom          
           xyz(i_atom+1,1)=i*a0;
           xyz(i_atom+1,2)=j*a0;
           xyz(i_atom+1,3)=k*a0;
           n = floor(rand*ele_num)+1;
           xyz(i_atom+1,4) = n;

%          [0.5,0.5,0] atom           
           xyz(i_atom+2,1)=i*a0+0.5;
           xyz(i_atom+2,2)=j*a0+0.5;
           xyz(i_atom+2,3)=k*a0;
           n = floor(rand*ele_num)+1;
           xyz(i_atom+2,4) = n;
           
%          [0.5,0,0.5] atom           
           xyz(i_atom+3,1)=i*a0+0.5;
           xyz(i_atom+3,2)=j*a0;
           xyz(i_atom+3,3)=k*a0+0.5;
           n = floor(rand*ele_num)+1;
           xyz(i_atom+3,4) = n;
           
%          [0,0.5,0.5] atom           
           xyz(i_atom+4,1)=i*a0;
           xyz(i_atom+4,2)=j*a0+0.5;
           xyz(i_atom+4,3)=k*a0+0.5;
           n = floor(rand*ele_num)+1;
           xyz(i_atom+4,4) = n;
                      
           i_atom = i_atom+4;
        end
    end
end

for i=1:length(k_norm)
    k_i = k_norm(i)*kn;
    for i_atom=1:length(xyz(:,1))
        n = xyz(i_atom,4);
        Re(n,i) = Re(n,i) + cos(k_i*xyz(i_atom,1:3)');
        Im(n,i) = Im(n,i) - sin(k_i*xyz(i_atom,1:3)');
        Re(ele_num+1,i) = Re(ele_num+1,i) + cos(k_i*xyz(i_atom,1:3)');
        Im(ele_num+1,i) = Im(ele_num+1,i) - sin(k_i*xyz(i_atom,1:3)');
    end
end

i = 1; j = 2;
X_Re = Re(i,:).*Re(j,:)+Im(i,:).*Im(j,:);
X_Im = Im(i,:).*Re(j,:)-Im(j,:).*Re(i,:);
X2 = X_Re.^2+X_Im.^2;
plot(k_norm,X2);