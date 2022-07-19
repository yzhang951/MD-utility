clearvars;

Ni_jpg = imread('Ni0.jpg');
Co_jpg = imread('Co0.jpg');
Cr_jpg = imread('Cr0.jpg');

Ni = double(Ni_jpg(:,:,2));
Co = double(Co_jpg(:,:,2));
Cr = double(Cr_jpg(:,:,2));

ROW = length(Ni(:,1));
COL = length(Ni(1,:));
aveNi = sum(sum(Ni))/ROW/COL;
aveCo = sum(sum(Co))/ROW/COL;
aveCr = sum(sum(Cr))/ROW/COL;
aveNiCoCr = aveNi+aveCo+aveCr;

Ni    = Ni/aveNiCoCr;
Co    = Co/aveNiCoCr;
Cr    = Cr/aveNiCoCr;

rmax = 800;    % unit pixel
rbin = 50;       % unit pixel
R    = 0:rbin:rmax;

S2Ni   = zeros(length(R),1);
numNi  = zeros(length(R),1);
S2Co   = zeros(length(R),1);
numCo  = zeros(length(R),1);
S2Cr   = zeros(length(R),1);
numCr  = zeros(length(R),1);





pixel2A = 106.692/2000;
theta = 0;
%theta = atan(sqrt(2));
% along 110 direction (x axis)
for row=1:ROW
    for col=1:COL
        for r=1:rmax
            rx = int32(r*cos(theta));
            ry = int32(r*sin(theta));
            if(col+rx>COL)
                continue;
            end
            if(row+ry>ROW)
                continue;
            end
            i = int32(r/rbin) + 1;
            if(i>length(R))
                continue;
            end
            numNi(i) = numNi(i) + 1;
            S2Ni(i)  = S2Ni(i) + Ni(row,col)*Ni(row+ry,col+rx);
            
            numCo(i) = numCo(i) + 1;
            S2Co(i)  = S2Co(i) + Co(row,col)*Co(row+ry,col+rx);
                        
            numCr(i) = numCr(i) + 1;
            S2Cr(i)  = S2Cr(i) + Cr(row,col)*Cr(row+ry,col+rx);
        end
    end
end

Ni_result = S2Ni.*(numNi.^(-1))-(sum(sum(Ni))/ROW/COL)^2;
Co_result = S2Co.*(numCo.^(-1))-(sum(sum(Co))/ROW/COL)^2;
Cr_result = S2Cr.*(numCr.^(-1))-(sum(sum(Cr))/ROW/COL)^2;
result = [R'*pixel2A Ni_result Co_result Cr_result];

clf;
hold on;
plot(R*pixel2A, Ni_result);
plot(R*pixel2A, Co_result);
plot(R*pixel2A, Cr_result);