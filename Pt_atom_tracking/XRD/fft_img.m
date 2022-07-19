% clearvars;
% 
% load Sq.dat;

for i = 1:18
    i1 = 1001*(i-1)+1;
    in = 1001*i;
    Im = Sq(i1:in,:)/300;
    output = [num2str(i) '.jpg'];
    imwrite(Im,output);
end