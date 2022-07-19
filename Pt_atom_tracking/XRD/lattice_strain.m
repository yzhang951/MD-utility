clear;
clf;
load Sq.dat;

%load fig2_full_Sq.dat;
%Sq = fig2_full_Sq;

q = 1.5:0.001:5.5-0.001;
hold on;
for i = 1:length(Sq(:,1))
    S = Sq(i,:);
    f = fit(q',S','gauss1');
    d1(i) = 1/(f.b1);
%     d2(i) = 1/(f.b2);
%    d(i) = (d1(ils)*f.a1 + d2(i)*f.a2)/(f.a1+f.a2);
    d(i) = d1(i);
    if(mod(i,10)==1)
        j = (i-1)/10+1;
        subplot(6,1,j);
        s = num2str(i);    
        plot(q,S);
        hold on;
        title(s);           
        xlim([4,4.8]);
        ylim([0,30]);
    end 
end

figure;
% t = [0 0.0167 0.033 0.05 0.0667 0.083];
plot(2*pi*d,'-ob');
x = 2*pi*d;
save fft.dat x -ascii;