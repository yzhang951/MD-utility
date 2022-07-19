clearvars;

load Sq.dat;
a = 1;
i1 = 1001*(a-1)+1;in = 1001*a;
Im = Sq(i1:in,:);

Im = Im/(max(max(Im)));
n = length(Im);
m = 7;
f = waitbar(0,'Please wait ...');

peak=zeros(100,2);
peak_num = 0;
% Find peak diffraction points
for i=1:n        
        
    waitbar(i/n,f,'Processing data');
    for j=1:n
        max_flag = 1;
        m0 = 0; %  zero momentum
        m1 = [0 0]'; %  first momentum

        for k=-m:m
            for l=-m:m
                i_near = i+k;
                j_near = j+l;
                if(i_near<=0||j_near<=0)
                    continue;
                end
                if(i_near>n||j_near>n)
                    continue;
                end
                
                m0 = m0 + Im(i_near,j_near);
                m1 = m1 + [k;l]*Im(i_near,j_near);
                
                if(k==0&&l==0)
                    continue;
                end
                
                if(Im(i_near,j_near)>Im(i,j))
                    max_flag = -1;
                end
            end
        end
        
        if(max_flag>0&&Im(i,j)>0.02)
            peak(peak_num+1,:) = [i;j] + m1/m0;
            peak_num = peak_num + 1;
        end
    end
end

peak = peak(1:peak_num,:);
peak = peak - (n+1)/2;
% Find primitive vector
% Sort these k-points by their norm
sort_flag = -1;
while sort_flag<0
    sort_flag = 1;
    for i = 1:(peak_num-1)
        if(norm(peak(i,:))>norm(peak(i+1,:)))
            sort_flag = -1;
            xy_tmp = peak(i+1,:);
            peak(i+1,:) = peak(i,:);
            peak(i,:) = xy_tmp;
        end
    end
end

a1 = peak(2,:);
a2 = peak(4,:);
theta = acosd(a1*a2'/norm(a1)/norm(a2));

close(f);
plot(peak(:,2),peak(:,1),'rx');
