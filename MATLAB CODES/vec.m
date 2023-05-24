
function [sv]=vec(v)

sv=[];
for i=1:length(v) % the length of a vector
    if i==1
        for j=0:v(i)
            sv=[sv;j,zeros(1,length(v)-1)];
        end
    else
        c=sv;
        for j=1:v(i)
            c=c+[zeros(size(c,1),i-1),ones(size(c,1),1),zeros(size(c,1),size(c,2)-i)];
            sv=[sv;c];
        end
    end
end