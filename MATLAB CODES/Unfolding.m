function [Ac,Mc] = Unfolding(mc,y,tpre1,tpos1)
%mc the basis marking in SCC reachable from a dangerous basis marking
%y the firing vector at mc
%N1 the transpose of incidence matrix of T1-induced subnet

N1=tpos1-tpre1;
Ac=[];
a=1;   %the index of mc
b=0;   %label the considered marking
Mc=[mc,y];
m1=size(N1,2); %the number of places in T1-induced subnet

while a<=size(Mc,1)
    for i=1:size(y,2)
        if Mc(a,m1+i)>0&&all(Mc(a,1:m1)-tpre1(i,:)>=0)
            %Amin is the successor node after firing ti in y at Mc(a,:)
            Amin=[Mc(a,1:m1)+N1(i,:),Mc(a,m1+1:size(Mc,2))];
            Amin(m1+i)=Amin(m1+i)-1; 
            %put the markings in Amin into the set Mc
            for k=1:size(Mc,1)
                temp=Amin(1:m1)-Mc(k,1:m1);
                tempx=temp>0;
                temp=temp(~isnan(temp)); %delete the NaN element in the vector temp
                if all (temp>=0)
                    Amin(tempx)=Inf;
                end
            end
            for j=1:size(Mc,1)
                if Amin(:,1:m1)==Mc(j,1:m1)
                    Ac=[Ac;a,i,j];
                    b=1;   %add the label ''old'' to a node
                    break
                end
            end
            if b==1
                b=0;
            else
                Mc=[Mc;Amin];
                Mc=unique(Mc,'row','stable');
                Ac=[Ac;a,i,size(Mc,1)];
            end         
        end
    end
    a=a+1;
end
Mc=Mc(:,1:m1);
end

