function [Ax,Mx]=BSONNODES(tpre,tpos,M01,L,T1,T2)
%the function of basis sonnodes
%first column: the index of current basis marking
%second column: the firing transition
%third column: the succeeding basis marking
%fourth column: the feasible event at the current basis marking

%M0 row vector
%T1 the set of transitions of the T1-induced subnet
%T2 the set of transitions of the T2-induced subnet

tpre1=tpre(T1,:);
tpos1=tpos(T1,:);
m=size(M01,2);  %the number of places
n2=size(T2,2);   %the number of transitions in the T2-induced subnet
n=size(tpos,1);  %the number of transitions in a fault pattern predictor net
N=tpos-tpre;

E1=zeros(1,size(L(:,T1),2));
for i=1:size(L,1)
    E1=E1+L(i,T1);
end
N1=tpos1-tpre1;
tnum=size(N,1);
[tnum1,~]=size(N1);
Mx=M01;
Ax=[];
a=1;
b=0;
Btemp1=[]; 
while a<=size(Mx,1)
        for i=1:tnum1
            if E1(i)==1
                Mbuf=Mx(a,:)-tpre1(i,:);
                us=zeros(1,tnum);
                while 1
                    if size(Mbuf,1)==0
                        break
                    end
                    if Mbuf>=0
                        break
                    end
                    j=1;
                    while j<=size(Mbuf,1)
                        for k=1:size(Mbuf,2)
                            if Mbuf(j,k)<0
%                                 l=intersect(find(E1==0)',find(N1(:,k)>0));
%                                 Mbuf=[Mbuf;bsxfun(@plus,N1(l,:),Mbuf(j,:))];
%                                 us=[us;repmat(us(j,:),[size(l,1),1])];
%                                 e=eye(tnum);
%                                 us(size(us,1)-size(l,1)+1:size(us,1),:)=us(size(us,1)-size(l,1)+1:size(us,1),:)+e(l,:);
                                for l=1:tnum1
                                    if N1(l,k)>0 && E1(l)==0 
                                        Mbuf=[Mbuf;Mbuf(j,:)+N1(l,:)];
                                        us=[us;us(j,:)];
                                        us(size(us,1),l)=us(size(us,1),l)+1;
                                    end
                                end
                                Mbuf(j,:)=[];
                                us(j,:)=[];
                                j=j-1;
                                break
                            end
                        end
                        j=j+1;
                    end
                end
                if size(Mbuf,1)==0
                    continue
                end
                %compute the minimal e-vectors of the T1-induced subnet
                us=[zeros(size(us,1),n2),us];
                us(:,n+1:n+n2)=[];
                Bmin=zeros(1,size(Mbuf,1));
                for j1=1:size(us,1)
                    for j2=1:size(us,1)
                        if j1==j2
                            continue
                        end
                        if us(j1,:)>=us(j2,:)
                            b=1;
                            break
                        end
                    end
                    if b==1
                        b=0;
                    else
                        Bmin(j1)=1;
                    end
                end
                 %Amin is the set of successor basis markings of M(a,:) by firing ti 
                Amin=[];
                for j=1:size(Mbuf,1)
                    if Bmin(j)==1
                        Amin=[Amin;Mbuf(j,:)+tpos1(i,:)];
                        Btemp1=[Btemp1;Mx(a,:),us(j,:)];
                    end
                end
                
                M2=Btemp1(:,1:m);
               for x=1:size(Btemp1,1)
                   M2=[M2;Btemp1(x,1:m)+Btemp1(x,m+n2+1:m+n)*N1];
                   d2=vec(Btemp1(x,m+1:m+n));      %find all the vectors that are less than a minimal e-vector
                   dx2=bsxfun(@plus,d2*N,Btemp1(x,1:m));
                   M2=[M2;dx2(all(dx2>=zeros(size(dx2,1),1),2),:)];
               end
               M2=unique(M2,'rows');
               
               %assign Inf to unbounded places
               M1=Mx(a,:);
               for q=1:size(Amin,1)
                    M1=[M1;Amin(q,:)];
                    d=vec(us(q,:));
                    dx=bsxfun(@plus,d*N,Mx(a,:));
                    M1=[M1;dx(all(dx>=zeros(size(dx,1),1),2),:)];
                    M1=unique(M1,'rows');
                    temp=repmat(M1,[size(M2,1),1])-repmat(M2,[size(M1,1),1]);
                    temp=unique(temp,'rows');
                    for z=1:size(temp,1)
                        tempx=temp(z,:)>0;
                        tempz=temp(z,:);
                        tempz=tempz(~isnan(tempz));
                        if all (tempz>=0)
                            Amin(q,tempx)=Inf;
                        end
                    end
                end
               
               for j=1:size(Amin,1)
                   for p=1:size(M1,1)
                       for v=1:size(M2,1)
                           temp=M1(p,:)-M2(v,:);
                           tempx=temp>0;
                           temp=temp(~isnan(temp));
                           if all (temp>=0)
                               Amin(j,tempx)=Inf;
                           end
                        end
                   end
               end
               for j=1:size(Amin,1)
                    [~,k]=ismember(Amin(j,:),Mx,'rows');
                    if k~=0
                        Ax=[Ax;a,i+n2,k,us(j,:)];
                    else
                        Mx=[Mx;Amin(j,:)];
                        Ax=[Ax;a,i+n2,size(Mx,1),us(j,:)];
                    end
                end
            end
        end
        a=a+1;
end
end


