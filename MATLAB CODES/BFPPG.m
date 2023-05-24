%Generate the basis fault pattern predictor graph
function [A,M,Ad]=BFPPG(Pre,Post,M0,L,T1,T2,Tf,Tf2)
%A: basis fault patteern predictor graph
%first column of A: the index of current basis marking
%second column of A: the firing transition
%third column of A: the succeeding marking
%fourth column to the last column of A: the mimial e-vecotrs
%M: the basis markings in the fault pattern predictor graph
%Ad: the subgraph strating from dangerous basis marking

%L1 a cell array each row contains indices of observable transitions with the same label of the plant net
%MO row vector
%TI the set of index of unobservable transitions (row vector)
%T1 the set of transitions of the T1-induced subnet of FPPN (row vector)
%T2 the set of transitions in T2-induced subnet of FPPN    (row vector)
%the index of T2 is before that of T1
%Tf the set of transitions before the fault place F (row vector)
%Tf2 the set of boundary transitions 

tpre=Pre';
tpos=Post';
m=size(M0,2);  %the number of places

%E2 is a row vector if ti is observable , then E2(i)=1 
E2=zeros(1,size(L(:,T2),2));
for i=1:size(L,1)
    E2=E2+L(i,T2);
end
N=tpos-tpre;
tpos2=tpos(T2,:);
tpre2=tpre(T2,:);
N2=tpos2-tpre2;     %the matrix of the T2-induced subnet
[tnum2,~]=size(N2);   %the number of transitions in the T2-induced subnet
[tnum,~]=size(N);   %the number of transitions in the FPPN
%Efx records the unobservable transitions in nonfailure subnet of T2-induced subnet
Efx=ones(1,tnum);      
Efx(:,setdiff(setdiff(T2,Tf),find(E2>0)))=0;
tnumx=size(find(Efx==0),2);  %the number of unobservable transitions in nonfailure subnet of T2-induced subnet
Nx=N(Efx==0,:);
M=M0;
Mnew=M0;
A=[];
Ad=[];
a=1;
b=0;
c=0;  %c=0: M is not a dangerous basis marking,c=1: M is a dangerous basis marking
Btemp=[];  %store the set of previously generated basis markings and the minimal e-vectors
while a<=size(Mnew,1)
    for i=1:size(Tf2,2)
        Mbuf=Mnew(a,:)-tpre(Tf2(i),:);
        us=zeros(1,tnumx);
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
                        l=find(Nx(:,k)>0);
                        Mbuf=[Mbuf;bsxfun(@plus,Nx(l,:),Mbuf(j,:))];
                        us=[us;repmat(us(j,:),[size(l,1),1])];
                        e=eye(tnumx);
                        us(size(us,1)-size(l,1)+1:size(us,1),:)=us(size(us,1)-size(l,1)+1:size(us,1),:)+e(l,:);
                        Mbuf(j,:)=[];
                        us(j,:)=[];
                        j=j-1;
                        break
                    end
                end
                j=j+1;
            end
        end
        if ~isempty(us)
            c=1;
            break;
        end
    end            
    if c~=1
        for i=1:tnum2
            if E2(i)==1      %for each observable transition compute e-vector of current basis marking M(a,:)
                Mbuf=Mnew(a,:)-tpre2(i,:);
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
                                l=intersect(find(Efx==0)',find(N2(:,k)>0));
                                Mbuf=[Mbuf;bsxfun(@plus,N2(l,:),Mbuf(j,:))];
                                us=[us;repmat(us(j,:),[size(l,1),1])];
                                e=eye(tnum);
                                us(size(us,1)-size(l,1)+1:size(us,1),:)=us(size(us,1)-size(l,1)+1:size(us,1),:)+e(l,:);
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
                us1=us;
                %compute the minimal e-vectors
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
                Amin=[Amin;bsxfun(@plus,Mbuf(Bmin==1,:),tpos2(i,:))];
                Btemp=[Btemp;repmat(Mnew(a,:),[size(find(Bmin==1),2),1]),us(Bmin==1,:)];

                M2=Btemp(:,1:m);
                for x=1:size(Btemp,1)
                    M2=[M2;Btemp(x,1:m)+Btemp(x,m+1:m+tnum2)*N2];
                    d2=vec(Btemp(x,m+1:m+tnum));
                    dx2=bsxfun(@plus,d2*N,Btemp(x,1:m));
                    M2=[M2;dx2(all(dx2>=zeros(size(dx2,1),1),2),:)];
                end
                M2=unique(M2,'rows');

                M1=Mnew(a,:);
                for q=1:size(Amin,1)
                    M1=[M1;Amin(q,:)];
                    d=vec(us(q,:));
                    dx=bsxfun(@plus,d*N,Mnew(a,:));
                    M1=[M1;dx(all(dx>=zeros(size(dx,1),1),2),:)];
                    M1=unique(M1,'rows');
                    temp=repmat(M1,[size(M2,1),1])-repmat(M2,[size(M1,1),1]);
                    temp=unique(temp,'rows');
                    for z=1:size(temp,1)
                        tempz=temp(z,:);
                        if isempty(find(tempz<0, 1))
                            tempx=temp(z,:)>0;
                            tempz=tempz(~isnan(tempz));
                            if all (tempz>=0)
                                Amin(q,tempx)=Inf;
                            end
                        end
                    end
                end
                for j=1:size(Amin,1)
                    [~,k]=ismember(Amin(j,:),Mnew,'rows');
                    if k~=0
                        A=[A;a,i,k,us1(j,:)];
                    else
                        Mnew=[Mnew;Amin(j,:)];
                        A=[A;a,i,max(size(Mnew,1),size(M,1)),us1(j,:)];
                    end
                end
            end
        end
        a=a+1;
    else
        M01=Mnew(a,:);
        [Ax,Mx]=BSONNODES(tpre,tpos,M01,L,T1,T2);
        if size(M,1)>=size(Mnew,1)
            M=[M;Mx(2:size(Mx,1),:)];
        else
            M=[Mnew;Mx(2:size(Mx,1),:)];
        end
        if ~isempty(Ax)
            for i=1:size(Ax,1)
                Ax(i,1)=find(ismember(M,Mx(Ax(i,1),:),'rows')==1,1,'last');
                Ax(i,3)=find(ismember(M,Mx(Ax(i,3),:),'rows')==1, 1, 'last' );
            end
            Ad=[Ad;Ax];
        end
        A=[A;Ax];
        a=a+1;
    end    
end
A=unique(A,'rows');
Ad=unique(Ad,'rows');
if isempty(Ad)
    M=Mnew;
end
end
