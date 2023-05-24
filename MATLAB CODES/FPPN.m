%Generated the fault pattern predictor net
function [Pre,Post,M0,L,T2,T1,Tf,Tf2] = FPPN(Pre1,Post1,M01,Pre2,Post2,M02,L1,L2,SF,F)

%Pre1 the pre-incidence matrix of a plant net
%Post1 the post-incidence matrix of a plant net
%M01 the initial marking of the a plant net
%L1 a cell array each row contains indices of observable transitions with the same label of the plant net
%Pre2 the pre-incidence matrix of a fault pattern net
%Post2 the post-incidence matrix of a fault pattern net
%M02 the initial marking of a fault pattern net
%L2 a cell array each row contains indices of observable transitions with the same label of a fault pattern net 
%L2 has the same number of rows as L1, empty matrix reprents no observable transition
%SF the bijective function [3 1;3 2]   (t3,N1)  (t3,N2)
%F  the index of the fault place

%Pref the pre-incidence matrix of a fault pattern labeled net
%Postf the post-incidence matrix of a fault pattern labeled net
%M0f the initial marking of a fault pattern labeled net
%Lf the matrix stroing labels of a fault pattern labeled net

m1=size(Pre1,1);   %the number of places in a plant net
n1=size(Post1,2);   %the number of transition in a plant net
m2=size(Pre2,1);  %the number of places in a fault pattern net
n2=size(Post2,2);  %the number of transitions in a fault pattern net
mf=m1+m2;           %the number of places in a fault pattern labeled net
Tsync=unique(SF(:,1))';  %the synchronization transitions

Pref=[Pre1,zeros(m1,n2);zeros(m2,n1),Pre2];
Postf=[Post1,zeros(m1,n2);zeros(m2,n1),Post2];
M0f=[M01,M02];

Lm1=zeros(size(L1,1),n1);
for i=1:size(Lm1,1)
    if ~isempty(L1{i})
    Lm1(i,L1{i})=1;
    end
end
Lm2=zeros(size(L2,1),n2);
for i=1:size(Lm2,1)
    if ~isempty(L2{i})
    Lm2(i,L2{i})=1;
    end
end
Lf=[Lm1,Lm2];

for i=1:n2
    Pref(1:m1,n1+i)=Pre1(:,SF(i,1));
    Postf(1:m1,n1+i)=Post1(:,SF(i,1));
end
Pref(m1+F,Tsync)=1;
Postf(m1+F,Tsync)=1;

%construct nonfailure subnet of fault pattern labeled net (n-FPLN)
Pren=Pref;
Postn=Postf;
Nf=Postf-Pref;
Ln=Lf;
Pren(:,[Tsync,find(Nf(m1+F,:)>0)])=[];
Postn(:,[Tsync,find(Nf(m1+F,:)>0)])=[];
Ln(:,[Tsync,find(Nf(m1+F,:)>0)])=[];
M0=[M0f,M0f];

ne=size(Ln,1);   %the number of events in n-FPLN
m=2*mf;  %the number of places of fault pattern predictor net
n=0;
for i=1:ne
    n=n+sum(Ln(i,:))^2;
end
nun=size(Ln,2)-sum(sum(Ln));%the number of unobservable transitions in n-FPLN
n=n+2*nun;
n=n+size(Tsync,2)+size(find(Nf(m1+F,:)>0),2);
n=n+size(Ln,2);

%Pre the pre-incidence matrix of a fault pattern predictor net (FPPN)
%Post the post-incidence matrix of an FPPN
%M0 the initial marking of an FPPN
%L the matrix storing labels of an FPPN

Pre=zeros(m,n);
Post=zeros(m,n);
L=zeros(2*ne,n);     
r=1;
no=[];    
for i=1:ne
    a=find(Ln(i,:)>0);   %the observable transitions related to event i of n-FPLN
    no=[no,a];       %the observable transitions of n-FPLN
    L(i,r:r+size(a,2)^2-1)=1;
    for j=1:size(a,2)
        for k=1:size(a,2)
            Pre(:,r)=[Pren(:,a(j));Pren(:,a(k))];
            Post(:,r)=[Postn(:,a(j));Postn(:,a(k))];
            r=r+1;
        end
    end    
end
nu=setdiff(1:size(Ln,2),no); %the unobservable trnasitions of n-FPLN
Prenu=Pren(:,nu); %the pre-incidence matrix of unobservable subnet of n-FPLN
Postnu=Postn(:,nu); %the post-incidence matrix of unobservable subnet of n-FPLN

for i=1:nun
    Pre(1:mf,r)=Prenu(:,i);
    Post(1:mf,r)=Postnu(:,i);
    r=r+1;
end

for i=1:nun
    Pre(mf+1:2*mf,r)=Prenu(:,i);
    Post(mf+1:2*mf,r)=Postnu(:,i);
    r=r+1;
end

Tf=[];
for i=1:size(Tsync,2)
    Pre(1:mf,r)=Pref(:,Tsync(i));
    Post(1:mf,r)=Postf(:,Tsync(i));
    Tf=[Tf,r];
    r=r+1;
end

a=find(Nf(m1+F,:)>0);
Tf2=[];
for i=1:size(a,2)
    Pre(1:mf,r)=Pref(:,a(i));
    Post(1:mf,r)=Postf(:,a(i));
    Tf=[Tf,r];
    Tf2=[Tf2,r];
    r=r+1;
end
T2=1:r-1;

for i=1:size(Pren,2)
    Pre(mf+1:2*mf,r)=Pren(:,i);
    Post(mf+1:2*mf,r)=Postn(:,i);
    r=r+1;
end
T1=size(T2,2)+1:r-1;
L(ne+1:2*ne,n-size(Ln,2)+1:n)=Ln;

end



