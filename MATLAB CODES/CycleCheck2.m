function [tf] = CycleCheck2(A,Ad,M,Pre,Post,T1)
%Check if there exists a cycle reachable from DBM for unbounded LPN
%A: the BFPPG 
%Ad: the subgraph of BFPPG starting from dangerous basis marking (DBM)
%Pre: the pre-incidence matrix of a fault pattern predictor net
%Post: the post-incidence matrix of a fault pattern predictor net (FPPN)

%tf=1: the unbounded LPN is predictable 
%tf=0: the unbounded LPN is not predictable

tf=1;

if ~isempty(Ad)
    tpre1=Pre(:,T1)';
    tpos1=Post(:,T1)';
    N1=Post(:,T1)-Pre(:,T1);
    m1=size(N1,1);  %the number of places in T1-induced subnet of FPPN
    n1=size(T1,2);   %the number of transitions in T1-induced subnet of FPPN
    n2=size(Pre,2)-n1;
    G=digraph(A(:,1)',A(:,3)');
    [bins,~]=conncomp(G);
    sd=unique(Ad(:,1)');
    td=unique(Ad(:,3)');
    vd=unique([sd,td]);   %all the nodes reachable from DBM
    binsd=bins(vd);     %the strongly connected component reachable from DBM
    a=unique(binsd);
    for i=1:size(a,2)
        v=find(bins==a(i));  %the indexs of nodes in a(i)th SCC
        SCC=A(intersect(find(ismember(A(:,1),v)==1),find(ismember(A(:,3),v)==1)),:);   %the a(i)th SCC
        GSCC=(digraph(SCC(:,1)',SCC(:,3)'));
        if ~isdag(GSCC)  %construct the PN related to the a(i)th SCC
            ys=zeros(1,n1);
            mc=M(SCC(1,1),:);
            for j=1:size(SCC,1)
                y=SCC(j,size(SCC,2)-n1+1:size(SCC,2));   %the minimal e-vector associated to (M,t) in T1-induced subnet
                y1=zeros(1,n1);
                y1(SCC(j,2)-n2)=1;
                y=y+y1;
                ys=ys+y;
            end
            [Ac,~]=Unfolding(mc,ys,tpre1,tpos1);
            Gc=digraph(Ac(:,1)',Ac(:,3)');
            ms=unique(Ac(:,1));
            ns=unique(Ac(:,2));
            ms1=size(ms,1);  %the number of places associted with the a(i)th SCC     
            ns1=size(ns,1);  %the number of transitions in the PN associated with the a(i)th SCC
            Pres=zeros(ms1,ns1);
            Posts=zeros(ms1,ns1);
            for k=1:size(Ac,1) 
                Pres(Ac(k,1),Ac(k,2)==ns)=1;
                Posts(Ac(k,3),Ac(k,2)==ns)=1;
            end
            %solve the LPP to check if there exists a repetitive cycle
            Aeq1=[zeros(ms1,n1),Posts-Pres];
            Aeq2=zeros(ns1,n1);
            for p=1:ns1
                Aeq2(p,ns(p))=1;
            end
            Aeq2=[Aeq2,-eye(ns1)];
            Aeq=[Aeq1;Aeq2];
            beq=zeros(ms1+ns1,1);
            An=[-N1,zeros(m1,ns1)];
            An=[An;zeros(1,n1),-ones(1,ns1)];
            bn=[zeros(m1,1);-1];
            %intcon=1:ns+n1;
            lb=zeros(ns1+n1,1);
            ub=Inf*ones(ns1+n1,1);
            f=[1;zeros(ns1+n1-1,1)];
            [x,~,flag,~]=linprog(f,An,bn,Aeq,beq,lb,ub);
            if flag==1
                [~,edgecycles]=allcycles(Gc,'MinCycleLength',sum(x(n1+1:n1+ns1)),'MaxCycleLength',sum(x(n1+1:n1+ns1)));  
                for q=1:size(edgecycles,1)
                    if size(find(x(n1+1:n1+ns1)==1),2)==size(edgecycles{q},2)
                        if all(find(x(n1+1:n1+ns1)==1)==edgecycles{q})
                            tf=0;
                            break;
                        end
                    end
                end    
            end
            if tf==0
                break;
            end
        end
    end
end

