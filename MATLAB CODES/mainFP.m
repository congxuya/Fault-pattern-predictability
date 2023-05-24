%Check fault pattern predictability for bounded and unbounded LPNs
function [A,M,tf] = mainFP(name1,m1,n1,name2,m2,n2,L1,L2,SF,F)

%name1 .pnt file name of a plant net
%m1 the number of places of a plant net
%n1 the number of transitions of a plant net
%name2 .pnt file name of a fault pattern net
%m2 the number of places of a fault pattern net
%n2 the number of transitions of a fault pattern net
%L1 a cell array each row contains indices of observable transitions with the same label of a plant net
%L2 a cell array each row contains indices of observable transitions with the same label of a fault pattern net
%L2 has the same number of rows as L1, empty matrix reprents no observable transition
%SF the bijective function e.g., [3 1;3 2]   (t3,N1)  (t3,N2)
%F the index of the fault place

%A: basis fault patteern predictor graph
%first column of A: the index of current basis marking
%second column of A: the firing transition
%third column of A: the succeeding marking
%fourth column to the last column of A: the mimial e-vecotrs
%M: the basis markings in a fault pattern predictor graph
%Ad: the subgraph strating from dangerous basis marking

%tf=1: the LPN is predictable 
%tf=0: the LPN is not predictable

[Pre1,Post1,M01]=LY_pnt2NW3(name1,m1,n1);
[Pre2,Post2,M02]=LY_pnt2NW3(name2,m2,n2);
[Pre,Post,M0,L,T2,T1,Tf,Tf2] = FPPN(Pre1,Post1,M01,Pre2,Post2,M02,L1,L2,SF,F);
[A,M,Ad]=BFPPG(Pre,Post,M0,L,T1,T2,Tf,Tf2);
[tf]=CycleCheck2(A,Ad,M,Pre,Post,T1);

end