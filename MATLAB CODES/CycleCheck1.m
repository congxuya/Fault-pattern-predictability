function [tf] = CycleCheck1(Ad)
%Check if there exists a cycle reachable from dangerous basis marking for bounded LPN
%Ad: the subgraph of BFPPG starting from dangerous basis marking
%tf=1: the bounded LPN is fault predictable w.r.t. a fault pattern
%tf=0: the bounded LPN is not fault predictable w.r.t. a fault pattern


if ~isempty(Ad)
    s=Ad(:,1)';
    t=Ad(:,3)';
    w=Ad(:,2)';
    G=digraph(s,t,w);
    tf=isdag(G);
else
    tf=1;
end
end

