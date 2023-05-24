%transform the pre, post incidence matrices and the initial marking to a .net file
function matrix_net(Pre,Post,M0)

filID=fopen('example.txt','w+');
fprintf(filID,'net example\n');

NumsTr=size(Pre,2);
NumsP=size(Pre,1);
NumsM=size(M0,2);
for i=1:NumsTr
    fprintf(filID,'tr t%d ',i-1);
    for j=1:NumsP
        if(Pre(j,i)==1)
           fprintf(filID,'p%d ',j-1); 
        end
        if(Pre(j,i)~=0 && Pre(j,i)~=1)
            fprintf(filID,'p%d*%d ',j-1,Pre(j,i));
        end
    end
    fprintf(filID,'-> ');
    for k=1:NumsP
       if(Post(k,i)==1)
          fprintf(filID,'p%d ',k-1); 
       end
       if(Post(k,i)~=0 && Post(k,i)~=1)
           fprintf(filID,'p%d*%d ',k-1,Post(k,i));
       end
    end
    fprintf(filID,'\n');
end
for pl=1:NumsM
   if(M0(1,pl)~=0) 
       fprintf(filID,'pl p%d (%d)\n',pl-1,M0(1,pl));
   end
end
fclose(filID);
end