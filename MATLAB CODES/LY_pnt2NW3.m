%transform the .pnt file to Pre and Post incidence matrices, and M0
function[Pre,Post,M0]=LY_pnt2NW3(name,m,n)
%name the name of the .pnt file
%m the number of places
%n the number of transitions
fid=fopen([name],'r');
line0=fgetl(fid); %skip the first line of the file
Post=zeros(m,n);
Pre=zeros(m,n);
line=fgetl(fid);%read the data of the second line of the file
M0=[];
k=1;
while(1)
clear L1 L2 c prestr poststr
c=strfind(line,',');%find the column place in each line
prestr(1:c-1)=line(1:c-1);%assign the characters before the comma to prestr
poststr(1:length(line)-c)=line(c+1:length(line));%assign the characters after the comma to poststr
L1=textscan(prestr,'%s');
r1=size(L1{1},1);   %the length of the string before comma (,)
if (poststr~='''')
    L2=textscan(poststr,'%s');
    r2=size(L2{1},1);   %the length of the string after comma(,)
else
    r2=0;    
end
p=str2num(L1{1}{1});%the index of the place
M0(1,k)=str2num(L1{1}{2});%M0 is the initial marking(row vector) k is the index of the place
k=k+1;
 for i=3:r1
      if isempty(strfind(L1{1}{i},':'))     %if the arc without the weight
          Post(p,str2num(L1{1}{i}))=1;
      else
          d=strfind(L1{1}{i},':');     %find the place of the colon (:)
          t=L1{1}{i}(1:d-1);    %assign the characters before the colon to t
          q=length(L1{1}{i});
          w=L1{1}{i}(d+1:q);    %assign the characters after the colon to w
          T=str2num(t);    %convert the characters before the colon to a numeric vector
          W=str2num(w);    %convert the characters after the colon to a numeric vector
          Post(p,T)=W;
      end
 end
 
 if r2~=0
    for i=1:r2
      if isempty(strfind(L2{1}{i},':'))     %if the arc without the weight
          Pre(p,str2num(L2{1}{i}))=1;
      else
          d=strfind(L2{1}{i},':');     %find the place of the colon (:)
          t2=L2{1}{i}(1:d-1);    %assign the characters before the colon to t
          q=length(L2{1}{i});
          w2=L2{1}{i}(d+1:q);    %assign the characters after the colon to w
          T=str2num(t2);    %convert the characters before the colon to a numeric vector
          W=str2num(w2);    %convert the characters after the colon to a numeric vector
          Pre(p,T)=W;
      end
      t=[];
      q=[];
    end
 end
  
  line=fgetl(fid);
  if ~isempty(strfind(line,'@'))%finish to read the intormation of the net structure
      break;
  end
end
fclose('all');
end





    









    



