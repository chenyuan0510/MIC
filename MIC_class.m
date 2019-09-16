function [MIC,mybestc]=MIC_class(data,B,c)
% n is the numbers of sample
% x*y<=n^0.6 ----> n^B
% c is max clumps (c*x)

n=length(data(:,1));
if size(data,2)<3
     temcount=tabulate(data(:,1));
%     log_v=-sum((temcount(:,2)/n).*mylog(temcount(:,2)/n));
    certain=length(temcount(:,1));
    max_seg=round(n^B/certain);
%     randnum=randperm(n)';
%     data=data(randnum,:);
    [~,I2]=sortrows(data,2);
    vector1=data(I2,1);
    vector1=int32(vector1);
    n=int32(n);
    avg=int32(n/(max_seg*c));
    best_c=getsuper2var(vector1,avg,n);
    len2=int32(length(best_c)-2);
    [mutual_I_2,tem_c]=getmutualI2var(vector1,best_c,int32(max_seg),int32(certain),n,len2);
%     [mutual_I_2,tem_c]=getmutualIclass(vector1,best_c,int32(max_seg),int32(certain),n,len2);
%     min_seg=2:length(mutual_I_2)+1;
%     min_seg=min([min_seg;zeros(1,length(mutual_I_2))+certain]);
%     mutual_I_2=mutual_I_2';%./log2(min_seg);
    [MIC,pos]=max(mutual_I_2);
    mybestc=tem_c(2:pos+1)';
end
end