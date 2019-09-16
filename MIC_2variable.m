function [MIC,bestc,mycertain,certain_seg,mutual_I_1,mutual_I_2]=MIC_2variable(data,B,c,n)
% n is the numbers of sample
% mutual_I_2 represent fixed the first column of data(x1)
% mutual_I_1 represent fixed the second column of data(x2)
% B=n^0.6
% c is max clumps (c*x)


mutual_I_2=zeros(round(B/2)-1,round(B/2)-1);
mutual_I_1=zeros(round(B/2)-1,round(B/2)-1);
last_c2=zeros(round(B/2)+1,round(B/2)+1);
last_c1=zeros(round(B/2)+1,round(B/2)+1);
n=int32(n);
% randnum=randperm(n)';
% data=data(randnum,:);
% [~,pos1]=sortrows(data,1);
% [~,pos2]=sortrows(data,2);
[value1,pos1]=sortrows(data,1);
[value2,pos2]=sortrows(data,2);
D1=nan(n,2);
D1=int32(D1);
for i=2:round(B/2)
    Q_x1=equipartitionYaxis2c(value1(:,1),int32(i),n);
    Q_x2=equipartitionYaxis2c(value2(:,2),int32(i),n);
%     Q_x=equipartitionYaxis(int32(i),n);
    sub_max_seg=int32(round(B/i));
    avg=int32(n/(round(B/i)*c));
    if max(Q_x1)==i
        D1(pos1,1)=Q_x1;
        vector1=D1(pos2,1);
        c_x2=getsuper2var(vector1,avg,n);
        len2=int32(length(c_x2)-2);
        [temI2,tem_c2]=getmutualI2var(vector1,c_x2,sub_max_seg,int32(i),n,len2);
        mutual_I_2(i-1,1:length(temI2))=temI2';%./log2(min_seg);
        last_c2(i-1,1:length(tem_c2))=tem_c2;
%     min_seg=2:length(temI2)+1;
%     min_seg=min([min_seg;zeros(1,length(temI2))+i]);
    end
    if max(Q_x2)==i
        D1(pos2,2)=Q_x2;
        vector2=D1(pos1,2);
        c_x1=getsuper2var(vector2,avg,n);
        len1=int32(length(c_x1)-2);
        [temI1,tem_c1]=getmutualI2var(vector2,c_x1,sub_max_seg,int32(i),n,len1);
    %     min_seg=2:length(temI1)+1;
    %     min_seg=min([min_seg;zeros(1,length(temI1))+i]);
        mutual_I_1(1:length(temI1),i-1)=temI1;%./log2(min_seg');
        last_c1(1:length(tem_c1),i-1)=tem_c1';
    end
end
MIC=max(max([mutual_I_1,mutual_I_2]));
position=find(mutual_I_2==MIC, 1);
if ~isempty(position)
    pos_row=mod(position(1)-1,round(B/2)-1)+1;
    pos_col=floor((position(1)-1)/(round(B/2)-1))+1;
    mycertain=1;
    certain_seg=pos_row+1;
    bestc=last_c2(pos_row,2:pos_col+1);
    position2=find(mutual_I_1==MIC, 1);
    if ~isempty(position2)
        pos_row2=mod(position2(1)-1,round(B/2)-1)+1;
        pos_col2=floor((position2(1)-1)/(round(B/2)-1))+1;
        mycertain2=2;
        certain_seg2=pos_col2+1;
        bestc2=last_c1(2:pos_row2+1,pos_col2)';
        if certain_seg2*(length(bestc2)+1)<certain_seg*(length(bestc)+1)
            mycertain=mycertain2;
            certain_seg=certain_seg2;
            bestc=bestc2;
        end
    end
else
    position=find(mutual_I_1==MIC, 1);
    pos_row=mod(position(1)-1,round(B/2)-1)+1;
    pos_col=floor((position(1)-1)/(round(B/2)-1))+1;
    mycertain=2;
    certain_seg=pos_col+1;
    bestc=last_c1(2:pos_row+1,pos_col)';
end
end
