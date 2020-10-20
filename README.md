# MIC
Maximal information coefficient (MIC)

**1） MATLAB is the tool of NDC;**

**2) Make sure that c++ has installed in your computer; In the matlab runtime environment, use the ```mex -setup``` command to set: <use of'Microsoft Visual C++' for C++ language compilation>**

**3) Running program “make.m” to compile the c files to mex files;**

> the first column of data is Y (dependent variable), the rest of the columns (X) independent variable;

    num=randperm(size(data,1));

    data=data(num',:);# scramble the samples

> while Y is numerical data

    sample_num=size(data,1);

    [MIC,bestc,mycertain,certain_seg,mutual_I_1,mutual_I_2]=MIC_2variable(data,B,c,n)

    [MIC,bestc,mycertain,certain_seg,mutual_I_1,mutual_I_2]=MIC_2variable(data,sample_num^0.55,5,sample_num);

> while Y is discrete data

    [MIC,~]=MIC_class(data,0.55,5)
    
## Contact me: chenyuan0510@126.com
