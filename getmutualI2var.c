#include "mex.h"
#include "math.h"
#include "stdlib.h"
#include "stdio.h"
#define NULL 0
#define LEN sizeof(struct clumb)
struct clumb
{int subc;
 struct clumb *next;
};
struct clumb *insert(struct clumb *head,struct clumb *stud)
{
 struct clumb *sub_node0,*sub_node1,*sub_node2;  
 sub_node1=head;
 sub_node0=stud;
 sub_node2=NULL;
 if(head==NULL){
     head=sub_node0;
     sub_node0->next=NULL;
 }else{
     while((sub_node0->subc>sub_node1->subc)&&(sub_node1->next!=NULL)){
         sub_node2=sub_node1;
         sub_node1=sub_node1->next;         
     }
     if(sub_node0->subc<=sub_node1->subc){
         if(head==sub_node1) head=sub_node0;
         else sub_node2->next=sub_node0;
         sub_node0->next=sub_node1;   
     }else{
         sub_node1->next=sub_node0;
         sub_node0->next=NULL;    
     }
 }
 return(head);
}
struct clumb *del(struct clumb *head,int num)
{
 struct clumb *sub_node1,*sub_node2;
 sub_node1=head;
//  sub_node1=sub_node1->next;
 if(head==NULL) mexPrintf("list is null\n");
 else
 {while(num!=sub_node1->subc && sub_node1->next!=NULL){
     sub_node2=sub_node1;
     sub_node1=sub_node1->next;  
 }
 if(num==sub_node1->subc){
     if(sub_node1==head) head=sub_node1->next;
     else sub_node2->next=sub_node1->next;
     mxFree(sub_node1);
 }
 else mexPrintf("not been found\n"); 
}
 return(head);
}
void print(struct clumb *head)
{
 struct clumb *node;
 node=head;
 if(head!=NULL)
     do
     {
     mexPrintf("%d\n",node->subc);
     node=node->next;
     } while(node!=NULL);
}

void release(struct clumb *head)
{
 struct clumb *node;
//  node=head;
 if(head!=NULL)
     do
     {
     node=head;
     head=head->next;
     mxFree(node);
     } while(head!=NULL);
}

void listtomatric(struct clumb *head,int *submat,int len_mat)
{
 struct clumb *node;
 int subcount=0;
 node=head;
 if(head!=NULL)
     do
     {
     if(subcount!=0) submat[subcount-1]=node->subc;
     subcount=subcount+1;
     node=node->next;
     } while(node!=NULL);
}
double myentropy(int mytable[],int tab_len,int num_sample)
{
    int i;
    double entr=0.0;
    for(i=0;i<tab_len;i++){
        if(mytable[i]!=0){
            entr=entr+(((double)mytable[i])/((double)num_sample))*(log(((double)mytable[i])/((double)num_sample))/log(2)); 
        }
    }
    entr=entr*(-1);
    return entr;
}
double mutual_I(int vector_x[],struct clumb *head,int certain,int num_sample,int len_seg)
{
    int i,j,subk=0,*subcertain,*temseg,*array_seg;
    int *array_seg_cer,*array_cer;
    double myI,H1,H2,H3;
    subcertain=(int *) mxCalloc(certain,sizeof(int));
    array_seg=(int *) mxCalloc(len_seg,sizeof(int));
    array_seg_cer=(int *) mxCalloc(certain*len_seg,sizeof(int));
    array_cer=(int *) mxCalloc(certain,sizeof(int));
    for(i=0;i<certain*len_seg;i++){
        if(i<certain){
            array_cer[i]=0;
            subcertain[i]=i+1;
        }
        if(i<len_seg) array_seg[i]=0;
        array_seg_cer[i]=0;
    }
    temseg=(int *) mxCalloc(len_seg,sizeof(int));
    listtomatric(head,temseg,len_seg);
    for(i=0;i<num_sample;i++){
        for(j=0;j<certain;j++){
            if(vector_x[i]==subcertain[j]){
                array_cer[j]=array_cer[j]+1;
                if(i<temseg[subk]){
                    array_seg_cer[subk+j*len_seg]=array_seg_cer[subk+j*len_seg]+1;
                    array_seg[subk]=array_seg[subk]+1;
                }else{
                    subk=subk+1;
                    array_seg[subk]=array_seg[subk]+1;
                    array_seg_cer[subk+j*len_seg]=array_seg_cer[subk+j*len_seg]+1;
                }
              break;
            }
        }
    }
    H1=myentropy(array_seg_cer,certain*len_seg,num_sample); 
    H2=myentropy(array_seg,len_seg,num_sample);
    H3=myentropy(array_cer,certain,num_sample);
//     if (H2<H3){
//        myI=(H2+H3-H1)/H2;
//     }else{
//        myI=(H2+H3-H1)/H3;
//     }
//     myI=2*(H2+H3-H1)/(H3+H2);
    if (len_seg<certain){
       myI=(H2+H3-H1)/(log(len_seg)/log(2));
    }else{
        myI=(H2+H3-H1)/(log(certain)/log(2));
    }
    return myI;
    mxFree(temseg);
    mxFree(subcertain);
    mxFree(array_cer);
    mxFree(array_seg);
    mxFree(array_seg_cer);
}
void mexFunction(int nlhs,mxArray *plhs[],int nrhs, const mxArray *prhs[])
{
    int *vector_x,*best_c,*out2;
    int num_sample,segm,certain,len_bestc,len_segm,i,j,count,num;
    double *myI,I0,tem_I;
    struct clumb *head_mybestc,*head_tembestc,*node1,*node2,*node_tem;
    vector_x=mxGetData(prhs[0]);
    best_c=mxGetData(prhs[1]);
    segm=*((int *)mxGetData(prhs[2]));
    certain=*((int *)mxGetData(prhs[3]));
    num_sample=*((int *)mxGetData(prhs[4]));
    len_bestc=*((int *)mxGetData(prhs[5]));
    if(len_bestc>segm-1) len_segm=segm-1;
    else len_segm=len_bestc;
    plhs[0]=mxCreateDoubleMatrix(len_segm,1,mxREAL);
    myI=mxGetPr(plhs[0]);
    plhs[1]=mxCreateNumericMatrix(len_segm+2,1,mxINT32_CLASS,mxREAL);
    out2=mxGetData(plhs[1]);
    head_mybestc=NULL;
    head_tembestc=NULL;
    node1=node2=(struct clumb *) mxMalloc(LEN);
    node1->subc=best_c[0];
    head_mybestc=node1;
    node2=node1;
    node1=(struct clumb *) mxMalloc(LEN);
    node1->subc=best_c[len_bestc+1];
    node2->next=node1;
    node1->next=NULL;  
    node1=node2=(struct clumb *) mxMalloc(LEN);
    for(i=1;i<len_bestc+1;i++){
        node1->subc=best_c[i];
        if(i==1) head_tembestc=node1;
        else node2->next=node1;
        node2=node1;
        node1=(struct clumb *) mxMalloc(LEN);
    }
    node2->next=NULL;    
    count=0;
    out2[0]=0;
    for(i=0;i<len_segm;i++){
        I0=0.0;
        node_tem=head_tembestc;
        if(head_tembestc!=NULL)
        do
        {
        j=node_tem->subc;
        node1=(struct clumb *) mxMalloc(LEN);
        node1->subc=j;
        head_mybestc=insert(head_mybestc,node1);
        tem_I=mutual_I(vector_x,head_mybestc,certain,num_sample,count+2);
        if (fabs(tem_I)>fabs(I0)){
            I0=tem_I;
            num=j;
        }
        head_mybestc=del(head_mybestc,j);
        node_tem=node_tem->next;
        }while(node_tem!=NULL);
        myI[count]=I0;
        count++; 
        node1=(struct clumb *) mxMalloc(LEN);
        node1->subc=num;
        out2[i+1]=num;
        head_mybestc=insert(head_mybestc,node1);
        head_tembestc=del(head_tembestc,num);
    }
    out2[len_segm+1]=num_sample;
    release(head_mybestc);
    release(head_tembestc);
    return;
}