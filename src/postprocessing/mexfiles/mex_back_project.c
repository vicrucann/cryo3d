#include <string.h>
#include <math.h>
#include "mex.h"

void mexFunction(int nlhs, mxArray *plhs[ ],int nrhs, const mxArray *prhs[ ]) 
{
    double *vol,*img;
    double *n,*n_z,*n_x,*n_y;
    double mask_r;
    double img_c_x,img_c_y;
    double vol_c_x,vol_c_y,vol_c_z;
    double img_r,z_lim;
    double z;
    double v_x,v_y,v_z;
    int i_v_x,i_v_y,i_v_z;
    int i1,j1,k1;
    double wt;
    double a_v_x,a_v_y,a_v_z;
    int v_size;
    int i,j,axis,n_axes;
    mwSignedIndex dims[3];
    
    
    /**********************************************/
       /*Retrieve the input data */
       /* The volume */
       img=mxGetPr(prhs[0]);
       v_size=(int)*(mxGetPr(prhs[1]));;
              /* The coordinate system */
       n=mxGetPr(prhs[2]);
       n_axes=(int)*(mxGetPr(prhs[3]));
              /* The radius */
       mask_r=*(mxGetPr(prhs[4]));
       

       
   /************************************************/
       
       /* Create and initialize the output volume */
       dims[0] = v_size;
       dims[1] = v_size;
       dims[2]= v_size;
       plhs[0] = mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL);
       vol=mxGetPr(plhs[0]);
       for (i=0;i<v_size*v_size*v_size;i++) *(vol+i)=0;
        
   /*************************************************/
        /* Create constants */
        
        img_c_x=((double)v_size)/2-0.5;
        img_c_y=((double)v_size)/2-0.5;
        
        vol_c_x=((double)v_size)/2-0.5;
        vol_c_y=((double)v_size)/2-0.5;
        vol_c_z=((double)v_size)/2-0.5;
        
    /**************************************************/
   for (axis=0;axis<n_axes;axis++)
   {    /* For each axis */
       n_x=n +9*axis;
       n_y=n+9*axis+3;
       n_z=n+9*axis+6;
       
       
        /* Back_project onto the image */
        for (i=0;i<v_size;i++)
            for (j=0;j<v_size;j++)
        {
            /* Calculate the radius in the image */
            img_r=sqrt((i-img_c_x)*(i-img_c_x)+(j-img_c_y)*(j-img_c_y));
            /* If the radius is less than mask_r */
            if (img_r < mask_r)
            { /* Radius is small, get the z limits */
                z_lim=sqrt(mask_r*mask_r-img_r*img_r);
         
                
                for (z=-z_lim;z<=z_lim;z++)
                {
                    /* Calculate an index into the volume */
                    /* Double indices */
               /*   v_x=((i-img_c_x)+vol_c_x);
                    v_y=((j-img_c_y)+vol_c_y);
                    v_z=(z+vol_c_z);
                */
                    v_x=(i-img_c_x)*(*n_x)+(j-img_c_y)*(*n_y)+z*(*n_z)+vol_c_x;
                    v_y=(i-img_c_x)*(*(n_x+1))+(j-img_c_y)*(*(n_y+1))+z*(*(n_z+1))+vol_c_y;
                    v_z=(i-img_c_x)*(*(n_x+2))+(j-img_c_y)*(*(n_y+2))+z*(*(n_z+2))+vol_c_z;
                    
                    /*Integer indices and fractional offsets */
                    i_v_x=(int)floor(v_x);
                    i_v_y=(int)floor(v_y);
                    i_v_z=(int)floor(v_z);
                    
                    a_v_x=v_x-(double)i_v_x;
                    a_v_y=v_y-(double)i_v_y;
                    a_v_z=v_z-(double)i_v_z;

                    /* Multi-linear interpolation */
                    for (i1=0;i1<=1;i1++)
                        for (j1=0;j1<=1;j1++)
                            for (k1=0;k1<=1;k1++) 
                            {
                                wt=((1-a_v_x)*(1-i1)+a_v_x*i1)
                                     *((1-a_v_y)*(1-j1)+a_v_y*j1)
                                     *((1-a_v_z)*(1-k1)+a_v_z*k1);
                        *(vol+(i_v_x+i1)+(i_v_y+j1)*v_size+(i_v_z+k1)*v_size*v_size) += 
                              (*(img+axis*v_size*v_size+i+j*v_size))*wt;
                            } /* End k1 loop */
                    
                } /* End z loop */
         
                
            } /* End img_r < mask_r condition */
        } /* End j loop */
   } /* End axis */  
}