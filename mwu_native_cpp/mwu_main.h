#ifndef _mwu_main_h_
#define _mwu_main_h_

#ifdef __cplusplus
extern "C" {
#endif

void run_mwu_cpp(double * Sigma,  
                 double * alpha,  
                 double * bsvm,   
                 int * posw,      
                 double * K,      
                 int * y,         
                 int n,           
                 int m,           
                 double eps,      
                 double ratio,    
                 double cutoff,   
                 double C,        
                 int norm1or2,    
                 int verbose      
                 );

/*
  void run_mwu_cpp(...);

  Output parameters:
  |==========|==========|=======================================|
  | Sigma    | double * | The kernel weights                    |       
  |----------|----------|---------------------------------------|
  | alpha    | double * | Support vector                        |       
  |----------|----------|---------------------------------------|
  | bsvm     | double * | Bias                                  |       
  |----------|----------|---------------------------------------|
  | posw     | int *    | Support indicators                    |       
  |----------|----------|---------------------------------------|

  Input parameters:
  |==========|==========|=======================================|
  | K        | double * | Kernel as one-dim array (columns,     |
  |          |          | then layers)                          |
  |----------|----------|---------------------------------------|
  | y        | int *    | Labels, +/-1                          |        
  |----------|----------|---------------------------------------|
  | n        | int      | Number of data points                 |        
  |----------|----------|---------------------------------------|
  | m        | int      | Number of kernels                     |        
  |----------|----------|---------------------------------------|
  | eps      | double   | Epsilon parameter                     |        
  |----------|----------|---------------------------------------|
  | ratio    | double   | Iteration multiplier                  |        
  |----------|----------|---------------------------------------|
  | cutoff   | double   | Exponentiation cutoff                 |        
  |----------|----------|---------------------------------------|
  | C        | double   | Margin parameter                      |        
  |----------|----------|---------------------------------------|
  | norm1or2 | int      | Is the soft margin 1-norm (1) or      |
  |          |          | 2-norm (2) or is it a hard margin (0) |                   
  |----------|----------|---------------------------------------|
  | verbose  | int      | Be noisy or not (boolean)             |        
  |----------|----------|---------------------------------------|
 */

enum kern_func_id {
  KERNEL_ID_LINEAR,    /* No parameters   */
  KERNEL_ID_POLY,      /* 1 param, degree */
  KERNEL_ID_GAUSSIAN,  /* 1 param, width  */
  KERNEL_ID_SIGMOID    /* No parameters   */
};

void run_mwu_dynamic(double * Sigma,                 
                     double * alpha,                 
                     double * bsvm,                  
                     int * posw,                     
                     enum kern_func_id * kern_funcs, 
                     double * kern_params,           
                     double * X,                     
                     int * y,                        
                     int n,                          
                     int m,                          
                     double eps,                     
                     double ratio,                   
                     double cutoff,                  
                     double C,                       
                     int norm1or2,                   
                     int verbose                     
                     );
/*
  void run_mwu_dynamic(...);
  
  Run MWU MKL algorithm on data X and labels y, and dynamically generate
  kernel entries as needed. The functions are specified as codes, and
  parameters are unlimited, provided that all the first parameters for the
  functions come first, then all the second parameters, etc. So if Pij is the
  jth parameter for the ith kernel function, the parameters would come in an
  array as

  P11, P21, ..., Pm1, P12, P22, ..., Pm2, P13, P23, ...

  The output parameters and remaining input parameters are the same as
  run_mwu_cpp above.
  
  Output parameters:
  |=============|=====================|=======================================|
  | Sigma       | double *            | The kernel weights                    |
  |-------------|---------------------|---------------------------------------|
  | alpha       | double *            | Support vector                        |
  |-------------|---------------------|---------------------------------------|
  | bsvm        | double *            | Bias                                  |
  |-------------|---------------------|---------------------------------------|
  | posw        | int *               | Support indicators                    |
  |-------------|---------------------|---------------------------------------|

  Input parameters:
  |=============|=====================|=======================================|
  | kern_funcs  | enum kern_func_id * | array of kernel function ids          |
  |-------------|---------------------|---------------------------------------|
  | kern_params | double *            | array of kernel parameters; all the   |
  |             |                     | first parameters for come first, then |
  |             |                     | the second, etc.                      |            
  |-------------|---------------------|---------------------------------------|
  | X           | double *            | Data as one-dim array (points as      |
  |             |                     | COLUMNS, and by columns first)        |         
  |-------------|---------------------|---------------------------------------|
  | y           | int *               | Labels, +/-1                          |
  |-------------|---------------------|---------------------------------------|
  | n           | int                 | Number of data points                 |
  |-------------|---------------------|---------------------------------------|
  | m           | int                 | Number of kernels                     |
  |-------------|---------------------|---------------------------------------|
  | eps         | double              | Epsilon parameter                     |
  |-------------|---------------------|---------------------------------------|
  | ratio       | double              | Iteration multiplier                  |
  |-------------|---------------------|---------------------------------------|
  | cutoff      | double              | Exponentiation cutoff                 |
  |-------------|---------------------|---------------------------------------|
  | C           | double              | Margin parameter                      |
  |-------------|---------------------|---------------------------------------|
  | norm1or2    | int                 | Is the soft margin 1-norm (1) or      |
  |             |                     | 2-norm (2) or is it a hard margin (0) |           
  |-------------|---------------------|---------------------------------------|
  | verbose     | int                 | Be noisy or not (boolean)             |
  |-------------|---------------------|---------------------------------------|
*/

#ifdef __cplusplus
}
#endif

#endif 
