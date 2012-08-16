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

#ifdef __cplusplus
}
#endif

#endif 
