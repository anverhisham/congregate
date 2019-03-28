 ### BRIEF DESCRIPTION
       Matlab code to fetch values from an array of cells (or an array of structures) and create a matrix.

 ### DETAILED DESCRIPTION
       Have you ever tried to create a submatrix from an array of cells? Or extract values of a specific member from an array of strctures? 
       If so, this function is for you!! Please see the examples below to get more info upon this function.

 ### EXAMPLES

       1. Extracting submatrix from an array of cells
           A = {11,12,13;21,22,23};
           congregate('A{1:2,2:3}')
               ans =
                   12    13
                   22    23

       2. Converting an array of structures to a 2x3 matrix
           S(1).a=[1,2,3]';
           S(2).a=[4,5,6]';
           congregate('S(:).a')
                 ans =
                      1     2     3
                      4     5     6

       3. Converting an array of structures of structures to a 2x2x3 matrix
           T(1).a(1).b=[111,112,113]';
           T(1).a(2).b=[121,122,123]';
           T(2).a(1).b=[211,212,213]';
           T(2).a(2).b=[221,222,223]';
           congregate('T(:).a(:).b')
                 ans(:,:,1) =
                    111   121
                    211   221
                 ans(:,:,2) =
                    112   122
                    212   222
                 ans(:,:,3) =
                    113   123
                    213   223


 ### SYNTAX AND OPTIONS
 #### Syntax 
       congregate('expression',isFillWithNan,squeezeLevel)
 #### Description
 #####      isFillWithNan
If set to true, then the absent fields would be regarded as NaN. 
Otherwise, if the field is absent then throws an error. 
Default value is true.       
       
#####       squeezeLevel 
Determines the level of squeezing of the input. Squeezing is defined as the removal of all singleton dimensions in a matrix/cell.
               'squeezeLevel' = 0 is the default, which implies no squeezing
               'squeezeLevel' = 1 implies moderate squeezing, where the content inside cells are squeezed.
               'squeezeLevel' = 2 implies high level squeezing, where the content inside both cells and matrices are squeezed.

 ### DETAILED NOTES
       1. The number of dimensions equals to the number of brackets appearing in 'inString' and the dimension of last field.
       2. Single ':' in the expression implies output is a column vector of cells. 
     

 
