 -BRIEF DESCRIPTION:
       Fetch values from an array of cells (or an array of structures) and create a matrix.

 -DETAILED DESCRIPTION:
       Have you ever tried to create a submatrix from an array of cells? Or extract values of a specific member from an array of strctures? 
       If so, this function is for you!! Please see the examples below to get more info upon this function.

 -EXAMPLES:

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


 -OPTIONS:
       1. 'isFillWithNan' = true (Default) means, absent fields would be regarded as NaN 
       2. 'squeezeLevel': determines the level of squeezing of the input. Squeezing is defined as the removal of all singleton dimensions in a matrix/cell.
               'squeezeLevel' = 0 is the default, which implies no squeezing
               'squeezeLevel' = 1 implies moderate squeezing, where the content inside cells are squeezed.
               'squeezeLevel' = 2 implies high level squeezing, where the content inside both cells and matrices are squeezed.

 -DETAILED NOTES:
       1. Here number of dimensions = number of brackets appearing in 'inString' + dimension of last field  
       2. Single ':' in structExprStr implies output is a column vector of cells  
       3. s{[1,2;3,4],[1,2]}.a will result in 1D output

 -PROGRAMMING NOTES:
       1. There are some design errors in MATLAB software, as follows
           a). MATLAB ignores, tralining singleton dimensions. 
           b). 'a(2:5)' and 'a{:}' are outputted as row-vectors, instead of column vecotors
           c). ndims of column vector is 2 (instead of 1)!!. 
       Solution:
           The problem (a) is rectified using 'nDimPerScope'.  To overcome (c), I implemented ndims_().

 -TROUBLESHOOTS:
       1. One way to check if this function works or not, is to execute your input by replacing all variable indices with '1's 
           i.e.  s{1}.s1(1).s2{1}.a

 -TODO:
       1. Include accepting 'end' keyword.
       2. Nonworking Input: congregate_('in{:,2,3}{1}.F')% 

 -VERSIONS:
       Version 1 release: 2019 January 31

 -AUTHORS:
       Anver Hisham <anverhisham@gmail.com>

 -LICENSE:
       This work is under a Creative Commons Share Alike (CC-SA) license.
