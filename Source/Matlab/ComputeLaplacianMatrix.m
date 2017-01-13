function [LapMatrix,Affmatrix] = ComputeLaplacianMatrix(MatrixSize,row_inds,col_inds,AffValue12)
% syntax:
%       [LapMatrix,Affmatrix] = ComputeLaplacianMatrix(data,row_inds,col_inds,AffValue12)
% description:
%       Transform neighbor information and affinity value to Affinity matrix and Compute Laplacian matrix
%
% Input arguments:
%       MatrixSize : Size of Affinity matrix and Lapalcian matrix 
%       row_inds : row index of affinity matrix    
%       col_inds : col index of affinity matrix
%       AffValue12 : affinity value of affinity matrix
% Output arguments:
%       LapMatrix :  Lapalacian matrix with (MatrixSize)x(MatrixSize)
%       sparse matrix
%       Affmatrix :  affinity matrix with (MatrixSize)x(MatrixSize) sparse
%       matrix
%  Jia-Hao Syu <jiahaushiu4@gmail.com>
%  December 2016
% 
% ex. [row_inds,col_inds,AffValue12]                            affinity matrix
%     [   1         2       0.1    ]                     [   0        0.1       0   ]  
%     [   2         1       0.1    ]        ====>        [  0.1        0       0.4  ]
%     [   2         3       0.4    ]                     [   0        0.4      0    ] 
%     [   3         2       0.4    ]
%     ==================================================================================
%           affinity matrix                         Laplacian matrix
%    [   0        0.1       0   ]               [  0.1   -0.1     0   ]
%    [  0.1        0       0.4  ]     =====>    [ -0.1    0.5   -0.4  ]
%    [   0        0.4       0  ]                [   0    -0.4    0.4  ]

[rowI,colI] = size(row_inds);
maxL = max(row_inds);

if maxL>5000 % sparse is faster in larger matrix 
   
        Affmatrix = sparse(row_inds,col_inds,AffValue12,MatrixSize,MatrixSize);
    
else         % C++ code is faster in dense matrix
    if rowI>colI
        Affmatrix  = IndexToMatrix(row_inds,col_inds,AffValue12,size(row_inds,1),MatrixSize);
        Affmatrix = sparse(Affmatrix);
    else
        Affmatrix  = IndexToMatrix(row_inds,col_inds,AffValue12,size(row_inds,2),MatrixSize);
        Affmatrix = sparse(Affmatrix);
    end

end

% Transform Affinity matrix to Symmetric affinity matrix
Affmatrix = (Affmatrix'+Affmatrix)/2;

DegreeVector = sum(Affmatrix,2); % degree vector

% Laplacia matrix = diagonal degree matrix - affinity matrix
LapMatrix = spdiags(DegreeVector,0,MatrixSize,MatrixSize)-Affmatrix;
