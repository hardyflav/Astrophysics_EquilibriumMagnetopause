
clear all
clc

DeltaX = 1 ;
DeltaY = 1 ;
DeltaZ = 1 ;

Nx = 10;
Ny = 10;
Nz = 10;

% Initialisation
U = ones( Nx, Ny, Nz );
Un = U;

% Boundary Conditions
U(end,:,:) =  U(end-1,:,:) + DeltaX .* BC_XRight;


U_Int = U_Large(2:end-1, 2:end-1, 2:end-1);

U_xp1 = U_Large(3:end, 2:end-1, 2:end-1);
U_xm1 = U_Large(1:end-2, 2:end-1, 2:end-1);

U_yp1 = U_Large(2:end-1, 3:end, 2:end-1);
U_ym1 = U_Large(2:end-1, 1:end-2, 2:end-1);

U_zp1 = U_Large(2:end-1, 2:end-1, 3:end);
U_zm1 = U_Large(2:end-1, 2:end-1, 1:end-2);

% Partial Derivatives
dU_x2 = (U_xp1 - 2*U_Int - U_xm1) ./ (DeltaX^2);
dU_y2 = (U_yp1 - 2*U_Int - U_ym1) ./ (DeltaY^2);
dU_z2 = (U_zp1 - 2*U_Int - U_zm1) ./ (DeltaZ^2);

% Boundary Conditions
( U_Int(end,:,:) - U_xm1(end,:,:) ) ./ (DeltaX) = BC_XRight;
( U_Int(end,:,:) - U_xm1(end,:,:) ) ./ (DeltaX) = BC_XLeft;

( U_Int(:,end,:) - U_ym1(:,end,:) ) ./ (DeltaY)  = BC_YTop;
( U_Int(:,2,:) - U_Int(:,1,:) ) ./ (DeltaY)  = BC_YBottom;







dU_y_Top = (U_Int - U_xm1) ./ (DeltaX);



B = {'DN' 'DN' 'P'}
Laplacian()

[~,~,A] = laplacian([10,10,10],{'DN' 'NN' 'NN'});
spy(A)

A*U