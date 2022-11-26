function cval = Module_TECH_convolvefamily2d(ind,val,dims)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% cval = convolvefamily2doriginal(ind,val,dims)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Convolves the characteristic function of union of grains whose
% level set representation is stored in "ind" and "val" input vars
% with the kernel "KERNEL" (global variable, assumed present).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global KERNEL Z; % These global variables *must* have been already
                 % allocated before calling this function! Both are
                 % of dimension "dims".

%% Convert level set data to volume of fluid representation:
[x,y] = ind2sub(dims,ind);

% This is C function
vf = loc_levset_to_volfluid(int32(x),int32(y),val,Z);

%% Carry out the convolution:
Ku = zeros(dims);
Ku(ind) = vf; % Characteristic function of the union of grains.
Ku = real(ifft2( fft2(Ku) .* KERNEL ));
cval = Ku; % Convolution values, returned.