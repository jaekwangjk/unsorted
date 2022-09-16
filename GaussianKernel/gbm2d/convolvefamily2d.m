function cval = convolvefamily2d(ind,val,dims)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% cval = convolvefamily2d(ind,val,dims)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Convolves the characteristic function of union of grains whose
% level set representation is stored in "ind" and "val" input vars
% with the kernel "KERNEL" (global variable, assumed present).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global KERNELalpha KERNELbeta Z; % These global variables *must* have been already
                      % allocated before calling this function! Both are
                      % of dimension "dims".

%% Convert level set data to volume of fluid representation:
[x,y] = ind2sub(dims,ind);
vf = loc_levset_to_volfluid(int32(x),int32(y),val,Z);

%% Carry out the convolution:
Kalphau = zeros(dims);
Kalphau(ind) = vf; % Characteristic function of the union of grains.
Kalphau = real(ifft2(fft2(Kalphau).*KERNELalpha));
Kbetau = zeros(dims);
Kbetau(ind) = vf; % Characteristic function of the union of grains.
Kbetau = real(ifft2(fft2(Kbetau).*KERNELbeta));
cval{1} = Kalphau; % Convolution values, returned.
cval{2} = Kbetau;