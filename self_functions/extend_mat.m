function [mat_out] = extend_mat(mat,dim)
%   Detailed explanation goes here
    sz = size(mat);
    sz(dim) = sz(dim)+2;
    mat_out = zeros(sz);
    if(dim==1)
        mat_out(2:end-1,:,:) = mat;
        mat_out(1,:,:) = mat_out(2,:,:);
        mat_out(end,:,:) = mat_out(end-1,:,:);
    end
end

