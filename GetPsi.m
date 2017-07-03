function psi_hat = GetPsi(q_hat,p)
persistent InvBT InvBC
if(isempty(InvBT))
    k = (2*pi/p.LX)*[0:p.N/2 -p.N/2+1:-1]';
    dX = 1i*repmat(k',[p.N 1]);
    dY = 1i*repmat(k,[1 p.N]);
    Laplacian = dX.^2+dY.^2;
    InvBT = 1./Laplacian; InvBT(1,1) = 0;
    InvBC = 1./(Laplacian-p.kd^2);InvBC(1,1) = 0;
end
% Invert for psi
q_bt = .5*(q_hat(:,:,1) + q_hat(:,:,2));
q_bc = .5*(q_hat(:,:,1) - q_hat(:,:,2));
psi_bt = InvBT.*q_bt;
psi_bc = InvBC.*q_bc;
psi_hat(:,:,2) = psi_bt-psi_bc;
psi_hat(:,:,1) = psi_bt+psi_bc;
%psi = real(ifft2(psi_hat));
