function psi_hat = GetPsi(q_hat,p)
persistent InvBT InvBC InvMat
if(p.model==1)
    if(isempty(InvBT))
        k = (2*pi/p.LX)*[0:p.N/2 -p.N/2+1:-1]';
        [KX,KY] = meshgrid(k,k);
        Laplacian = -(KX.^2+KY.^2);
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
else
    if(isempty(InvMat))
        k = (2*pi/p.LX)*[0:p.N/2 -p.N/2+1:-1]';
        [KX,KY] = meshgrid(k,k);
        K2 = KX.^2+KY.^2;
        K = sqrt(K2);
        InvMat = zeros([p.N p.N 2 2]);
        InvMat(:,:,1,1) = 1./(K.*tanh(K));
        InvMat(:,:,1,2) = -1./(K.*sinh(K));
        InvMat(:,:,2,1) = 1./(K.*sinh(K));
        InvMat(:,:,2,2) = -1./(K.*tanh(K));
        InvMat(1,1,:,:) = 0;
    end
    % Invert for psi
    psi_hat(:,:,1) = InvMat(:,:,1,1).*q_hat(:,:,1) + InvMat(:,:,1,2).*q_hat(:,:,2);
    psi_hat(:,:,2) = InvMat(:,:,2,1).*q_hat(:,:,1) + InvMat(:,:,2,2).*q_hat(:,:,2);
end
