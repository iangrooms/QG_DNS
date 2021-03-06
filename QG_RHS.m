function RHS = QG_RHS(q_hat,p)
% Function takes Fourier coefficients of PV (q_hat) and struct containing
% parameters (p) and evaluates RHS of 2-layer QG equations except for
% high-k dissipation. Returns Fourier coefficients of RHS.
% Jacobian is dealiased via 3/2 rule.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent dX dY DX DY Laplacian InvBT InvBC
if isempty(DX)
    k = (2*pi/p.LX)*[0:p.N/2 -p.N/2+1:-1]';
    dX = 1i*repmat(k',[p.N 1]);
    dY = 1i*repmat(k,[1 p.N]);
    Laplacian = dX.^2+dY.^2;
    k = (2*pi/p.LX)*[0:p.N/2-1 0 -p.N/2+1:-1]';
    dX = 1i*repmat(k',[p.N 1]);
    dY = 1i*repmat(k,[1 p.N]);
    % For the dealiased jacobian:
    k = (2*pi/p.LX)*[0:.75*p.N-1 0 -.75*p.N+1:-1]';
    DX = 1i*repmat(k',[1.5*p.N 1 2]);
    DY = 1i*repmat(k,[1 1.5*p.N 2]);
    clear k
end

% Invert for psi
psi_hat = GetPsi(q_hat,p);

RHS = zeros([p.N p.N 2]);
if( p.model==1 )
    % beta plus mean shear
    RHS(:,:,1) =-dX.*q_hat(:,:,1)-(p.kb^2 + p.kd^2)*dX.*psi_hat(:,:,1);
    RHS(:,:,2) = dX.*q_hat(:,:,2)-(p.kb^2 - p.kd^2)*dX.*psi_hat(:,:,2);
else
    % mean buoyancy gradient plus mean shear
    RHS(:,:,1) =-dX.*q_hat(:,:,1)+2*dX.*psi_hat(:,:,1);
    RHS(:,:,2) = dX.*q_hat(:,:,2)+2*dX.*psi_hat(:,:,2);
end

% Bottom drag
u =-real(ifft2(dY.*psi_hat(:,:,2)));
v = real(ifft2(dX.*psi_hat(:,:,2)));
U = sqrt(u.^2+v.^2);
RHS(:,:,2) = RHS(:,:,2) - p.cd*dX.*fft2(U.*v)+p.cd*dY.*fft2(U.*u);
clear u v U
RHS(:,:,2) = RHS(:,:,2) - p.r*Laplacian.*psi_hat(:,:,2);

% The following code computes the Jacobian J[psi,q], dealiased using a 3/2-rule.
% I.e., if you set N=256 then the code uses 3*N/2=384 Fourier modes (in each
% direction) to compute the Jacobian.
    % physical space, 3/2 grid; factor of (9/4) scales fft
    Psi_hat = zeros([1.5*p.N 1.5*p.N 2]);
    Psi_hat(1:p.N/2+1,1:p.N/2+1,:) = (9/4)*psi_hat(1:p.N/2+1,1:p.N/2+1,:);
    Psi_hat(1:p.N/2+1,p.N+2:1.5*p.N,:) = (9/4)*psi_hat(1:p.N/2+1,p.N/2+2:p.N,:);
    Psi_hat(p.N+2:1.5*p.N,1:p.N/2+1,:) = (9/4)*psi_hat(p.N/2+2:p.N,1:p.N/2+1,:);
    Psi_hat(p.N+2:1.5*p.N,p.N+2:1.5*p.N,:) = (9/4)*psi_hat(p.N/2+2:p.N,p.N/2+2:p.N,:);
    Q_hat = zeros([1.5*p.N 1.5*p.N 2]);
    Q_hat(1:p.N/2+1,1:p.N/2+1,:) = (9/4)*q_hat(1:p.N/2+1,1:p.N/2+1,:);
    Q_hat(1:p.N/2+1,p.N+2:1.5*p.N,:) = (9/4)*q_hat(1:p.N/2+1,p.N/2+2:p.N,:);
    Q_hat(p.N+2:1.5*p.N,1:p.N/2+1,:) = (9/4)*q_hat(p.N/2+2:p.N,1:p.N/2+1,:);
    Q_hat(p.N+2:1.5*p.N,p.N+2:1.5*p.N,:) = (9/4)*q_hat(p.N/2+2:p.N,p.N/2+2:p.N,:);
    % calculate u.gradq on 3/2 grid
    u = real(ifft2(-DY.*Psi_hat));
    v = real(ifft2( DX.*Psi_hat));
    qx= real(ifft2( DX.*Q_hat));
    qy= real(ifft2( DY.*Q_hat));
    jaco_real = u.*qx+v.*qy;
    % fft, 3/2 grid; factor of (4/9) scales fft
    Jaco_hat = (4/9)*fft2(jaco_real);
    % reduce to normal grid
    jaco_hat = zeros([p.N p.N 2]);
    jaco_hat(1:p.N/2+1,1:p.N/2+1,:) = Jaco_hat(1:p.N/2+1,1:p.N/2+1,:);
    jaco_hat(1:p.N/2+1,p.N/2+2:p.N,:) = Jaco_hat(1:p.N/2+1,p.N+2:1.5*p.N,:);
    jaco_hat(p.N/2+2:p.N,1:p.N/2+1,:) = Jaco_hat(p.N+2:1.5*p.N,1:p.N/2+1,:);
    jaco_hat(p.N/2+2:p.N,p.N/2+2:p.N,:) = Jaco_hat(p.N+2:1.5*p.N,p.N+2:1.5*p.N,:);
    
% Put it all together
RHS = RHS - jaco_hat;
