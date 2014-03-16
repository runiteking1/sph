function check_derivation

%@T
% \section{Sanity checks}
%
% I fairly regularly make typographical and copying errors when I do
% algebra and implement it in code.  In order to stay sane when I
% actually write something somewhat complicated, I find it helpful to
% put together little test scripts to check my work numerically.
% For your edification, in this section I give my MATLAB test script 
% corresponding to the derivation in these notes.  The test script 
% is done in MATLAB.
%
% I begin by implementing the functions $f(q)$, the normalizing constants,
% and the kernel functions for each of the three kernels.
%@c
fp6 = @(q) (1-q.^2).^3;
fsp = @(q) (1-q).^3;
fvi = @(q) q.^2 - 0.5*q.^3 + 0.5./q - 1;
Cp6 = 64*pi/315;
Csp = pi/15;
Cvi = 2*pi/15;
Wp6 = @(r,h) 1/Cp6/h^3 * fp6( norm(r)/h );
Wsp = @(r,h) 1/Csp/h^3 * fsp( norm(r)/h );
Wvi = @(r,h) 1/Cvi/h^3 * fvi( norm(r)/h );

%@T
%
% I computed the normalization constants analytically, but I'm
% prone to algebra mistakes when I compute integrals by hand.
% Let's check against MATLAB's [[quad]] function.
%@c
fprintf('Relerr for normalization constants:\n');
nerr_p6 = quad( @(q) 4*pi*q.^2.*fp6(q)/Cp6, 0,     1 ) - 1;
nerr_sp = quad( @(q) 4*pi*q.^2.*fsp(q)/Csp, 0,     1 ) - 1;
nerr_vi = quad( @(q) 4*pi*q.^2.*fvi(q)/Cvi, 1e-12, 1 ) - 1;
fprintf('  Cp6: %g\n', nerr_p6);
fprintf('  Csp: %g\n', nerr_sp);
fprintf('  Cvi: %g\n', nerr_vi);

%@T
%
% Now check that I did the calculus right for the gradient and
% Laplacian of the $\Wps$ kernel, the gradient of the pressure kernel,
% and the Laplacian of the viscosity kernel
%@c
h = rand(1);
r = rand(3,1)*h/4;
q = norm(r)/h;
r2 = r'*r;
h2 = h^2;
dr = norm(r)*1e-4;

gWp6_fd = fd_grad(@(r) Wp6(r,h), r, dr);
gWsp_fd = fd_grad(@(r) Wsp(r,h), r,dr);
lWp6_fd = fd_laplace(@(r) Wp6(r,h), r, dr);
lWvi_fd = fd_laplace(@(r) Wvi(r,h), r, dr);

gWp6_ex = -(945/32/pi)/h^5 *(1-q^2)^2 * r;
gWsp_ex = -45/pi/h^5*(1-q)^2/q * r;
lWp6_ex = (945/32/pi)/h^5 * (1-q^2)*(7*q^2-3);
lWvi_ex = 45/pi/h^5 * (1-q);

fprintf('Check kernel derivatives:\n');
fprintf('  grad Wp6:  %g\n', norm(gWp6_fd-gWp6_ex)/norm(gWp6_ex));
fprintf('  grad Wsp:  %g\n', norm(gWsp_fd-gWsp_ex)/norm(gWsp_ex));
fprintf('  lapl Wp6: %g\n', (lWp6_fd-lWp6_ex)/lWp6_ex);
fprintf('  lapl Wvi: %g\n', (lWvi_fd-lWvi_ex)/lWvi_ex);

%@T
%
% Now check that $\fWvi(q)$ satisfies the boundary conditions
% \begin{align*}
%   f(1) &= 0 \\
%   f'(1) &= 0 \\
% \end{align*}
% The first two conditions we check directly.
%@c
dq   = 1e-4;
fprintf('Relerr for viscosity kernel checks:\n');
fprintf('  fvi (1): %g\n', fvi(1) );
fprintf('  dfvi(1): %g\n', fd_deriv(fvi,1,dq) );


%@T
%
% Now, let me check that I did the algebra right in getting the
% condensed formula for the interaction forces.
%@c

% Set up random parameter choices
r_ij = rand(3,1);
v_ij = rand(3,1);
k    = rand(1);
rho0 = rand(1);
rhoi = rand(1);
rhoj = rand(1);
mass = rand(1);
mu   = rand(1);
q    = norm(r_ij)/h;

% Compute pressures via equation of state
Pi = k*(rhoi-rho0);
Pj = k*(rhoj-rho0);

% Differentiate the kernels 
Wsp_x = -45/pi/h^5*(1-q)^2/q * r_ij;
LWvi   = 45/pi/h^5*(1-q);

% Do the straightforward computation
fpressure  = -mass*(Pi+Pj)/2/rhoj * Wsp_x;
fviscous   = -mu*mass*v_ij/rhoj * LWvi;
finteract1 = fpressure + fviscous;

% Do the computation based on my condensed formula
finteract2 = 45*mass/pi/h^5/rhoj * (1-q) * ...
    ( k/2*(rhoi+rhoj-2*rho0)*(1-q)/q * r_ij - mu * v_ij );

% Compare
fprintf('Relerr in interaction force check:\n');
fprintf('  fint:  %g\n', norm(finteract1-finteract2)/norm(finteract1));

%@T
%
% Of course, all the above is supported by a number of little
% second-order accurate finite difference calculations.
%@c

function fp   = fd_deriv(f,r,h)
  fp = ( f(r+h)-f(r-h) )/2/h;

function fpp = fd_deriv2(f,r,h)
  fpp = ( f(r+h)-2*f(r)+f(r-h) )/h/h;

function del2f = fd_laplace_radial(f,r,h)
  del2f = fd_deriv2(f,r,h) + 2*fd_deriv(f,r,h)/r;

function del2f = fd_laplace(f,r,h)
  e1 = [1; 0; 0];
  e2 = [0; 1; 0];
  e3 = [0; 0; 1];
  del2f = (-6*f(r)+...
           f(r+h*e1)+f(r+h*e2)+f(r+h*e3)+...
           f(r-h*e1)+f(r-h*e2)+f(r-h*e3) )/h/h;

function gradf = fd_grad(f,r,h)
  e1 = [1; 0; 0];
  e2 = [0; 1; 0];
  e3 = [0; 0; 1];
  gradf = [f(r+h*e1)-f(r-h*e1); 
           f(r+h*e2)-f(r-h*e2);
           f(r+h*e3)-f(r-h*e3)] / 2 / h;

  
