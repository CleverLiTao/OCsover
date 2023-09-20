function sol=OCsolver_IntegralConstraint_FixEnd(x0,xgoal,horizon,dyn,cost,constraint_cost,d)
import casadi.*

% setup the dimensions of the state and input
n=dyn.numel_in('X');
m=dyn.numel_in('U');

% setup the bounds of the state and input
controlLowerBound   =   -2000*ones(m,1);
controlUpperBound   =   2000*ones(m,1);
stateLowerBound     =   -1000*ones(n,1);
stateUpperBound     =   1000*ones(n,1);


% Start with an empty NLP
w={};
w0 = [];
lbw = [];
ubw = [];
J0 = 0;
J1 = 0;
cnt={};
lbg = [];
ubg = [];



% "Lift" initial conditions
Xk = MX.sym('x0', n);
w = {w{:}, Xk};
lbw = [lbw; x0];
ubw = [ubw; x0];
w0 = [w0; x0];

% Formulate the NLP for the path
for k=0:horizon
    % New NLP variable for the control
    Uk = MX.sym(['U_' num2str(k)],m);
    w = {w{:}, Uk};
    lbw = [lbw(:); controlLowerBound];
    ubw = [ubw(:); controlUpperBound];
    w0 = [w0(:);  zeros(m,1)];
    
    % Integrate till the end of the interval
    res = dyn('X', Xk, 'U', Uk);
    Xkp = res.f;
    Ck  = cost('X', Xk, 'U', Uk);
    J0=J0+Ck.c;
    Consk = constraint_cost('X', Xk, 'U', Uk);
    J1 = J1 + Consk.cons;
    
    % New NLP variable for state at the next time step
    Xk = MX.sym(['X_' num2str(k+1)], n);
    w = {w{:}, Xk};
    lbw = [lbw(:);  stateLowerBound];
    ubw = [ubw(:);  stateUpperBound];
    w0 = [w0(:); zeros(n,1)];
    
    % Add equality constraint
    cnt = [cnt{:}; Xkp-Xk];
    lbg = [lbg(:);zeros(n,1)];
    ubg = [ubg(:);zeros(n,1)];
end

% adding the final contstraint
cnt=[cnt{:};Xk-xgoal];
lbg = [lbg(:);zeros(n,1)];
ubg = [ubg(:);zeros(n,1)];

% Add integral constraints
cnt = [cnt{:}; J1];
lbg = [lbg(:);0];
ubg = [ubg(:);d];


% Create an NLP solver
prob = struct('f', J0, 'x', vertcat(w{:}), 'g', vertcat(cnt{:})    );
solver = nlpsol('solver', 'ipopt', prob);


% Solve the NLP
solution = solver('x0', w0, 'lbx', lbw, 'ubx', ubw,...
            'lbg', lbg, 'ubg', ubg);

states_inputs=full(solution.x);
states_inputs=reshape(states_inputs(n+1:end),[m+n,horizon+1]);
lambdas=full(solution.lam_g);
beta = lambdas(end);
lambdas=reshape(lambdas(1:end-1),[n,horizon+2]);
g = full(solution.g);
cons = g(end);


sol.u=states_inputs(1:m,:);
sol.x=[x0,states_inputs(m+1:end,:)];
sol.lambda=lambdas;
sol.beta = beta;
sol.cost = solution.f;
sol.constraint = cons;
        
end



