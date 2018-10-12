function [V_a, V, success, cvx_status_VSDP] = getVSDP(Yp, Yq, P_l, Q_l, P_pv, MagPcc, ...
    Nnode, Node_PV, Node_noPV, Node_batt, ...
    Ninverters, alpha, PB,Q)

cvx_begin sdp quiet

variable V(Nnode,Nnode) hermitian

% Constant
minimize  real(trace(Yp(:,1:Nnode)*V))

subject to

% PCC
indeces = 1:Nnode;
real(trace(Yp(:,indeces)*V)) <= 300000;
real(trace(Yp(:,indeces)*V)) >= -300000;
real(trace(Yq(:,indeces)*V)) <= 3000000;
real(trace(Yq(:,indeces)*V)) >= -300000;


% nodes with inverters
for nn = 1:Ninverters
    node = Node_PV(nn);
    indeces = (node-1)*Nnode+1:node*Nnode;
    batt_idx = find(Node_batt == node);
    if(batt_idx > 0) % battery here
        real(trace(Yp(:,indeces)*V)) + P_l(node) ...
            - (1-alpha(nn))*P_pv(nn) + PB(batt_idx) == 0;
    else
        real(trace(Yp(:,indeces)*V)) + P_l(node) ...
            - (1-alpha(nn))*P_pv(nn) == 0;
    end
     real(trace(Yq(:,indeces)*V)) + Q_l(node) - Q(nn) == 0;

end

% nodes with no inverters
for nn = 1:length(Node_noPV)
    node = Node_noPV(nn);
    indeces = (node-1)*Nnode+1:node*Nnode;
    batt_idx = find(Node_batt == node);
    if(batt_idx > 0) % battery here
        real(trace(Yp(:,indeces)*V)) + P_l(node) + PB(batt_idx) == 0;
    else
        real(trace(Yp(:,indeces)*V)) + P_l(node) == 0;
    end
    real(trace(Yq(:,indeces)*V)) + Q_l(node) == 0;

end

% V has to be positive semidefinite
V >= 0;

% voltage magnitude at the PCC
V(1,1) == (MagPcc)^2;

cvx_end
cvx_status_VSDP = cvx_status;

    [U,L] = eig(V);
    if sum(sum(L > 0.01)) == 1 % the rank of V is 1
        % problem solved
        [sigma,ii] = sort(diag(L),'descend');
        success = 1;
        V_vec_no_service = sqrt(abs(sigma(1))).*U(:,ii(1));
        V_a = diag(sqrt(V));
    else
        success = 0;
        V_a = 0;
    end


end