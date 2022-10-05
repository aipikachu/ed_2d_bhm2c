function [ham_mat] = hamiltonian_2d_bhm2c(basis,ham_elems,...
    Jx_up,Jx_dn,Jy_up,Jy_dn,U_uu,U_dd,U_ud,mu_up,mu_dn)


%%
ham_mat = struct();


%%
n_bs = basis.n_bs;
L = basis.L;
Lx = basis.Lx;
Ly = basis.Ly;

BDCx = ham_elems.BDCx;
BDCy = ham_elems.BDCy;


%% boundary condition setting
if strcmp(BDCy,'obc') && strcmp(BDCx,'obc')
    BDC = 'oo';
elseif strcmp(BDCy,'obc') && strcmp(BDCx,'pbc')
    BDC = 'op';
elseif strcmp(BDCy,'pbc') && strcmp(BDCx,'obc')
    BDC = 'po';
else
    BDC = 'pp';
end


%% hubbard parameters check
% 01_1: Jx_up
if numel(Jx_up) == 1
    Jx_up_lt = Jx_up;
    Jx_up_tag = 1;
else
    switch BDC
        case 'oo'
            if sum(size(Jx_up) == [max(Ly-1,1),max(Lx-1,1)],2) == 2
                Jx_up_lt = reshape(Jx_up',1,[]);
                Jx_up_tag = 0;
            else
                error('Error! Invalid input parameter of "Jx_up"!')
            end
        case 'op'
            if sum(size(Jx_up) == [max(Ly-1,1),Lx],2) == 2
                Jx_up_lt = reshape(Jx_up',1,[]);
                Jx_up_tag = 0;
            else
                error('Error! Invalid input parameter of "Jx_up"!')
            end
        case 'po'
            if sum(size(Jx_up) == [Ly,max(Lx-1,1)],2) == 2
                Jx_up_lt = reshape(Jx_up',1,[]);
                Jx_up_tag = 0;
            else
                error('Error! Invalid input parameter of "Jx_up"!')
            end
        case 'pp'
            if sum(size(Jx_up) == [Ly,Lx],2) == 2
                Jx_up_lt = reshape(Jx_up',1,[]);
                Jx_up_tag = 0;
            else
                error('Error! Invalid input parameter of "Jx_up"!')
            end
    end
end


% 01_2: Jx_dn
if numel(Jx_dn) == 1
    Jx_dn_lt = Jx_dn;
    Jx_dn_tag = 1;
else
    switch BDC
        case 'oo'
            if sum(size(Jx_dn) == [max(Ly-1,1),max(Lx-1,1)],2) == 2
                Jx_dn_lt = reshape(Jx_dn',1,[]);
                Jx_dn_tag = 0;
            else
                error('Error! Invalid input parameter of "Jx_dn"!')
            end
        case 'op'
            if sum(size(Jx_dn) == [max(Ly-1,1),Lx],2) == 2
                Jx_dn_lt = reshape(Jx_dn',1,[]);
                Jx_dn_tag = 0;
            else
                error('Error! Invalid input parameter of "Jx_dn"!')
            end
        case 'po'
            if sum(size(Jx_dn) == [Ly,max(Lx-1,1)],2) == 2
                Jx_dn_lt = reshape(Jx_dn',1,[]);
                Jx_dn_tag = 0;
            else
                error('Error! Invalid input parameter of "Jx_dn"!')
            end
        case 'pp'
            if sum(size(Jx_dn) == [Ly,Lx],2) == 2
                Jx_dn_lt = reshape(Jx_dn',1,[]);
                Jx_dn_tag = 0;
            else
                error('Error! Invalid input parameter of "Jx_dn"!')
            end
    end
end


% 01_3: Jy_up
if numel(Jy_up) == 1
    Jy_up_lt = Jy_up;
    Jy_up_tag = 1;
else
    switch BDC
        case 'oo'
            if sum(size(Jy_up) == [max(Ly-1,1),max(Lx-1,1)],2) == 2
                Jy_up_lt = reshape(Jy_up',1,[]);
                Jy_up_tag = 0;
            else
                error('Error! Invalid input parameter of "Jy_up"!')
            end
        case 'op'
            if sum(size(Jy_up) == [max(Ly-1,1),Lx],2) == 2
                Jy_up_lt = reshape(Jy_up',1,[]);
                Jy_up_tag = 0;
            else
                error('Error! Invalid input parameter of "Jy_up"!')
            end
        case 'po'
            if sum(size(Jy_up) == [Ly,max(Lx-1,1)],2) == 2
                Jy_up_lt = reshape(Jy_up',1,[]);
                Jy_up_tag = 0;
            else
                error('Error! Invalid input parameter of "Jy_up"!')
            end
        case 'pp'
            if sum(size(Jy_up) == [Ly,Lx],2) == 2
                Jy_up_lt = reshape(Jy_up',1,[]);
                Jy_up_tag = 0;
            else
                error('Error! Invalid input parameter of "Jy_up"!')
            end
    end
end


% 01_4: Jy_dn
if numel(Jy_dn) == 1
    Jy_dn_lt = Jy_dn;
    Jy_dn_tag = 1;
else
    switch BDC
        case 'oo'
            if sum(size(Jy_dn) == [max(Ly-1,1),max(Lx-1,1)],2) == 2
                Jy_dn_lt = reshape(Jy_dn',1,[]);
                Jy_dn_tag = 0;
            else
                error('Error! Invalid input parameter of "Jy_dn"!')
            end
        case 'op'
            if sum(size(Jy_dn) == [max(Ly-1,1),Lx],2) == 2
                Jy_dn_lt = reshape(Jy_dn',1,[]);
                Jy_dn_tag = 0;
            else
                error('Error! Invalid input parameter of "Jy_dn"!')
            end
        case 'po'
            if sum(size(Jy_dn) == [Ly,max(Lx-1,1)],2) == 2
                Jy_dn_lt = reshape(Jy_dn',1,[]);
                Jy_dn_tag = 0;
            else
                error('Error! Invalid input parameter of "Jy_dn"!')
            end
        case 'pp'
            if sum(size(Jy_dn) == [Ly,Lx],2) == 2
                Jy_dn_lt = reshape(Jy_dn',1,[]);
                Jy_dn_tag = 0;
            else
                error('Error! Invalid input parameter of "Jy_dn"!')
            end
    end
end


% 02_1: U_uu
if numel(U_uu) == 1
    U_uu_lt = U_uu;
    U_uu_tag = 1;
elseif sum(size(U_uu) == [Ly,Lx],2) == 2
    U_uu_lt = reshape(U_uu',1,[]);
    U_uu_tag = 0;
else
    error('Error! Invalid input parameter "U_uu"!')
end

% 02_2: U_dd
if numel(U_dd) == 1
    U_dd_lt = U_dd;
    U_dd_tag = 1;
elseif sum(size(U_dd) == [Ly,Lx],2) == 2
    U_dd_lt = reshape(U_dd',1,[]);
    U_dd_tag = 0;
else
    error('Error! Invalid input parameter "U_dd"!')
end

% 02_3: U_ud
if numel(U_ud) == 1
    U_ud_lt = U_ud;
    U_ud_tag = 1;
elseif sum(size(U_ud) == [Ly,Lx],2) == 2
    U_ud_lt = reshape(U_ud',1,[]);
    U_ud_tag = 0;
else
    error('Error! Invalid input parameter "U_ud"!')
end


% 03_1: mu_up
if numel(mu_up) == 1
    mu_up_lt = mu_up;
    mu_up_tag = 1;
elseif sum(size(mu_up) == [Ly,Lx],2) == 2
    mu_up_lt = reshape(mu_up',1,[]);
    mu_up_tag = 0;
else
    error('Error! Invalid input parameter "mu_up"!')
end

% 03_2: mu_dn
if numel(mu_dn) == 1
    mu_dn_lt = mu_dn;
    mu_dn_tag = 1;
elseif sum(size(mu_dn) == [Ly,Lx],2) == 2
    mu_dn_lt = reshape(mu_dn',1,[]);
    mu_dn_tag = 0;
else
    error('Error! Invalid input parameter "mu_dn"!')
end


%% hamiltonian generation
ham = sparse(n_bs,n_bs);


% 01_1: Onsite interaction, (up,up);
if U_uu_tag
    ham = ham + U_uu_lt * ham_elems.ham_elems_U_uu.sumAllSite;
else
    for kk = 1:L
        field_cur = ['site_',num2str(kk)];
        ham = ham + U_uu_lt(kk) * ham_elems.ham_elems_U_uu.(field_cur);
    end
end


% 01_2: Onsite interaction, (dn,dn);
if U_dd_tag
    ham = ham + U_dd_lt * ham_elems.ham_elems_U_dd.sumAllSite;
else
    for kk = 1:L
        field_cur = ['site_',num2str(kk)];
        ham = ham + U_dd_lt(kk) * ham_elems.ham_elems_U_dd.(field_cur);
    end
end


% 01_3: Onsite interaction, (up,dn);
if U_ud_tag
    ham = ham + U_ud_lt * ham_elems.ham_elems_U_ud.sumAllSite;
else
    for kk = 1:L
        field_cur = ['site_',num2str(kk)];
        ham = ham + U_ud_lt(kk) * ham_elems.ham_elems_U_ud.(field_cur);
    end
end


% 02_1: chemical potential, |up>
if mu_up_tag
    ham = ham + mu_up_lt * ham_elems.ham_elems_mu_up.sumAllSite;
else
    for kk = 1:L
        field_cur = ['site_',num2str(kk)];
        ham = ham + mu_up_lt(kk) * ham_elems.ham_elems_mu_up.(field_cur);
    end
end


% 02_2: chemical potential, |dn>
if mu_dn_tag
    ham = ham + mu_dn_lt * ham_elems.ham_elems_mu_dn.sumAllSite;
else
    for kk = 1:L
        field_cur = ['site_',num2str(kk)];
        ham = ham + mu_dn_lt(kk) * ham_elems.ham_elems_mu_dn.(field_cur);
    end
end


%%
% 03: tunneling xdir, |up> and |dn>
switch BDCx
    case 'obc'
        % 01: for 'Jx_up' term
        if Jx_up_tag
            ham = ham + Jx_up_lt * ham_elems.ham_elems_Jx_up.sumAllSite;            
        else
            for kk = 1:Ly
                for jj = 1:Lx-1
                    ia = (kk-1)*Lx + jj;
                    ib = ia + 1;
                    field_cur = ['site_',num2str(ia),'_',num2str(ib)];
                    
                    ham = ham + Jx_up_lt(kk) * ham_elems.ham_elems_Jx_up.(field_cur);                    
                end
            end
        end
        
        % 02: for 'Jx_dn' term
        if Jx_dn_tag
            ham = ham + Jx_dn_lt * ham_elems.ham_elems_Jx_dn.sumAllSite; 
        else
            for kk = 1:Ly
                for jj = 1:Lx-1
                    ia = (kk-1)*Lx + jj;
                    ib = ia + 1;
                    field_cur = ['site_',num2str(ia),'_',num2str(ib)];
                    
                    ham = ham + Jx_dn_lt(kk) * ham_elems.ham_elems_Jx_dn.(field_cur);                    
                end
            end
        end
    
    case 'pbc'
        % 01: for 'Jx_up' term
        if Jx_up_tag            
            ham = ham + Jx_up_lt * ham_elems.ham_elems_Jx_up.sumAllSite;            
        else
            for kk = 1:Ly
                for jj = 1:Lx
                    ia = (kk-1)*Lx + jj;
                    ib = (kk-1)*Lx + mod(jj,Lx)+1;
                    field_cur = ['site_',num2str(ia),'_',num2str(ib)];
                    
                    ham = ham + Jx_up_lt(kk) * ham_elems.ham_elems_Jx_up.(field_cur);                    
                end
            end
        end
        
        % 02: for 'Jx_dn' term
        if Jx_dn_tag
            ham = ham + Jx_dn_lt * ham_elems.ham_elems_Jx_dn.sumAllSite; 
        else
            for kk = 1:Ly
                for jj = 1:Lx
                    ia = (kk-1)*Lx + jj;
                    ib = (kk-1)*Lx + mod(jj,Lx)+1;
                    field_cur = ['site_',num2str(ia),'_',num2str(ib)];
                    
                    ham = ham + Jx_dn_lt(kk) * ham_elems.ham_elems_Jx_dn.(field_cur);                    
                end
            end
        end
end


% 04: tunneling ydir, |up> and |dn>
switch BDCy
    case 'obc'
        % 01: for 'Jy_up' term
        if Jy_up_tag            
            ham = ham + Jy_up_lt * ham_elems.ham_elems_Jy_up.sumAllSite;            
        else
            for kk = 1:Ly-1
                for jj = 1:Lx
                    ia = (kk-1)*Lx + jj;
                    ib = ia + Lx;
                    field_cur = ['site_',num2str(ia),'_',num2str(ib)];
                    
                    ham = ham + Jy_up_lt(kk) * ham_elems.ham_elems_Jy_up.(field_cur);                    
                end
            end
        end
        
        % 02: for 'Jy_dn' term
        if Jy_dn_tag
            ham = ham + Jy_dn_lt * ham_elems.ham_elems_Jy_dn.sumAllSite; 
        else
            for kk = 1:Ly-1
                for jj = 1:Lx
                    ia = (kk-1)*Lx + jj;
                    ib = ia + Lx;
                    field_cur = ['site_',num2str(ia),'_',num2str(ib)];
                    
                    ham = ham + Jy_dn_lt(kk) * ham_elems.ham_elems_Jy_dn.(field_cur);                    
                end
            end
        end
    
    case 'pbc'
        % 01: for 'Jy_up' term
        if Jy_up_tag            
            ham = ham + Jy_up_lt * ham_elems.ham_elems_Jy_up.sumAllSite;            
        else
            for kk = 1:Ly
                for jj = 1:Lx
                    ia = (kk-1)*Lx + jj;
                    ib = mod(kk,Ly)*Lx + jj;
                    field_cur = ['site_',num2str(ia),'_',num2str(ib)];
                    
                    ham = ham + Jy_up_lt(kk) * ham_elems.ham_elems_Jy_up.(field_cur);                    
                end
            end
        end
        
        % 02: for 'Jy_dn' term
        if Jy_dn_tag
            ham = ham + Jy_dn_lt * ham_elems.ham_elems_Jy_dn.sumAllSite; 
        else
            for kk = 1:Ly
                for jj = 1:Lx
                    ia = (kk-1)*Lx + jj;
                    ib = mod(kk,Ly)*Lx + jj;
                    field_cur = ['site_',num2str(ia),'_',num2str(ib)];
                    
                    ham = ham + Jy_dn_lt(kk) * ham_elems.ham_elems_Jy_dn.(field_cur);                    
                end
            end
        end
end


%%
ham_mat.ham_elems = ham_elems;
ham_mat.ham = ham;


