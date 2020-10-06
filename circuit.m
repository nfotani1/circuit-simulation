classdef circuit < handle
    % V1.3: added delete function (destructor) which closes the log file
    % V1.2: add I_NL to keep track of currents through nonlinear components
    % Note: V1.1 was incorrect--had the init. conds. calc. flaw.
    
    properties
        N_nodes; % number of nodes, including the ground node
        alpha; % 0.5 for trapezoidal method; 1.0 for backward Euler
        Dt; % timestep
        
        GR; % admittance matrix for the resistors
        
        GC; % admittance matrix for the capacitors
        AC; % incidence matrix for the capacitors
        V_C_init; % Iniial voltages across the capacitors
        reduced_AC; % capacitor incidence matrix w/o cycles
        reduced_V_C_init; % remove corresponding in initial voltages
        
        AE; % incidence matrix for the voltage sources
        E; % column vector containing voltage source voltages
        E_func; % column cell array holding handles to voltage source functions of time
        
        I_nonlin_fns; % cell array: each column contains handles to functions
        % defining one nonlinear component
        I_nonlin_fn_groups; % cell array table: each row contains a different type
        % of nonlinear component, followed by a row array of columns
        % in I_nonlin_fns and ANL that contain that type of nonlinear
        % component
        q_array; % each column contain auxiliary variables for one nonlinear component
        ANL; % incidence matrix for the nonlinear components in the circuit
        I_NL; % column vector of currents through the nonlinear components
        
        GA; % left-hand side matrix
        V; % column vector containing nodal voltages and some currents
        GB; % right-hand side matrix
        
        t; % current time
        
        fp_log; % file pointer of a log
    end
    
    methods
        function c = circuit(N_nodes,alpha,Dt)
            
            fprintf('Circuit class V1.3: contains I_NL\n');
            c.N_nodes = N_nodes;
            c.alpha = alpha;
            c.Dt = Dt;
            
            % Note: Dimensions of the following may change
            % as the code is run:
            
            c.GR = sparse(N_nodes,N_nodes);
            
            c.GC = sparse(N_nodes,N_nodes);
            c.AC = sparse(N_nodes,0);
            c.V_C_init = zeros(0,1);
            
            c.AE = sparse(N_nodes,0);
            c.E_func = cell(0,1);
            c.E = zeros(N_nodes,1);
            
            c.ANL = sparse(N_nodes,0);
            % c.I_nonlin_fns = cell(0,3);
            
            c.V = zeros(N_nodes,1);
            
            c.t = 0.0;
            c.fp_log = fopen('circuit.log','w');
            
        end
        
        function delete(c)
            fclose(c.fp_log);
        end
        
        function add_resistor(c,node1,node2,R)
            % Add a resistor with resistance R between nodes node1 and
            % node2 by adding the resistor to the resistor admittance
            % matrix:
            
            c.GR(node1,node1) = c.GR(node1,node1) + 1/R;
            c.GR(node1,node2) = c.GR(node1,node2) - 1/R;
            c.GR(node2,node1) = c.GR(node2,node1) - 1/R;
            c.GR(node2,node2) = c.GR(node2,node2) + 1/R;
            
        end
        
        function add_capacitor(c,node1,node2,C,V_C_initial)
            % Add a capacitor with capacitance C between nodes
            % node1 and node2 by adding the capacitor to the
            % capacitance admittance matrix:
            
            c.GC(node1,node1) = c.GC(node1,node1) + C/c.Dt;
            c.GC(node1,node2) = c.GC(node1,node2) - C/c.Dt;
            c.GC(node2,node1) = c.GC(node2,node1) - C/c.Dt;
            c.GC(node2,node2) = c.GC(node2,node2) + C/c.Dt;
            
            % Also construct an admittance matrix for the capacitors,
            % for use in determining the initial conditions:
            AC_cols = size(c.AC,2);
            c.AC(node1,AC_cols+1) = 1;
            c.AC(node2,AC_cols+1) = -1;
            
            % Add the initial voltage to V_C_init:
            c.V_C_init(end+1,1) = V_C_initial;
            
        end
        
        function add_voltage_source(c,node1,node2,E_of_t_handle)
            % Add a voltage source between node1 and node2 by adding
            % the component to the voltage source incidence matrix AE.
            % Store the the handle of the function defining the
            % time-varying voltage, E_of_t_handle, in the cell array
            % E_func.
            
            % Add the voltage source to the voltage source incidence
            % matrix:
            AE_cols = size(c.AE,2);
            c.AE(node1,AE_cols+1) = 1;
            c.AE(node2,AE_cols+1) = -1;
            
            % Add the function handle of the function
            % defining the voltage source voltage vs. time
            % to the list of voltage source functions:
            c.E_func{end+1} = E_of_t_handle;
            
        end
        
        function add_I_nonlin(c,node1,node2,f_handle,g_handle,Df_handle,q_init_conds)
            % Add a nonlinear component directing current from node1 to node2
            % defined by the functions pointed to f_handle (which calculates
            % the current), g_handle (which updates the auxiliary variables q),
            % and Df_handle (which calculates the partial derivatives,
            % including df/dV), with initial conditions for q stored in the -row-
            % array, q_init_conds.
            
            % f_handle should be a handle to function of the form,
            % f(V,q,params), where V can be a column vector of several
            % component voltages (allows for vectorization),
            % q is an array, each of whose rows corresponds to one component,
            % and whose columns correspond to different auxiliary
            % variables, and p is a structure of parameters.
            % Similarly for g_handle.
            % Df_handle should point to a function of the form,
            % [df/dv,df/dq_1,...df/dq_m] = Df(V,q,params)
            % Note that df/dV is the first column in the returned matrix.
            
            % The length of the -row- array q_init_conds
            % also tells the routines how many auxiliary variables there are
            % for this nonlinear component and for corresponding the functions
            % f, g and Df.
            %
            % Each row of the q_array contains the auxiliary variables
            % for one nonlinear component.
            % The rows of auxiliary variables are stored in the same order
            % as the nonlinear components were added.
            
            % I_nonlin_fns is a table, stored as a cell array.
            % Each -row- in the table contains the function handles for
            % the f, g, and Df functions for one nonlinear component:
            c.I_nonlin_fns = [c.I_nonlin_fns; {f_handle,g_handle,Df_handle}];
            n_comp = size(c.I_nonlin_fns,1);
            
            % Add a -column- to the nonlinear component incidence
            % matrix for this component:
            c.ANL(node1,n_comp) = 1;
            c.ANL(node2,n_comp) = -1;
            
            % I_nonlin_fn_groups is another table stored as a cell array.
            % Each row of this table contains a unique combination
            % of f, g and Df functions, followed by a row array of all
            % the row numbers in I_nonlin_fns in which that combination
            % of f, g and Df occurs.
            % Here, we update this table for the nonlinear component we
            % are adding:
            N_groups = size(c.I_nonlin_fn_groups,1);
            if (N_groups==0)
                c.I_nonlin_fn_groups = {f_handle,g_handle,Df_handle,n_comp};
            else
                for i_group = 1:N_groups
                    % If f, g and Df for this nonlinear element match an existing pair
                    % in I_nonlin_fn_groups, just add the column number of the
                    % entry just added to I_nonlin_fns:
                    if (isequal(c.I_nonlin_fn_groups{i_group,1},f_handle) && ...
                            isequal(c.I_nonlin_fn_groups{i_group,2},g_handle) && ...
                            isequal(c.I_nonlin_fn_groups{i_group,3},Df_handle))
                        c.I_nonlin_fn_groups{i_group,4} = [c.I_nonlin_fn_groups{i_group,4},n_comp];
                        break;
                        % Otherwise, no matching entry is found, add a new entry
                        % to I_nonlin_fn_groups:
                    elseif (i_group == N_groups)
                        c.I_nonlin_fn_groups = [c.I_nonlin_fn_groups;{f_handle,g_handle,Df_handle,n_comp}];
                    end
                end
            end
            
            % Store the initial conditions for the auxiliary variables for
            % the new component as a new row in q_array.
            % If the number of initial conditions is not equal to the
            % number of rows in q_array, pad either q_array or the
            % initial conditions with NaN's:
            no_of_q_init_conds = length(q_init_conds);
            [n_rows,n_cols] = size(c.q_array);
            if no_of_q_init_conds > n_cols
                c.q_array = [c.q_array,nan(n_rows,no_of_q_init_conds-n_cols)];
            elseif length(q_init_conds)<n_cols
                q_init_conds = [q_init_conds,nan(1,n_cols-no_of_q_init_conds)];
            end
            c.q_array = [c.q_array; q_init_conds];
            
        end
        
        function prepare_matrices(c)
            
            % Remove any capacitors from the graph of capacitors
            % and voltage sources that form a cycle.
            % Also remove the initial voltages on these capacitors
            % as redundant or inconsistent.
            % Record in the log if this happens:
            AE_cols = size(c.AE,2);
            AC_cols = size(c.AC,2);
            c.reduced_AC = sparse(c.N_nodes,0);
            c.reduced_V_C_init = [];
            creates_cycle(1); % reset cycle seeking
            for ie = 1:AE_cols
                cycle_created = creates_cycle (0,c.AE(:,ie));
                if (cycle_created)
                    fprintf(c.fp_log,'ERROR: Voltage source %i forms a cycle with',ie);
                    fprintf(c.fp_log,' previously defined voltage sources.\n');
                end
            end
            for ic = 1:AC_cols
                cycle_created = creates_cycle (0,c.AC(:,ic));
                if (cycle_created)
                    fprintf(c.fp_log,'Capacitor %i forms a cycle with',ic);
                    fprintf(c.fp_log,' previously defined capacitors and voltage sources.\n');
                    fprintf(c.fp_log,' So the initial voltage value %f ignored\n',c.V_C_init(ic));
                else
                    c.reduced_AC = [c.reduced_AC,c.AC(:,ic)];
                    c.reduced_V_C_init(end+1,1) = c.V_C_init(ic);
                end
            end
            
            % Remove the row and/or column associated with electrical
            % ground from all the relevant matrices and vectors:
            c.GR = c.GR(1:(end-1),1:(end-1));
            c.GC = c.GC(1:(end-1),1:(end-1));
            c.AE = c.AE(1:(end-1),:);
            c.AC = c.AC(1:(end-1),:);
            c.reduced_AC = c.reduced_AC(1:(end-1),:);
            c.ANL = c.ANL(1:(end-1),:);
            
            % Construct the matrices required to march forwards in time:
            N_E_sources = size(c.AE,2);
            c.GA = [c.GC + c.alpha*c.GR, c.alpha*c.AE;...
                (c.AE)',             zeros(N_E_sources,N_E_sources)];
            c.GB = [c.GC - (1-c.alpha)*c.GR, -(1-c.alpha)*c.AE;...
                zeros(N_E_sources,c.N_nodes-1), zeros(N_E_sources,N_E_sources)];
            
            c.E = zeros(c.N_nodes-1+N_E_sources,1);
            
        end
        
        function I = I_nonlin(c,params)
            % Calculate the currents flowing through
            % all the nonlinear elements by evaluating the functions
            % whose handles are stored in the cell array I_nonlin_fns{:,1}.
            % The functions are assumed to be of the form f(V_drop,q,params),
            % where V_drop is a -column- array voltage drops across the components,
            % q is a set of row vectors of auxiliary variables,
            % one for each nonlinear component, and
            % params is a set of parameters used in defining the the function f.
            
            N_comps = size(c.I_nonlin_fns,1);
            N_groups = size(c.I_nonlin_fn_groups,1);
            
            % Vectorized, group version:
            
            I = zeros(N_comps,1);
            
            for i_group = 1:N_groups
                comp_indices = c.I_nonlin_fn_groups{i_group,4};
                V_comp = c.ANL(:,comp_indices)'*c.V(1:(c.N_nodes-1));
                I(comp_indices) = c.I_nonlin_fn_groups{i_group,1}(V_comp,c.q_array(comp_indices,:),params);
            end
            
            c.I_NL = I; % also store the nonlinear currents here, for ease of access
            
        end
        
        function dIdV = dI_nonlin_dV(c,params)
            % Calculate partial derivatives of I_nonlin with respect to V
            % for use in calculating the initial conditions.
            
            N_comps = size(c.I_nonlin_fns,1);
            N_groups = size(c.I_nonlin_fn_groups,1);
            
            % Vectorized, group version:
            
            dIdV = zeros(N_comps,1);
            
            for i_group = 1:N_groups
                comp_indices = c.I_nonlin_fn_groups{i_group,4};
                V_comp = c.ANL(:,comp_indices)'*c.V(1:(c.N_nodes-1));
                DI = c.I_nonlin_fn_groups{i_group,3}(V_comp,c.q_array(comp_indices,:),params);
                dIdV(comp_indices) = DI(:,1); % dI/dV should be in the first column
            end
            
            
        end
        
        function DI = DI_nonlin(c,params)
            % Calculate ALL partial derivatives of I_nonlin.
            
            N_comps = size(c.I_nonlin_fns,1);
            N_groups = size(c.I_nonlin_fn_groups,1);
            N_q = size(c.q_array,2);
            
            % Vectorized, group version:
            
            DI = zeros(N_comps,1+N_q);
            
            for i_group = 1:N_groups
                comp_indices = c.I_nonlin_fn_groups{i_group,4};
                V_comp = c.ANL(:,comp_indices)'*c.V(1:(c.N_nodes-1));
                DI(comp_indices,:) = c.I_nonlin_fn_groups{i_group,3}(V_comp,c.q_array(comp_indices,:),params);
            end
            
        end
              
        function calc_initial_conditions(c,params,t0,V_augmented_guess)
            err_tol = 1.e-7;
            max_iter = 400;
            AE_cols = size(c.AE,2);
            AC_cols = size(c.reduced_AC,2);
            if (nargin<=2)
                t0 = 0;
            end
            if (nargin==4)
                V_augmented = V_augmented_guess;
            else
                V_augmented = zeros(c.N_nodes-1 + AE_cols + AC_cols,1);
            end
            
            EE = zeros(length(c.E_func),1);
            for i = 1:length(c.E_func)
                EE(i) = c.E_func{i}(t0);
            end
            
            fprintf('Circuit V1.2: Correction to calc_initial_condition_in place\n');
            % Newton-Raphson iteration loop:
            for j = 1:max_iter
                INL = c.I_nonlin(params);
                nNL_comps = length(INL);
                diag_dINL_dV = spdiags(c.dI_nonlin_dV(params),0,nNL_comps,nNL_comps);
                G_augmented = [c.GR+c.ANL*diag_dINL_dV*(c.ANL)',c.AE,c.reduced_AC;...
                    c.AE',zeros(AE_cols,AE_cols+AC_cols);...
                    c.reduced_AC',zeros(AC_cols,AE_cols+AC_cols)];
                f = [ [c.GR,c.AE,c.reduced_AC]*V_augmented + c.ANL*INL;...
                    c.AE'*V_augmented(1:(c.N_nodes-1))- EE;...
                    c.reduced_AC'*V_augmented(1:(c.N_nodes-1))-c.reduced_V_C_init];
                delta_V_augmented = - G_augmented \ f;
                V_augmented = V_augmented + delta_V_augmented;
                c.V = V_augmented(1:((c.N_nodes-1)+AE_cols)); 
                % fprintf('Iter no. %i: Error = %e\n',j,norm(delta_V_augmented));
                if (norm(delta_V_augmented)<err_tol)
                    fprintf('Initial conditions converged after %i iterations.\n',j);
                    break;
                end
                if (j==max_iter)
                    fprintf('Warning: Initial conditions did not converge after %i iterations\n',...
                        max_iter);
                end
            end
      
            
            c.t = t0;
        end
        
        function q_advance(c,params)
            % Advance all the auxilary dynamical variables associated
            % with the nonlinear currents.
                        
            N_groups = size(c.I_nonlin_fn_groups,1);
            
            % Vectorized group version:
            
            for i_group = 1:N_groups
                comp_indices = c.I_nonlin_fn_groups{i_group,4};
                V_comp = c.ANL(:,comp_indices)'*c.V(1:(c.N_nodes-1));
                c.q_array(comp_indices,:) = c.I_nonlin_fn_groups{i_group,2}(V_comp,c.q_array(comp_indices,:),params);
            end
            
        end
        
        function advance_circuit(c,params)
            % Advance the circuit from time t to t+Dt:
            AE_cols = size(c.AE,2);
            for i = 1:length(c.E_func)
                c.E(c.N_nodes+i-1) = c.E_func{i}(c.t+c.Dt);
            end
            INL = c.I_nonlin(params);
            c.V = c.GA \ (c.GB*c.V + c.E - [c.ANL*INL;zeros(AE_cols,1)]);
            c.q_advance(params);
            c.t = c.t + c.Dt;
        end
        
    end
    
end


