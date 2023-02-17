function func = per_element_energy(phi,x,E,varargin)
  % 
  % Inputs:
  %   phi  function with prototype
  %      phi(x(E(1,:),:,…),varargin{1}(1,:,…),…)
  %   x  n by ... list of nodal data
  %   E  m by k list of element indices
  %  


  name = tempname;
  constants = {};
  nodal = [];
  use_file = true;
  % Map of parameter names to variable names
  params_to_variables = containers.Map( ...
    {'Name','Constants','Nodal','UseFile'}, ...
    {'name','constants','nodal','use_file'});
  v = 1;
  while v <= numel(varargin)
    param_name = varargin{v};
    if isKey(params_to_variables,param_name)
      assert(v+1<=numel(varargin));
      v = v+1;
      % Trick: use feval on anonymous function to use assignin to this workspace
      feval(@()assignin('caller',params_to_variables(param_name),varargin{v}));
    else
      error('Unsupported parameter: %s',varargin{v});
    end
    v=v+1;
  end
  if ~iscell(constants) && isnumeric(constants)
    constants = {constants};
  end
  if isempty(nodal) 
    nodal = false(size(constants));
  end

  [m,k,n,size_x] = sizes(x,E);

  % Important that this is a row vector so that functions are vectorized
  sym_x_vec = sym('x',[1 prod(size_x)],'real');

  % set up symbolic per-element constant parameters
  sym_c_vec = {};
  sym_c = {};
  for ci = 1:numel(constants)
    % check that we have one (row of) constants per element (or total and
    % we'll broadcast)
    size_ci = size(constants{ci});
    if nodal(ci)
      assert(size_ci(1) == n || size_ci(1) == 1);
      size_ci(1) = k;
    else
      assert(size_ci(1) == m || size_ci(1) == 1);
      size_ci(1) = 1;
    end
    sym_ci = sym(['c_' num2str(ci)],size_ci,'real');
    sym_c = {sym_c{:} sym_ci};
    sym_c_vec = {sym_c_vec{:} reshape(sym_ci,1,[])};
  end

  sym_f = phi(reshape(sym_x_vec,size_x),sym_c{:});
  sym_dfdx = gradient(sym_f,sym_x_vec);
  %sym_d2fdx2 = hessian(sym_f,sym_x_vec);
  % Similar to `hessian` but vector output
  hess = @(sf,sX) cell2sym(arrayfun(@(g) gradient(g,sX),gradient(sf,sX),'UniformOutput',false));
  sym_d2fdx2 = hess(sym_f,sym_x_vec);
  
  
  aux = [name '_sym.m'];
  if use_file
    fprintf('generating file...\n');
    m_f = matlabFunction(sym_f,sym_dfdx,sym_d2fdx2,'vars',{sym_x_vec,sym_c_vec{:}},'File',aux);
    if patch_final_reshape(aux);
      warning('patch final reshape detected and applied...');
    end
  else
    fprintf('generating m_f...\n');
    m_f = matlabFunction(sym_f,'vars',{sym_x_vec,sym_c_vec{:}});
    fprintf('generating m_dfdx...\n');
    m_dfdx = matlabFunction(sym_dfdx,'vars',{sym_x_vec,sym_c_vec{:}});
    fprintf('generating m_d2fdx2...\n');
    m_d2fdx2 = matlabFunction(sym_d2fdx2,'vars',{sym_x_vec,sym_c_vec{:}});
  end
  func = @scatter_gather;

  psd_project = true;

  function [m,k,n,size_x] = sizes(x,E)
    % Number of elements
    m = size(E,1);
    % Size of elements
    k = size(E,2);
    % Number of nodes
    n = size(x,1);
    % set up symbolic variables of differentiation
    size_x = size(x);
    size_x(1) = k;
  end

  function [f,G,H] = scatter_gather(x,E,varargin)
    % Recompute these rather than capture
    [m,k,n,size_x] = sizes(x,E);

    % scattering 
    nc = prod(size_x(2:end));
    I = sub2ind(size(x),repmat(E,1,nc),repmat(repelem(1:nc,k),m,1));

    effective_constants = {};
    for ci = 1:numel( varargin)
      if nodal(ci)
        constant_ci = varargin{ci}(I);
      else
        constant_ci = varargin{ci};
      end
      effective_constants = {effective_constants{:} constant_ci};
    end
    
    if use_file
      switch nargout
      case {0,1}
        [f] = m_f(x(I),effective_constants{:});
      case 2
        [f,dfdx] = m_f(x(I),effective_constants{:});
      case 3
        [f,dfdx,d2fdx2] = m_f(x(I),effective_constants{:});
      end
    else
      f = m_f(x(I),effective_constants{:});
      if nargout>1
        dfdx = m_dfdx(x(I),effective_constants{:});
        if nargout>2
          d2fdx2 = m_d2fdx2(x(I),effective_constants{:});
        end
      end
    end
    f = sum(f);
    if nargout>1
      G = reshape(full_sparse(I,ones(size(I)),reshape(dfdx,size(I)),numel(x),1),size(x));

      if nargout>2
        d2fdx2 = reshape(d2fdx2,size(I,1),[]);
        if psd_project
          d2fdx2 = psd_project_rows(d2fdx2);
        end
        HI = repmat(I,[1 size(I,2)]);
        HJ = repelem(I,1,size(I,2)); 
        H = sparse(HI,HJ,d2fdx2,numel(x),numel(x));
      end
    end
  end

end
