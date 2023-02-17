%%[V,F] = torus(8);
%%[V,F] = subdivided_sphere(0);
total_area_func = [];
for r = [0.8 0.9]
  N = 200;
  [V,F] = cylinder_mesh(r,N,'Stacks',round(N/pi/r/2));
  
  U = V;
  if isempty(total_area_func)
    total_area_func = per_element_energy(@triangle_area,U,F,'Name','triangle_area','UseFile',true);
  end
  
  O = outline(F); 
  b = unique(O);
  bc = V(b,:);
  U0 = V;
  U = V;
  
  M = repdiag(massmatrix(V,F),size(V,2));
  h = 1/6;
  
  U(b,:) = bc;
  clf;
  CM = oklab2rgb([0.74*ones(1,256);0.125*[cos((1:256)/256*2*pi);sin((1:256)/256*2*pi)]]');
  colormap(CM);
  C = interp1(linspace(caxis*[1;0],caxis*[0;1],size(colormap,1))',colormap,U(:,3));
  tsh = tsurf(F,U,fphong,'FaceVertexCData',C,fsoft,falpha(1,0.1));
  l = add_lights();
  set(gca,'Visible','off','Position',[0 0 1 1]);
  set(gcf,'Color',CM(1,:));
  ssh = add_shadow(tsh,l{5});
  hold on;
  plot_edges(U,O,'-k','LineWidth',2);
  hold off;
  axis equal;
  camproj('persp');
  
  for outer = 1:10
    U1 = U0;
    U0 = U;
    f = @(U) total_area_func(U,F) + ...
        0.5 * 1/h^2*(U(:)- 2*U0(:) + U1(:))'*M*(U(:) - 2*U0(:) + U1(:));
    for newton_iter = 1:1
      [f_a,G_a,H_a] = total_area_func(U,F);
      G = h^2*G_a(:) + M*(U(:)-2*U0(:)+U1(:));
      H = h^2*H_a    + M;
    
      dU = reshape(min_quad_with_fixed( ...
        0.5*H, ...
        G, ...
        reshape(b+[0:size(U,2)-1]*size(U,1),[],1),bc(:)*0),size(U));
      if norm(dU)<1e-5
        break;
      end
      [t,U,fU] = backtracking_line_search(f,U,reshape(G,size(U)),dU,0.3,0.5,30);
      %tsh.Vertices = U;
      %drawnow;
      if t == 0
        break;
      end
      if newton_iter == 100
        warning('newton did not converge');
      end
    end
    tsh.Vertices = U;
    %tsh.FaceVertexCData = C;
    %apply_ambient_occlusion(tsh,'SoftLighting',false,'AddLights',false);
    hold on;
    delete(ssh{1});
    ssh = add_shadow(tsh,l{5});
    hold off;
    drawnow;
  end
end

function f = triangle_area(V)
  N = cross(V(2,:)-V(1,:),V(3,:)-V(1,:),2);
  f = sqrt(sum(N.^2,2))/2;
end
