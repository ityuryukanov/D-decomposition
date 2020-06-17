function stbl=face_check(regs, P, Q, L, R, colors)
stbl=repmat(struct('x',0,'y',0,'hole',0), length(regs), 1);
k=1;
for i=1:1:length(regs)
  % Polygon's bounding box:
  x_min=min(regs(i).x);
  x_max=max(regs(i).x);
  y_min=min(regs(i).y);
  y_max=max(regs(i).y);
  x_offset=1e-3*min(abs(x_max-x_min));
  y_offset=1e-3*min(abs(y_max-y_min));
  Ngrid=7;  % initial density of testing grid NgridxNgrid
  fl=true;
  N_iter=0;
  
  while fl
    % Setup testing grid:
    xgv=linspace(x_min-x_offset, x_max+x_offset, Ngrid);
    ygv=linspace(y_min-y_offset, y_max+y_offset, Ngrid);
    [X,Y] = meshgrid(xgv,ygv);
    X=reshape(X,1,Ngrid*Ngrid);
    Y=reshape(Y,1,Ngrid*Ngrid);
    
    % Check which points are "in":
    IN=inpoly([X;Y], [regs(i).x; regs(i).y]);
    %{
    plot(X(IN), Y(IN), 'b.', X(~IN), Y(~IN), 'r.', regs(i).x, regs(i).y, 'c.', 'Linewidth', 2); hold on;
    plot(regs(i).x, regs(i).y, 'c-', 'Linewidth', 2);
    axis equal
    %}
    
    X=X(IN);
    Y=Y(IN);
    if length(X)>14      
      d=p_poly_dist(X, Y, regs(i).x, regs(i).y);
      [dd, j]=max(d);      
      Np=clp_check_poly(Y(j), X(j), P, Q, L, R);
      %{
      %c=abs(ceil(length(colors)*rand(1)));
      fill(regs(i).x, regs(i).y, colors{Np+1});
      hold on;
      %plot(X(j), Y(j), 'm.', 'MarkerSize', 17);
      text(X(j), Y(j), num2str(Np), 'FontSize', 9, 'Color', 'k',...
        'FontWeight', 'Bold', 'HorizontalAlignment', 'center',...
        'VerticalAlignment', 'middle');  
      %close;
      %}
      fl=false;
      if Np==0
        stbl(k)=regs(i);
        k=k+1;
      end
    else
      Ngrid=floor(1.3*Ngrid); % increase the density
    end
    
    N_iter=N_iter+1;
    if N_iter>20
      error('face_check.m: Iteration limit exceeded. No suitable testing point was found.\nThe tested region is likely to be degenerate.%s', ' ');
    end
  end
  
end
stbl(k:end)=[];

end