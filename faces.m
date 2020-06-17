function regs=faces(edges,k1_min, k1_max, k2_min, k2_max, colors)
% This function finds all faces of a planar embedding of a graph.
% nargcheck and other input checks can be included as well.

regs=repmat(struct('x',0,'y',0,'hole',0), length(edges), 1);

% BE(1:2,:) - coordinates of edges' endpoints
% BE(3,:) - edges' numbers
% BE(4,:) - the "Begin->End was traversed" flag of an edge
% BE(5,:) - the "End->Begin was traversed" flag of an edge
BE=zeros(5,2*length(edges));   % BE="begin-end". This is an array of BE-coordinates of all edges
for i=1:1:length(edges)
    BE(1:2,2*i-1)=[edges(i).x(1); edges(i).y(1)];
    BE(1:2,2*i)=[edges(i).x(end); edges(i).y(end)];
    BE(3,2*i-1)=i;
    BE(3,2*i)=i;
end

% % Connect loops (i.e. segments with BEGIN==END):
k=1;   % index for (found) closed contours
% for j=1:2:length(BE(1,:))
%     if isequal([BE(1,j),BE(2,j)],[BE(1,j+1),BE(2,j+1)])
%         regs(k).x=edges(BE(3,j)).x;
%         regs(k).y=edges(BE(3,j)).y;
%         CCW = sum(diff(regs(k).x).*(regs(k).y(1:end-1)+regs(k).y(2:end))) < 0; % clockwise-test
%         if CCW
%             BE(5,j)=1;   % if B->E is counterclockwise, check "E->B was traversed", since E->B is then clockwise!
%             BE(5,j+1)=1; % the cycle edge was "checked"  
%         else
%             BE(4,j)=1;   % if B->E is clockwise, check "B->E was traversed"!
%             BE(4,j+1)=1; % the cycle edge was "checked"   
%         end
%         k=k+1;
%     end
% end

% Connect circles: %%???
% Nothing is done now, but later it's needed to divide the initial list of
% edges to the lists of edges of connected (or even 2-connected) components.
% And then apply loops detection and faces traversal to every component.
% If a component consists of only one edge, it's a circle.

% Connect closed contours:
% "Open contours" are impossible in D-decomposition (if everything was done right)

j=1;
% j=find(~(BE(4,:)|BE(5,:)),1);
% if isempty(j)
%     j=find(~(BE(4,:)&BE(5,:)),1);
%     if (BE(4,j)==1)&&(mod(j,2)==0)
%         j=j-1;   % go different direction than before
%     elseif (BE(5,j)==1)&&(mod(j,2)==1)
%         j=j+1;   % go different direction than before
%     end
% end
Niter_max=1000;      % since we use 2 while cycles, a "hang-protection" was implemented
N_unseen=Inf;     % initially set N_unseen>0 (any integer)
while N_unseen>0
    start=BE(3,j);      % "starting" segment
    regs(k).x=edges(start).x;
    regs(k).y=edges(start).y;
    fl_E=false; fl_B=false;
    if (mod(j,2)==0)
        BE(4,j)=1;      % the starting edge was "checked"
        BE(4,j-1)=1;    
        fl_E=true;      % contour grows at its END 
    else
        BE(5,j)=1;      % the starting edge was "checked"
        BE(5,j+1)=1;    
        fl_B=true;      % contour grows at its BEGIN
    end

    Niter=0;
    while true
        a=(BE(1,:)==BE(1,j))&(BE(2,:)==BE(2,j))&(BE(3,:)~=BE(3,j));
        j_fnd=find(a);
                
        % Choose the next segment of the same face:
        if length(j_fnd)>1
            vec_prev=[edges(BE(3,j)).x; edges(BE(3,j)).y];
            if (mod(j,2)==0)
                END=size(vec_prev,2);
                vec_prev=fliplr(vec_prev(:,END-1:END));
            elseif (mod(j,2)==1)
                vec_prev=vec_prev(:,1:2);
            end
            vec_prev=vec_prev(:,2)-vec_prev(:,1); % vector coordinates
            vec_next=zeros(2*length(j_fnd),1);
            for i=1:1:length(j_fnd)
                vec=[edges(BE(3,j_fnd(i))).x; edges(BE(3,j_fnd(i))).y];
                if (mod(j_fnd(i),2)==0)
                    END=size(vec,2);
                    vec=fliplr(vec(:,END-1:END));
                elseif (mod(j_fnd(i),2)==1)
                    vec=vec(:,1:2);
                end
                vec=vec(:,2)-vec(:,1);
                vec_next(2*i-1:2*i,:)=vec;
            end
            Angle=zeros(1,length(j_fnd));
            for i=1:1:length(j_fnd)
                vec=vec_next(2*i-1:2*i,:);
                Angle(i) = atan2(vec_prev(1)*vec(2)-vec(1)*vec_prev(2),vec_prev(1)*vec(1)+vec_prev(2)*vec(2));
                Angle(i) = mod(180/pi * Angle(i), 360);
            end
            [ang, i]=min(Angle);
            j_fnd=j_fnd(i);
        end
        
        if BE(3,j_fnd)==start
            if ~isequal([regs(k).x(1);regs(k).y(1)], [regs(k).x(end);regs(k).y(end)])
                error('ERROR(faces.m): The contour was closed in a wrong way.\nTerminating...%s',1);
            end
            c1=(regs(k).x==k2_min)|(regs(k).x==k2_max);
            c2=(regs(k).y==k1_min)|(regs(k).y==k1_max);
            if all(c1|c2)&&(length(regs(k).x)~=5) %%??? IF EVEN CRB IS TOTALLY INSIDE, there will be AT LEAST 2 RRBS,
                regs(k)=[]; %%??? SO THAT THERE WILL be a connection to the "plotting window"
                k=k-1;
            end
            k=k+1;   % take another contour
            if N_unseen>0
                j=find(~(BE(4,:)|BE(5,:)),1);   % the start of a new contour
                if isempty(j)
                    j=find(~(BE(4,:)&BE(5,:)),1);
                    if (BE(4,j)==1)&&(mod(j,2)==0)
                        j=j-1;   % go different direction than before
                    elseif (BE(5,j)==1)&&(mod(j,2)==1)
                        j=j+1;   % go different direction than before
                    end
                end
            end
            break   % contour is closed
        end
        
        % Connect the found edge to the contour:
        aX=edges(BE(3,j_fnd)).x;
        aY=edges(BE(3,j_fnd)).y;
        if fl_B
            if isequal([aX(end);aY(end)],[regs(k).x(1);regs(k).y(1)])
                regs(k).x=[aX(1:end-1), regs(k).x];
                regs(k).y=[aY(1:end-1), regs(k).y];
            elseif isequal([aX(1);aY(1)],[regs(k).x(1);regs(k).y(1)])
                regs(k).x=[fliplr(aX(2:end)), regs(k).x];
                regs(k).y=[fliplr(aY(2:end)), regs(k).y];
            else
                error('ERROR(faces.m): Something went wrong, the segments don''t match');
            end
        end
        if fl_E
            if isequal([aX(1);aY(1)],[regs(k).x(end);regs(k).y(end)])
                regs(k).x=[regs(k).x, aX(2:end)];
                regs(k).y=[regs(k).y, aY(2:end)];
            elseif isequal([aX(end);aY(end)],[regs(k).x(end);regs(k).y(end)])
                regs(k).x=[regs(k).x, fliplr(aX(1:end-1))];
                regs(k).y=[regs(k).y, fliplr(aY(1:end-1))];
            else
                error('ERROR(faces.m): Something went wrong, the segments don''t match');
            end
        end
        if mod(j_fnd,2)==0
            BE(5,j_fnd)=1;
            j_fnd=j_fnd-1;   % go to the same segment's another end
            BE(5,j_fnd)=1;
        elseif mod(j_fnd,2)==1
            BE(4,j_fnd)=1;
            j_fnd=j_fnd+1;   % go to the same segment's another end
            BE(4,j_fnd)=1;
        end
        N_unseen=sum(BE(4,:)==0)+sum(BE(5,:)==0);
        
        j=j_fnd;   % take the next vertex
%         plot(regs(k).x, regs(k).y, '.-', 'Color', colors{k}, 'LineWidth', 2); hold on;
        
        Niter=Niter+1;   % check the maximum iteration limit
        if Niter>Niter_max
            error('ERROR(faces.m): Maximal iteration number was exceeded. The contours were not identified! Terminating...')
        end
    end
end
regs(k:end)=[];
end