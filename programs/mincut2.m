function [mv, mw] = mincut2(g, s, t)
%MINCUT performs Min Cut algorithm
%   input -
%           g       symetric adjacent matrix of graph. gi,j is weight of
%                   the edge connecting nodes i,j.
%           s       a list of nodes that are forced to be kept in one side
%                   of the cut. default is node 1
%           t       a lust of nodes that are forced to be in the other side
%                   of the cut. 
%   output -
%           mv      mincut belonging vector of nodes, have the same value
%                   belong the same partition.
%           mw      mincut weight.
% Note:     The function ask for the connectivity of graph !!
% Example:
%           g = zeros(7,7);  % creat a graph
%           g(1,[2,4]) = [3 2];  g(2,[3,4]) = [2 1];
%           g(3,[4,6]) = [2 2];  g(4,5) = 1;
%           g(5,[6 7]) = [3 2];  g(6,7) = 1;
%           g = g + g';      % keep symmetry
%           [mv, mw] = mincut(g, [6 7]);
%   result: mv = [1 1 1 1 0 0 0]' , mw = 3
%           means the graph was split to two groups: {1,2,3,4}, {5,6,7}
%           minimum cut weight value is 3
% Complexity:
%           O(|V|*|E|+(log|V|)*|V|^2)
% Reference:
%           M. Stoer and F. Wagner. "A Simple Min Cut Algorithm". 1997
% Creat by: HuYo, 2011-07-13

    
    g = double(g);      N = length(g);    % number of nodes
    
    
    if nargin == 1,     s = 1;  t = 2;     end
    if nargin == 2,     t = setdif(1:N, s);     t = t(1);   end

    
    % remove self edges and inf-value
    g(1:N+1:end) = 0;   
    g(isinf(g))  = 0; 
    
    % predefine
    mv = (1:N)';        mw = inf;

    % merge all source nodes to one, so they'll be unbreakable
    st = [s,t];
    st = sort(st);
    ni = 1:N;
    
    s = sort(s);  
    t = sort(t); 
    u = max(find(st==s(1)), find(st==t(1))) +1;
    
    % merge all sink and source nodes into 2 unbreakable vertices (source
    % and sink)
    
    for n = length(st):-1:u            % descending order is VITAL!!!
        if any(s == st(n)),
            [g, ni] = merge2nodes(g, s(1), st(n), ni);
        else 
            [g, ni] = merge2nodes(g, t(1), st(n), ni);
        end
    end
    
    if any(s == st(u-1)),
        mv(s) = s(1);    % update belonging vector
        for n = u-2:-1:2,
            [g, ni] = merge2nodes(g, t(1), st(n), ni);
        end
        mv(t) = t(1);
    else
        mv(t) = t(1);    % update belonging vector
        for n = u-2:-1:2,
            [g, ni] = merge2nodes(g, s(1), st(n), ni);
        end
        mv(s) = s(1);
    end
    
    pokracuj = 1;
    
    while length(ni)>1 && pokracuj
        % minimum cut phase
        [vs, vt] = mincutphase( g );
        st = sum(g(vt, setdiff(1:length(g),vt)));

        if st<mw  && numel(setdiff([mv(s); mv(t)] ,[])) > numel(setdiff([mv(s); mv(t)],mv(ni(vt)))),           % record
            mw = st;
            end_mv = mv;                % end_mv - the last result, 
            end_mv(mv==mv(ni(vt))) = 0; % one partition is 0
        end
       
        if numel(setdiff( [ mv(s(1)); mv(t(1))], [mv(ni(vt)); mv(ni(vs))])) == 0,
            pokracuj = 0;
        else 
            % update belonging vector
            mv(mv == mv(ni(vt))) = mv(ni(vs));

            % merge t to s
            [g, ni] = merge2nodes(g, vs, vt, ni);
        end
    end
    mv = end_mv;  mv(mv>0) = 1;         % another partiotion is 1
    
% +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+ %

function [ug, uni] = merge2nodes(g, n1, n2, ni)
%MERGE2NODES takes n1 and n2 as nodes in the given graph g and merges n2
% into n1 !
% ni  - node index, records original node mark
% ug  - the updated graph
% uni - updated node index


    % weight changed
    g(n1,:) = g(n1,:) + g(n2,:);
    g(:,n1) = g(n1,:);
    g(n1,n1)= 0;    % self-edge

    selectN = true(1, size(g,2));
    selectN(n2) = false;
    ug  = g(selectN, selectN);
    uni = ni(selectN);
    
% +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+ %

function [vs, vt] = mincutphase( g )
%MINCUTPHASE maximum adjacency search or maximum cardinality search
% yields the s-t-cut.
% vs and vt is the two vertices added last!

    N  = length(g);
    if N == 2,                  % only two nodes 
        vs = 1;     vt = 2;     
        return;     
    end

    s  = randi(N);              % start node by random selection
    cs = s;                     % index of start node in current graph
    ni = 1:N;
    for n = 1:N-3
        [nul, t] = max(g(cs,:));% most tightly connected with cs
        [g, ni]  = merge2nodes(g, cs, t, ni);
        cs = find(ni==s);       % find start node current site
    end

    % ni only has the last two-nodes and start node
    [nul, vs] = max(g(cs,:));   vs = ni(vs);
    vt = setdiff(ni,[vs, ni(cs)]);  
    