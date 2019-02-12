function polygons = mesh_xsections( verts, faces, planes, precision, display )
% GETXSECTIONPOLYGONS finds cross-sections of a mesh by a family of
% planes and returns the resulting polygons
%
% OUTPUT:
%   polygons - a cell array of cells holding one or more x-section
%              polygon matrics for each plane
% INPUT:
%   verts - N x 3 matrix of all the vertices in the mesh
%   faces - M x 3 matrix of all the faces in the mesh
%   planes - a structure defining the x-section planes as a matrix 
%            of K x 3 normals planes.n and K x 3 orgin points plane.r
%   precision - vertices closer than this will be fused
%   display - 0, no plot, 1 - vertices, 2 - edges, 3 - faces
%
% Copyright: Yury Petrov, 2019
%

if nargin < 5 || isempty( display )
    display = 0;
end
if nargin < 4 || isempty( precision )
    precision = 1e-6; % vertices closer than this in space will be fused
end
if nargin < 3 || isempty( planes )
    planes.n = [ 1 0 0; 0 1 0; 0 0 1 ]; % 3 cardinal crossections
    planes.r = 1e-6 + [ 0 0 0; 0 0 0; 0 0 0 ]; % all plane origins to 0
end
nplanes = size( planes.n, 1 );
if nargin < 2 || isempty( faces ) || isempty( verts )
    error( 'Mesh Vertices and Faces input parameters needed!' );
end

% correct the mesh for vertices pointing to the same location (happens a
% lot in real-world meshes for some reason)
cnt = 0;
for i = 1 : size( verts, 1 )
    d = sum( ( verts( i, : ) - verts( i + 1 : end, : ) ).^2, 2 );
    zi = i + find( d < precision^2 );
    zi( zi == i ) = []; % exclude self distance
    if ~isempty( zi ) % replace indices for repeated locations by i
        for j = 1 : length( zi )
            faces( faces == zi( j ) ) = i;
            cnt = cnt + 1;
        end
    end
end
fprintf( 'Corrected the mesh for %d vertex replicas\n', cnt );

% check for duplicated faces
sf =  sort( faces, 2 );
[ ~, ia ] = unique( sf, 'rows' );
dups = setdiff( 1 : size( faces, 1 ), ia );
if ~isempty( dups )
    fprintf( 'Mesh topology: Found duplicate faces, removing %d duplicates!\n', length( dups ) );   
    faces( dups, : ) = [];
end

% plot the input mesh
if display == 1
    plot3( verts( :, 1 ), verts( :, 2 ), verts( :, 3 ), 'k.' )
elseif display == 2
    patch( 'Faces', faces, 'Vertices', verts, 'FaceAlpha', 0, 'EdgeAlpha', 1 ); % display the mesh edges
elseif display == 3
    patch( 'Faces', faces, 'Vertices', verts, 'FaceColor', 'y', 'FaceAlpha', 0.75, 'EdgeAlpha', 0.75 ); % display the mesh
end
hold on;

% find cross-sections
warning( 'off', 'MATLAB:triangulation:PtsNotInTriWarnId' ); % turn off Matlab warning about verts not accounted by faces
TR = triangulation( faces, verts ); % create the Matlab triangulation object for the mesh
E = TR.edges; % get all the edges in the mesh
nedges = size( E, 1 );
nverts = size( verts, 1 );
er = verts( E( :, 1 ), : );
ed = verts( E( :, 2 ), : ) - er; % edge vectors
el = sqrt( sum( ed.^2, 2 ) ); % edge lengths
en = ed ./ repmat( el, 1, 3 ); % normalized edge directions
polygons = cell( 1, nplanes ); % pre-allocate for speed

for s = 1 : nplanes
    % check that the plane doesn't intersect mesh vertices
    pn = repmat( planes.n( s, : ), nverts, 1 );
    pr = repmat( planes.r( s, : ), nverts, 1 );
    d = dot( pn, verts - pr, 2 );
    zinds = find( abs( d ) < precision );
    if ~isempty( zinds )
        fprintf( [ 'Plane crosses vertices [ ' repmat( '%d ', 1, length( zinds ) ) '] within precision %e!\n'], zinds, precision );
        verts( zinds, : ) = verts( zinds, : ) + ...
            repmat( planes.n( s, : ), length( zinds ), 1 ) * 2 * precision; % shift the vertices along the plane normal by 2 * precision
        % recalculate edges
        er = verts( E( :, 1 ), : );
        ed = verts( E( :, 2 ), : ) - er; % edge vectors
        el = sqrt( sum( ed.^2, 2 ) ); % edge lengths
        en = ed ./ repmat( el, 1, 3 ); % normalized edge directions
    end
    
    % distance to the plane along the edge rays
    pn = repmat( planes.n( s, : ), nedges, 1 );
    pr = repmat( planes.r( s, : ), nedges, 1 );
    d = dot( pn, pr - er, 2 ) ./ dot( en, pn, 2 );
    
    % find distances smaller than the edge length
    % inti = d > -precision  & d < el + precision; % logical indices of the edges intersecting the plane within the precision
    inti = d > 0  & d < el; % logical indices of the edges intersecting the plane within the precision
    
    if ~sum( inti ) == 0 % found some intersections
        rinter = er( inti, : ) + repmat( d( inti ), 1, 3 ) .* en( inti, : ); % intersection vectors
        % plot3( rinter( :, 1 ), rinter( :, 2 ), rinter( :, 3 ), 'r.', 'MarkerSize', 13 );
        
        einter = E( inti, : );  % interaction edges
        finter = edgeAttachments( TR, einter ); % faces attached to the interaction edges
        
        % check basic topology
        for edge = 1 : length( finter )
            nfaces = length( finter{ edge } );
            if nfaces > 2
                error( 'Mesh topology: edge %d has more than 2 faces attached to it!', edge );
            elseif nfaces == 0
                error( 'Mesh topology: edge %d has no faces attached to it!', edge );
            elseif nfaces == 1
                fprintf( 'Mesh topology: edge %d has only one face attached to it, folding over\n', edge );
                finter{ edge } = [ finter{ edge } finter{ edge } ]; % fold the face on itself at the mesh edge
                % r = [ verts( einter( edge, 1 ), : ); verts( einter( edge, 2 ), : ) ];
                % plot3( r( :, 1 ),  r( :, 2 ),  r( :, 3 ), 'g-', 'LineWidth', 2 );
                % disp( einter( edge, : ) );
            end
        end
        
        % check if every face appears twice now
        finter = cell2mat( finter );
        [ uf, ~, ic ] = unique( finter(:) );
        counts = accumarray( ic, 1 );
        de = find( counts == 1 );
        if ~isempty( de )
            fprintf( 'Singleton faces found!' );
            for i = 1 : length( de )
                sf = uf( de( i ) ); % the singleton face
                fprintf( '\nSingleton face: %d\n', sf );
                fv = faces( sf, : ); % vertices of this face
                fe(1) = find( E( :, 1 ) == fv(1) & E( :, 2 ) == fv(2) | E( :, 1 ) == fv(2) & E( :, 2 ) == fv(1) );
                fe(2) = find( E( :, 1 ) == fv(2) & E( :, 2 ) == fv(3) | E( :, 1 ) == fv(3) & E( :, 2 ) == fv(2) );
                fe(3) = find( E( :, 1 ) == fv(1) & E( :, 2 ) == fv(3) | E( :, 1 ) == fv(3) & E( :, 2 ) == fv(1) );
                fprintf( 'Edge numbers:\t' );
                fprintf( '%d ', fe );
                fprintf( '\nEdge lengths:\t' );
                fprintf( '%.6f ', el( fe ) );
                fprintf( '\nEdge interxs:\t' );
                fprintf( '%.6f ', d( fe ) );
                r = [ verts( fv, : ); verts( fv(1), : ) ];
                plot3( r( :, 1 ),  r( :, 2 ),  r( :, 3 ), 'g-', 'LineWidth', 2 );
            end
            fprintf( '\n ' );
        end
        
        % form the loops
        loops = getLoops( finter );
    else
        loops = {};
    end
    
    % Store the results as vertex polygons rather than index loops.
    if ~isempty( loops )
        polys = cell( size( loops ) );
        for p = 1 : length( loops )
            if ~isempty( loops{ p } )
                polys{ p } = rinter( loops{ p }, : );
            end
            polygons{ s } = polys;
        end
    end
end

        
function loops = getLoops( inds )
%
% Find closed loops in Nx2 array of interger links. Returns arrays of 
% indices to inds forming the loops. If any row of inds has the same 
% integer in both columns, this is treated as a cut in the loop, 
% and the loop is closed to the next cut found in it.
%
% Copyright: Yury Petrov, 2019
%

loops = cell( 1, ceil( size( inds, 1 ) / 3 ) ); % preallocate the max number of space loops possible
looped = zeros( size( inds, 1 ), 1 ); % marks indices as looped already
newloop = true;
cnt = 1; % loop counter
i = 1;   % loop index counter
cut = 0; % cut location in the loop, 0 means no cut
iv = inds( i, 1 );
if iv == inds( i, 2 ) % if happened to start the loop with a cut
    cut = 1;
end
while ~all( looped )
    if ~looped( i )
        looped( i ) = 1;  % mark as being in a loop already
        if newloop
            loops{ cnt } = i; % start a new loop
            newloop = false;
        end
        sm = sum( inds == iv, 2 );
        ni = find( sm ); % at least one of the two indices matches
        ni = ni( ni ~= i ); % skip the self link
        nv = inds( ni, inds( ni, : ) ~= iv ); % index value at the new location
        if length( nv ) == 1
            if ni ~= loops{ cnt }(1) % if a new index is not the beginning of the loop
                loops{ cnt } = [ loops{ cnt } ni ]; % append the new index to the loop
            else % close the loop
                if length( loops{ cnt } ) > 1  % loops should be at least 2 elements long
                    cnt = cnt + 1;
                end
                newloop = true;
                ni = find( looped == 0, 1 ); % find the next unexplored link
                nv = inds( ni, 1 );
                if nv == inds( i, 2 ) % happens to be a cut
                    cut = ni;
                end
            end
        elseif length( nv ) > 1
            error( 'Link topology: index %d connected to more than 2 indices!\n', ni );
        else % no connection
            if sm( ni ) == 2 % link with both elements equal (loops on itself), cut loop
                looped( ni ) = 1; % mark as looped
                loops{ cnt } = [ loops{ cnt } ni ]; % append the new index to the loop
                if cut == 0 % first time a cut was encountered
                    cut = length( loops{ cnt } ); % store the cut location in the loop
                    ns = loops{ cnt }(1); % start from the beginning of the loop again
                    nv = inds( ns, 2 ); % but from the second link this time
                    ni = find( sum( inds == nv, 2 ) ); % at least one of the two indices matches
                    ni = ni( ni ~= ns ); % skip the self link
                    nv = inds( ni, 1 );
                else % the other edge found, close the cut loop
                    if length( loops{ cnt } ) > 1  % loops should be at least 2 elements long
                        if cut ~= 1 % started the loop from the cut edge
                            tmp = loops{ cnt }( cut + 1 : end ); % the second part of the loop
                            loops{ cnt } = [ fliplr( tmp ) loops{ cnt }( 1 : cut ) ]; % flip and prepend to the first part
                        end
                        cnt = cnt + 1;
                    end
                    newloop = true;
                    ni = find( looped == 0, 1 ); % find the next unexplored link
                    nv = inds( ni, 1 );
                    if nv == inds( i, 2 ) % happens to be a cut
                        cut = 1;
                    else
                        cut = 0;
                    end
                end
            else % a face is only encountered once along a loop, topolology defect!
                error( 'Link topology: index %d connected to a singleton index!\n', i );
            end
        end
        i = ni; % interate to the new index
        iv = nv;
    else
        error( 'Link topology: there is some problem with links!' );
    end
end
loops = loops( ~cellfun( @isempty, loops ) ); % remove empty cells
