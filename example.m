%
% Make cross-sections of 3D objects defined in the Wavefront .obj file
%

% read the data file (Wavefront .obj file saved into a matlab structure, e.g., using read_wobj)
load( 'stomach+duodenum+pancreas.mat', 'O' );
% load( '12140_Skull_v3_L2.mat', 'O' );

% count the number of objects in the file
nobj = 0;
for i = 1 : length( O.objects )
    if strcmp( O.objects( i ).type, 'o' ) || strcmp( O.objects( i ).type, 'g' )
        nobj = nobj + 1; 
    end
end
fprintf( 'Found %d objects in the input file\n', nobj );

verts = O.vertices;
fprintf( 'Found %d vertices in the input file\n', size( verts, 1 ) );

% define cross-section planes
planes.n = [ 1 0 0; 0 1 0; 0 0 1 ]; % 3 cardinal crossections
planes.r = [ 0 0 2; 0 0 2; 0 0 2 ]; % plane origins at [ 0, 0, 2 ]

figure( 'Position', [ 100, 100, 800, 800 ] ), hold on;
axis vis3d equal;
view( -45, 45 );
xlabel( 'X', 'FontSize', 14 );
ylabel( 'Y', 'FontSize', 14 );
zlabel( 'Z', 'FontSize', 14 );

for i = 1 : length( O.objects )
    if strcmp( O.objects( i ).type, 'o' ) || strcmp( O.objects( i ).type, 'g' ) % start of the new object or group
        objectName = O.objects( i ).data;
        disp( [ 'Processing object ' objectName ] ); % print out object name
    elseif strcmp( O.objects( i ).type, 'f' ) % object faces block
        faces = O.objects( i ).data.vertices;
        patch( 'Faces', faces, 'Vertices', verts, 'FaceAlpha', 0, 'EdgeAlpha', 1 ); % display the mesh
        
        % get intersection polygons for all planes
        polygons = mesh_xsections( verts, faces, planes, [], 2 );
        
        % draw the object polygons
        for s = 1 : numel( polygons )
            if ~isempty( polygons{ s } )
                for p = 1 : numel( polygons{ s } )
                    patch( polygons{ s }{ p }( :, 1 ), ...
                           polygons{ s }{ p }( :, 2 ), ...
                           polygons{ s }{ p }( :, 3 ), 'y', 'FaceAlpha', 0.5 );
                    plot3( polygons{ s }{ p }( :, 1 ), ...
                           polygons{ s }{ p }( :, 2 ), ...
                           polygons{ s }{ p }( :, 3 ), 'r*' );
                end
            end
        end
    end
end
