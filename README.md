## mesh_xsections
<img width="540" alt="screen shot 2019-02-07 at 7 13 20 pm" src="https://user-images.githubusercontent.com/46988982/52511437-d6d67f80-2bb4-11e9-970b-a81049337552.png"><img width="740" alt="screen shot 2019-02-08 at 2 53 51 pm" src="https://user-images.githubusercontent.com/46988982/52511637-bbb83f80-2bb5-11e9-953d-f5a14f1c1d9f.png">

mesh_xsections() is a Matlab function returning cross-sections of a triangulation mesh with a set of planes. Each cross-section is a cell array with one or more polygons (consecutive points in 3D). The function handles meshes with duplicated vertices (by fusing them) and open boundary edges (by cutting across the edge).
