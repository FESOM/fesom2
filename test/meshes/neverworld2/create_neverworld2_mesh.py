#!/usr/bin/env python3
"""
Python translation of the MATLAB NeverWorld mesh generator.
Indices remain 1-based to match Fortran/FESOM conventions.
Plotting code is removed as requested.
"""

import numpy as np
from scipy.interpolate import CubicSpline
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.tri as mtri

def main():
    print(" Create neverworld2 mesh.")
    #___________________________________________________________________________
    # Mesh parameters
    outputflag = 1
    Lx         = 60.0         # degrees
    Ly         = 70.0         # domain from -Ly to Ly
    Lyperiodic = [-60.0, -40.0]   # Latitude range for periodic boundary
    
    
    cyclic_len = Lx

    meshtype   = 1            # 1 = lon/lat grid
    dxm        = 1.0          # triangle side at equator / Resolution
    zonal      = 0            # 0 = castellated in x, 1 = castellated in y
    
    # define layer thicknesses of each layer
    hthick = np.array([25, 50, 100, 125, 150, 175, 200, 225, 250, 300, 350, 400, 500, 550, 600])        
    max_depth  = -hthick.sum()  #-4000.0
    
    
    do_plot    = True

    #___________________________________________________________________________
    # Grid generation
    xx = np.arange(0, Lx + dxm, dxm)
    
    if meshtype == 1:
        yy = np.arange(-Ly, Ly + dxm, dxm)
        
    else:
        # Mercator version
        ylist = []
        y = 0
        while y < Ly:
            y = y + dxm * np.cos(np.pi * (y + dxm/2) / 180)
            ylist.append(y)
        ylist = np.array(ylist)
        ylist = ylist * Ly / ylist[-1]
        
        # make ycoordinates sym,etric around equator
        yy = np.concatenate((-ylist[::-1], [0], ylist))
        
    nx = len(xx)
    ny = len(yy)

    #___________________________________________________________________________
    # Create the coordinate arrays
    # create a regular 2D array of points, but shift their x (or y) coordinates
    # in every second row.
    xcoord, ycoord = np.meshgrid(xx, yy)

    # Castellated / shifted grid adjustment
    if zonal == 1:
        # shift x for even rows
        xcoord[1::2, :] += 0.5 * dxm
    else:
        # shift y for alternating columns
        if meshtype == 1:
            ycoord[:, 1::2] += 0.5 * dxm
        else:
            ycoord[:, 1::2] += 0.5 * dxm * np.cos(np.pi * ycoord[:, 1::2] / 180)


    #___________________________________________________________________________
    # fortran based node indices (1-based) --> here shape [ny, nx]
    nodnum = np.arange(0, nx * ny).reshape(ny, nx)

    # reshape x/y coordinates to 1D arrays
    xcoord = xcoord.ravel()
    ycoord = ycoord.ravel()
    
    # nodind (vertical wall index)
    nodind = np.zeros(ny*nx)
    
    #___________________________________________________________________________
    # Build triangulation elem array [3, n2de]
    aux_tri = []
    #  castellated in y
    if zonal == 1:
        for n in range(nx - 1):
            for nn in range(0, ny - 1, 2):
                # two triangles in "odd" row
                aux_tri.append([nodnum[nn, n], nodnum[nn+1, n], nodnum[nn, n+1]])
                aux_tri.append([nodnum[nn+1, n], nodnum[nn+1, n+1], nodnum[nn, n+1]])

            for nn in range(1, ny - 1, 2):
                # two triangles in "even" row
                aux_tri.append([nodnum[nn, n], nodnum[nn+1, n], nodnum[nn+1, n+1]])
                aux_tri.append([nodnum[nn, n], nodnum[nn+1, n+1], nodnum[nn, n+1]])
                
    # castellated in x
    else:
        # first set: odd columns
        for n in range(0, nx - 1, 2):
            for nn in range(ny - 1):
                aux_tri.append([nodnum[nn, n], nodnum[nn+1, n], nodnum[nn, n+1]])
                aux_tri.append([nodnum[nn+1, n], nodnum[nn+1, n+1], nodnum[nn, n+1]])

        # second set: even columns
        for n in range(1, nx - 1, 2):
            for nn in range(ny - 1):
                aux_tri.append([nodnum[nn, n], nodnum[nn+1, n+1], nodnum[nn, n+1]])
                aux_tri.append([nodnum[nn, n], nodnum[nn+1, n], nodnum[nn+1, n+1]])
    
    # mesh is now given by xcoord, ycoord, and tri.
    tri = np.array(aux_tri, dtype=int).T
    n2de = tri.shape[1]
    del(aux_tri)

    #___________________________________________________________________________
    # Cyclic reduction (topological merge)
    # Make mesh cyclic between 40 and 60 S
    # Cyclic reduction:
    # The respective nodes in xcoord, ycoord at 0 and 60 E are 
    # equivalent, and we eliminate the eastern.
    idx_left  = np.where((ycoord >= Lyperiodic[0]) & (ycoord <= Lyperiodic[1]) & (xcoord == 0 ))[0]
    idx_right = np.where((ycoord >= Lyperiodic[0]) & (ycoord <= Lyperiodic[1]) & (xcoord == Lx))[0]
    if len(idx_left) != len(idx_right): raise RuntimeError("Cyclic merge mismatch.")

    # Replace references in tri (still 1-based!), replace right node index with 
    # left one in elem list --> close periodicity
    for L, R in zip(idx_left, idx_right): tri[tri == R] = L

    # Identify used vertices
    keep     = np.unique(tri.flatten())
    
    # Create mapping dictionary between old and new vertex indices
    mapping  = {old_index: new_index for new_index, old_index in enumerate(keep)}

    # compute right boundary reduced triangulation
    tri      = np.vectorize(mapping.get)(tri)
    xcoord   = xcoord[keep]
    ycoord   = ycoord[keep]
    nodind   = nodind[keep]
    n2dn     = len(xcoord)
    n2de     = tri.shape[1]
    
    xc       = xcoord[tri]
    isnotpbnd= (xc.max(axis=0)-xc.min(axis=0)) <= cyclic_len/2.0
    print(' --> number of vertices: ', n2dn )
    print(' --> number of elements: ', n2de )
    # --> mesh definiiton is finished until here     
    
    #___________________________________________________________________________
    # Topography generation   
    depth = max_depth * np.ones_like(xcoord)

    
    #___________________________________________________________________________
    # Middle topographic ridge profile @30°: (cubic spline)
    # S = ... 
    # row1:	polynomial value at x=10
    # row2:	polynomial value at x=30
    # row3:	derivative at x=10
    # row4:	derivative at x=30
    # The cubic polynomial has the form:
    #       P(x)=a0+a1*x+a2*x^2+a3*x^3

    S = np.array([[1, 10, 100, 1000],
                  [1, 30, 900, 27000],
                  [0, 1, 20, 300],
                  [0, 1, 60, 2700]])
    
    # We want to enforce:
    # P(10) = 0
    # P(30) = 2000
    # P'(10) = 0 (flat slope at start)
    # P'(30) = 0 (flat slope at peak)
    Sbnd = np.array([0, 2000, 0, 0])
    
    aa = np.linalg.solve(S, np.array([0, 2000, 0, 0]))
    for i in range(n2dn):
        di = xcoord[i]
        if di > 30:
            di = 30 - (di - 30)
        if di >= 10:
            f = aa[0] + aa[1]*di + aa[2]*di**2 + aa[3]*di**3
        else:
            f = 0
        depth[i] += f
        
    #___________________________________________________________________________
    # do topographic shelf
    S = np.array([[1.0, 2.5, 2.5**2, 2.5**3],
                  [1.0, 5.0, 25.0  , 125],
                  [0.0, 1.0, 5.0   , 3.0*2.5**2],
                  [0.0, 1.0, 10.0  , 75.0]])
    aa = np.linalg.solve(S, np.array([3800.0, 0.0, 0.0, 0.0]))

    for i in range(n2dn):
        f = 0
        x = xcoord[i]
        y = ycoord[i]

        if not (y >= Lyperiodic[0]) & (y <= Lyperiodic[1]):
            di = min(abs(x), abs(x - 60), abs(y - 70), abs(y + 70))
        else:
            di = min(np.sqrt(x**2 + (y + 60)**2),
                     np.sqrt(x**2 + (y + 40)**2),
                     np.sqrt((x - 60)**2 + (y + 60)**2),
                     np.sqrt((x - 60)**2 + (y + 40)**2))
        
        if di <= 2.5:
            f = 3800
        elif di <= 5:
            f = aa[0] + aa[1]*di + aa[2]*di**2 + aa[3]*di**3
        
        depth[i] = max_depth + max(depth[i] - max_depth, f)
        
    #___________________________________________________________________________
    # semicircle
    for i in range(n2dn):
        di = np.sqrt(xcoord[i]**2 + (ycoord[i] + 50)**2)
        f = 2000 if abs(di - 10) <= 1.5 else 0
        depth[i] = max_depth + max(depth[i] - max_depth, f)

    #___________________________________________________________________________
    # triangle depth
    depth_elem = np.min(depth[tri], axis=0)  # tri is 1-based
    
    
    #___________________________________________________________________________
    #
    if do_plot: 
        fig, axes = plt.subplots(1, 2, figsize=(14, 7))
        triang_A = mtri.Triangulation(xcoord, ycoord, tri[:,isnotpbnd].T)
        h1 = axes[0].tripcolor(triang_A, depth, linewidth=0.01, color="black", vmin=max_depth, vmax=0,
                               cmap=matplotlib.colormaps["GnBu"].resampled(32))
        axes[0].set_title("Mesh (@vertices)")
        axes[0].set_aspect("equal", "box")
        axes[0].set_xlabel("Longitude")
        axes[0].set_ylabel("Latitude")
        
        triang_B = mtri.Triangulation(xcoord, ycoord, tri[:,isnotpbnd].T)
        h2 = axes[1].tripcolor(triang_B, depth_elem[isnotpbnd], linewidth=0.01, color="black", vmin=max_depth, vmax=0,
                               cmap=matplotlib.colormaps["GnBu"].resampled(32))
        axes[1].set_title("Mesh (@elements)")
        axes[1].set_aspect("equal", "box")
        axes[1].set_xlabel("Longitude (wrapped)")
        
        fig.colorbar(h2, )

        plt.tight_layout()
        plt.savefig("mesh_neverworld2.png", dpi=150)
        
        
        from mpl_toolkits.mplot3d import Axes3D    
        fig = plt.figure(figsize=(12, 8))
        ax  = fig.add_subplot(111, projection='3d')
        
        # Important: depth should be negative (bathymetry)
        surf = ax.plot_trisurf(triang_A, depth, cmap=matplotlib.colormaps["GnBu"].resampled(32), linewidth=0.1,
                                antialiased=True, shade=True)

        fig.colorbar(surf, ax=ax, shrink=0.5, label="Depth (m)")
        ax.set_xlabel("Longitude")
        ax.set_ylabel("Latitude")
        ax.set_zlabel("Depth (m)")
        ax.set_title("NeverWorld 3D Bathymetry")

        # Make the vertical exaggeration smaller
        ax.set_box_aspect((1, 1, 0.2))   # X:Y:Z = 1:1:0.2 → compressed Z
        ax.view_init(elev=60, azim=-90)
        plt.tight_layout()
        plt.savefig("mesh_neverworld2_3d.png", dpi=150)
    

    
    #___________________________________________________________________________
    # Wind stress (Cubic spline)
    yi     = np.array([-70, -45, -15, 0, 15, 45, 70])
    taui   = np.array([0., 0.2, -0.1, -0.02, -0.1, 0.1, 0])
    yw     = np.arange(-70, 71)
    spline = CubicSpline(yi, taui)
    tauw   = spline(yw)
    tau_node = spline(ycoord)
    if do_plot: 
        fig = plt.figure(figsize=(7, 14))
        ax = plt.gca()
        plt.plot(tauw, yw)
        ax.set_xlabel(" Wind-stress")
        ax.set_ylabel("Latitude")
        ax.set_title("NeverWorld2 wind-stress profile ")
        
        # Make the vertical exaggeration smaller
        plt.tight_layout()
        plt.savefig("windstress_neverworld2.png", dpi=150)
    
        
    #___________________________________________________________________________
    # Output files (Fortran-style)    
    if outputflag == 1:
        #_______________________________________________________________________
        # write nodal information
        print(' --> Write file:')
        print('    - nod2d.out')
        with open("nod2d.out", "w") as f:
            f.write(f"{n2dn}\n")
            for i in range(n2dn):
                f.write(f"{i+1:8d} {xcoord[i]:8.4f} {ycoord[i]:8.4f} {int(nodind[i]):8d}\n")

        #_______________________________________________________________________
        # write elemental information 
        print('    - elem2d.out')
        tri = tri + 1 # (fortran style indices)
        with open("elem2d.out", "w") as f:
            f.write(f"{n2de}\n")
            for e in range(n2de):
                f.write(f"{tri[0,e]:8d} {tri[1,e]:8d} {tri[2,e]:8d}\n")

        #_______________________________________________________________________
        # write vertical information 
        zbar   = np.zeros(hthick.size+1)
        for i in range(len(hthick)):
            zbar[i+1] = zbar[i] + hthick[i]
        
        print('    - aux3d.out')        
        with open("aux3d.out", "w") as f:
            f.write(f"{len(zbar)}\n")
            for z in zbar:
                f.write(f"{z}\n")        
        
        print('    - depth_zlev.out (alternativ)')
        with open("depth_zlev.out", "w") as f:
            f.write(f"{len(zbar)}\n")
            for z in zbar:
                f.write(f"{z}\n")
                
        print('    - depth@node.out (alternativ)')        
        with open("depth@node.out", "w") as f:
            for d in depth:
                f.write(f"{d:7.1f}\n")
                
        print('    - depth@elem.out (alternativ)')        
        with open("depth@elem.out", "w") as f:        
            for d in depth_elem:
                f.write(f"{d:7.1f}\n")
            for d in depth:
                f.write(f"{d:7.1f}\n")
                
        #_______________________________________________________________________        
        # write windstress information 
        print('    - windstress@node.out')        
        with open("windstress@node.out", "w") as f:
            for tw in tau_node:
                f.write(f"{tw}\n")

    print(" Mesh generation completed.")


if __name__ == "__main__":
    main()
