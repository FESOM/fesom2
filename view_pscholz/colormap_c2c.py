# Patrick Scholz, 14.12.2017
def colormap_c2c(cmin,cmax,cref,cnumb,cname,cstep=[]):
    
    import numpy as np
    from matplotlib.colors import LinearSegmentedColormap
    ##import cmocean
    # cmin ... value of minimum color
    # cmax ... value of maximum color
    # cref ... value of reference color
    # cnumb... minimum number of colors that should be in colormap
    # cname... name of colormap ---> see colormap definition
    
    # colormaps are:
    # --> blue2red
    # --> red2blue
    # --> grads
    # --> rainbow
    # --> heat
    # --> jet
    # --> jet_w
    # --> hsv
    # --> gnuplot
    # --> arc
    # --> wbgyr
    # --> odv
    # --> odv_w
    
    # get colormap matrix from predifined cmaps
    # --> import matplotlib.pyplot as plt
    # --> cmap = plt.get_cmap('rainbow',9)
    # --> cmaplist = [cmap(i) for i in range(cmap.N)]
    
    #___________________________________________________________________________
    # calculate optimal color step size
    if not cstep:
        cdelta   = (cmax-cmin)/cnumb
        cstep_all= np.array([0.1, 0.2, 0.25, 0.5, 1.0, 2.0, 2.5, 5.0, 10.0, 20.0, 25.0, 50.0])
        cstep_all=cstep_all*(10**np.floor(np.log10( np.abs(cdelta) )))
        cstep_i  = np.squeeze(np.where(cstep_all<=cdelta))
        cstep_i  = cstep_i[-1] 
        cstep    = cstep_all[cstep_i]
        
    #print('[cmin,cmax,cref]=',cmin,cmax,cref)
    #print('cstep  = ',cstep)
    #print('cdelta = ',cdelta)
    
    #___________________________________________________________________________
    # calculate colormap levels
    #print(np.arange(cref-cstep,cmin-cstep,-cstep))
    
    clevel   = np.concatenate((np.sort(np.arange(cref-cstep,cmin-cstep,-cstep)),np.arange(cref,cmax+cstep,cstep)))
    #print(clevel)
    if np.abs(clevel.min())>1.0e-15:
        #clevel   = np.around(clevel, -np.int32(np.floor(np.log10(np.abs( clevel.min() ))-2) ) )
        clevel   = np.around(clevel, -np.int32(np.floor(np.log10(np.abs( cstep ))-2) ) )
    
    clevel   = np.unique(clevel)
    #print(clevel)
    #print(clevel[:-1]-clevel[1:])
    #print(clevel)
    cdelta2  = clevel[-1]-clevel[0]
    if cmin==0.0 and clevel[0]<cmin:
        #clevel = clevel[1:]
        clevel[0]=0.0
        clevel   = np.unique(clevel)
        cdelta2  = clevel[-1]-clevel[0]
    
        
    #___________________________________________________________________________
    # different colormap definitions
    cmap_def = []
    #-->________________________________________________________________________
    if cname=='blue2red':
        cmap_def = [(0.0                                            , [0.0, 0.19, 1.0]),
                    (((cref-clevel[0])/2)/cdelta2                    , [0.0, 0.72, 1.0]),
                    ((cref-clevel[0])/cdelta2                        , [1.0, 1.0 , 1.0]),
                    ((cref-clevel[0]+(clevel[-1]-cref)/2)/cdelta2    , [1.0, 0.6 , 0.0]),
                    ((clevel[-1]-clevel[0])/cdelta2                    , [1.0, 0.19, 0.0] )]
    if cname=='red2blue':
        cmap_def = [(0.0                                            , [1.0, 0.19, 0.0]),
                    (((cref-clevel[0])/2)/cdelta2                    , [1.0, 0.6 , 0.0]),
                    ((cref-clevel[0])/cdelta2                        , [1.0, 1.0 , 1.0]),
                    ((cref-clevel[0]+(clevel[-1]-cref)/2)/cdelta2    , [0.0, 0.72, 1.0]),
                    ((clevel[-1]-clevel[0])/cdelta2                    , [0.0, 0.19, 1.0])]    
    #-->________________________________________________________________________
    elif cname=='green2orange':    
        cmap_def = [(0.0                                            , [0.2196,    0.4196,      0.0]),
                    (((cref-clevel[0])*0.3333)/cdelta2                , [0.6039,    0.8039,      0.0]),
                    (((cref-clevel[0])*0.6666)/cdelta2                , [0.8000,    1.0000,      0.0]),
                    ((cref-clevel[0])/cdelta2                        , [1.0000,    1.0000,   1.0000]),
                    ((cref-clevel[0]+(clevel[-1]-cref)*0.5)/cdelta2 , [1.0000,    0.6000,      0.0]),
                    ((clevel[-1]-clevel[0])/cdelta2                    , [0.6000,    0.2000,      0.0])]
    elif cname=='orange2green':    
        cmap_def = [(0.0                                            , [0.6000,    0.2000,      0.0]),
                    (((cref-clevel[0])*0.5)/cdelta2                    , [1.0000,    0.6000,      0.0]),
                    ((cref-clevel[0])/cdelta2                        , [1.0000,    1.0000,   1.0000]),
                    ((cref-clevel[0]+(clevel[-1]-cref)*0.3333)/cdelta2,[0.8000,    1.0000,      0.0]),
                    ((cref-clevel[0]+(clevel[-1]-cref)*0.6666)/cdelta2,[0.6039,    0.8039,      0.0]),
                    ((clevel[-1]-clevel[0])/cdelta2                    , [0.2196,    0.4196,      0.0])]
        
    #-->________________________________________________________________________
    elif cname=='grads':    
        cmap_def = [(0.0                                            , [0.6275, 0.0   , 0.7843]),
                    (((cref-clevel[0])*0.25)/cdelta2                , [0.1176, 0.2353, 1.0000]),
                    (((cref-clevel[0])*0.50)/cdelta2                , [0.0   , 0.6275, 1.0000]),
                    (((cref-clevel[0])*0.75)/cdelta2                , [0.0   , 0.8627, 0.0   ]),
                    ((cref-clevel[0])/cdelta2                        , [1.0   , 1.0   , 1.0   ]),
                    ((cref-clevel[0]+(clevel[-1]-cref)*0.25)/cdelta2, [0.9020, 0.8627, 0.1961]),
                    ((cref-clevel[0]+(clevel[-1]-cref)*0.50)/cdelta2, [0.9412, 0.5098, 0.1569]),
                    ((cref-clevel[0]+(clevel[-1]-cref)*0.75)/cdelta2, [0.9804, 0.2353, 0.2353]),
                    ((clevel[-1]-clevel[0])/cdelta2                    , [0.9412, 0.0   , 0.5098])]
    #-->________________________________________________________________________
    elif cname=='rainbow':    
        cmap_def = [(0                                                , [0.5   , 0.0   , 1.0   ]),
                    (((cref-clevel[0])*0.25)/cdelta2                , [0.25  , 0.3826, 0.9807]),
                    (((cref-clevel[0])*0.50)/cdelta2                , [0.0   , 0.7071, 0.9238]),
                    (((cref-clevel[0])*0.75)/cdelta2                , [0.25  , 0.9238, 0.8314]),
                    ((cref-clevel[0])/cdelta2                        , [0.5   , 1.0   , 0.7071]),
                    ((cref-clevel[0]+(clevel[-1]-cref)*0.25)/cdelta2, [0.75  , 0.9238, 0.5555]),
                    ((cref-clevel[0]+(clevel[-1]-cref)*0.50)/cdelta2, [1.0   , 0.7071, 0.3826]),
                    ((cref-clevel[0]+(clevel[-1]-cref)*0.75)/cdelta2, [1.0   , 0.3826, 0.1950]),
                    ((clevel[-1]-clevel[0])/cdelta2                    , [1.0   , 0.0   , 0.0   ])]
    #-->________________________________________________________________________
    elif cname=='heat':    
        cmap_def = [(0.0                                            , [1.0   , 1.0   , 1.0]),
                    (((cref-clevel[0])*0.25)/cdelta2                , [1.0   , 0.75  , 0.5]),
                    (((cref-clevel[0])*0.50)/cdelta2                , [1.0   , 0.5   , 0.0]),
                    (((cref-clevel[0])*0.75)/cdelta2                , [0.9375, 0.25  , 0.0]),
                    ((cref-clevel[0])/cdelta2                        , [0.75  , 0.0   , 0.0]),
                    ((cref-clevel[0]+(clevel[-1]-cref)*0.25)/cdelta2, [0.5625, 0.0   , 0.0]),
                    ((cref-clevel[0]+(clevel[-1]-cref)*0.50)/cdelta2, [0.375 , 0.0   , 0.0]),
                    ((cref-clevel[0]+(clevel[-1]-cref)*0.75)/cdelta2, [0.1875, 0.0   , 0.0]),
                    ((clevel[-1]-clevel[0])/cdelta2                    , [0.0   , 0.0   , 0.0])]
    #-->________________________________________________________________________
    elif cname=='jet':    
        cmap_def = [(0.0                                            , [0.0   , 0.0   , 0.5]),
                    (((cref-clevel[0])*0.25)/cdelta2                , [0.0   , 0.0   , 1.0]),
                    (((cref-clevel[0])*0.50)/cdelta2                , [0.0   , 0.5   , 1.0]),
                    (((cref-clevel[0])*0.75)/cdelta2                , [0.0806, 1.0   , 0.8870]),
                    ((cref-clevel[0])/cdelta2                        , [0.4838, 1.0   , 0.4838]),
                    ((cref-clevel[0]+(clevel[-1]-cref)*0.25)/cdelta2, [0.8870, 1.0   , 0.0806]),
                    ((cref-clevel[0]+(clevel[-1]-cref)*0.50)/cdelta2, [1.0   , 0.5925, 0.0]),
                    ((cref-clevel[0]+(clevel[-1]-cref)*0.75)/cdelta2, [1.0   , 0.1296, 0.0]),
                    ((clevel[-1]-clevel[0])/cdelta2                    , [0.5   , 0.0   , 0.0])]
    #-->________________________________________________________________________
    elif cname=='jet_w':    
        cmap_def = [(0.0                                            , [0.0   , 0.0   , 0.5]),
                    (((cref-clevel[0])*0.25)/cdelta2                , [0.0   , 0.0   , 1.0]),
                    (((cref-clevel[0])*0.50)/cdelta2                , [0.0   , 0.5   , 1.0]),
                    (((cref-clevel[0])*0.75)/cdelta2                , [0.0806, 1.0   , 0.8870]),
                    ((cref-clevel[0])/cdelta2                        , [1.0   , 1.0   , 1.0]),
                    ((cref-clevel[0]+(clevel[-1]-cref)*0.25)/cdelta2, [0.8870, 1.0   , 0.0806]),
                    ((cref-clevel[0]+(clevel[-1]-cref)*0.50)/cdelta2, [1.0   , 0.5925, 0.0]),
                    ((cref-clevel[0]+(clevel[-1]-cref)*0.75)/cdelta2, [1.0   , 0.1296, 0.0]),
                    ((clevel[-1]-clevel[0])/cdelta2                    , [0.5   , 0.0   , 0.0])]
    #-->________________________________________________________________________
    elif cname=='hsv':    
        cmap_def = [(0.0                                            , [1.0   , 0.0   , 0.0]),
                    (((cref-clevel[0])*0.25)/cdelta2                , [1.0   , 0.7382, 0.0]),
                    (((cref-clevel[0])*0.50)/cdelta2                , [0.5236, 1.0   , 0.0]),
                    (((cref-clevel[0])*0.75)/cdelta2                , [0.0   , 1.0   , 0.2148]),
                    ((cref-clevel[0])/cdelta2                        , [0.0   , 1.0   , 0.9531]),
                    ((cref-clevel[0]+(clevel[-1]-cref)*0.25)/cdelta2, [0.0   , 0.3085, 1.0]),
                    ((cref-clevel[0]+(clevel[-1]-cref)*0.50)/cdelta2, [0.4291, 0.0   , 1.0]),
                    ((cref-clevel[0]+(clevel[-1]-cref)*0.75)/cdelta2, [1.0   , 0.0   , 0.8320]),
                    ((clevel[-1]-clevel[0])/cdelta2                    , [1.0   , 0.0   , 0.0937])]
    #-->________________________________________________________________________
    elif cname=='gnuplot':    
        cmap_def = [(0.0                                            , [1.0   , 1.0   , 1.0]),
                    (((cref-clevel[0])*0.25)/cdelta2                , [1.0   , 1.0   , 0.0]),
                    (((cref-clevel[0])*0.50)/cdelta2                , [0.9354, 0.6699, 0.0]),
                    (((cref-clevel[0])*0.75)/cdelta2                , [0.8660, 0.4218, 0.0]),
                    ((cref-clevel[0])/cdelta2                        , [0.7905, 0.2441, 0.0]),
                    ((cref-clevel[0]+(clevel[-1]-cref)*0.25)/cdelta2, [0.7071, 0.125 , 0.0]),
                    ((cref-clevel[0]+(clevel[-1]-cref)*0.50)/cdelta2, [0.6123, 0.0527, 0.7071]),
                    ((cref-clevel[0]+(clevel[-1]-cref)*0.75)/cdelta2, [0.5   , 0.0156, 1.0]),
                    ((clevel[-1]-clevel[0])/cdelta2                    , [0.3535, 0.0019, 0.7071])]
    #-->________________________________________________________________________
    elif cname=='arc':    
        cmap_def = [(0.0                                            , [1.0000,    1.0000,    1.0000]),
                    (((cref-clevel[0])*0.25)/cdelta2                , [0.6035,    0.8614,    0.7691]),
                    (((cref-clevel[0])*0.50)/cdelta2                , [0.2462,    0.7346,    0.4610]),
                    (((cref-clevel[0])*0.75)/cdelta2                , [0.2980,    0.7399,    0.2196]),
                    ((cref-clevel[0])/cdelta2                        , [0.7569,    0.8776,    0.0754]),
                    ((cref-clevel[0]+(clevel[-1]-cref)*0.25)/cdelta2, [0.9991,    0.9390,    0.0017]),
                    ((cref-clevel[0]+(clevel[-1]-cref)*0.50)/cdelta2, [0.9830,    0.7386,    0.0353]),
                    ((cref-clevel[0]+(clevel[-1]-cref)*0.75)/cdelta2, [0.9451,    0.2963,    0.1098]),
                    ((clevel[-1]-clevel[0])/cdelta2                    , [0.9603,    0.4562,    0.5268])]
    #-->________________________________________________________________________
    elif cname=='wbgyr':    
        cmap_def = [(0.0                                            , [1.0000,    1.0000,    1.0000]),
                    (((cref-clevel[0])*0.5)/cdelta2                    , [0.2000,    0.6000,    1.0000]),
                    ((cref-clevel[0])/cdelta2                        , [0.0   ,    1.0000,    0.6000]),
                    ((cref-clevel[0]+(clevel[-1]-cref)*0.5)/cdelta2    , [1.0000,    1.0000,       0.0]),
                    ((clevel[-1]-clevel[0])/cdelta2                    , [1.0000,       0.0,       0.0])]        
    #-->________________________________________________________________________
    elif cname=='rygbw':    
        cmap_def = [(0.0                                            , [1.0000,       0.0,       0.0]),
                    (((cref-clevel[0])*0.5)/cdelta2                    , [1.0000,    1.0000,       0.0]),
                    ((cref-clevel[0])/cdelta2                        , [0.0   ,    1.0000,    0.6000]),
                    ((cref-clevel[0]+(clevel[-1]-cref)*0.5)/cdelta2    , [0.2000,    0.6000,    1.0000]),
                    ((clevel[-1]-clevel[0])/cdelta2                    , [1.0000,    1.0000,    1.0000])]        
    
    #-->________________________________________________________________________
    elif cname=='odv':    
        cmap_def = [(0.0                                                , [0.9373,    0.7765,    0.9373]),
                    (((cref-clevel[0])*0.25)/cdelta2                , [0.7804,    0.3647,    0.7490]),
                    (((cref-clevel[0])*0.50)/cdelta2                , [0.1922,    0.2235,    1.0000]),
                    (((cref-clevel[0])*0.75)/cdelta2                , [0.4824,    0.9686,    0.8706]),
                    ((cref-clevel[0])/cdelta2                        , [0.4980,    1.0000,    0.4980]),
                    ((cref-clevel[0]+(clevel[-1]-cref)*0.25)/cdelta2, [1.0000,    0.7843,    0.1373]),
                    ((cref-clevel[0]+(clevel[-1]-cref)*0.50)/cdelta2, [1.0000,       0.0,       0.0]),
                    ((cref-clevel[0]+(clevel[-1]-cref)*0.75)/cdelta2, [0.8392,    0.0627,    0.1922]),
                    ((clevel[-1]-clevel[0])/cdelta2                    , [1.0000,    0.7765,    0.5804])]
    #-->________________________________________________________________________
    elif cname=='odv_w':    
        cmap_def = [(0.0                                            , [0.9373,    0.7765,    0.9373]),
                    (((cref-clevel[0])*0.25)/cdelta2                , [0.7804,    0.3647,    0.7490]),
                    (((cref-clevel[0])*0.50)/cdelta2                , [0.1922,    0.2235,    1.0000]),
                    (((cref-clevel[0])*0.75)/cdelta2                , [0.4824,    0.9686,    0.8706]),
                    ((cref-clevel[0])/cdelta2                        , [1.0   ,    1.0000,    1.0   ]),
                    ((cref-clevel[0]+(clevel[-1]-cref)*0.25)/cdelta2, [1.0000,    0.7843,    0.1373]),
                    ((cref-clevel[0]+(clevel[-1]-cref)*0.50)/cdelta2, [1.0000,       0.0,       0.0]),
                    ((cref-clevel[0]+(clevel[-1]-cref)*0.75)/cdelta2, [0.8392,    0.0627,    0.1922]),
                    ((clevel[-1]-clevel[0])/cdelta2                    , [1.0000,    0.7765,    0.5804])]
            
    #-->________________________________________________________________________
    elif cname=='wvt':    
        cmap_def = [(0.0                                              , [255.0/255, 255.0/255, 255.0/255]),
                    (((cref-clevel[0])*0.3333)/cdelta2                  , [255.0/255, 255.0/255, 153.0/255]),
                    (((cref-clevel[0])*0.6666)/cdelta2                  , [255.0/255, 204.0/255,  51.0/255]),
                    ((cref-clevel[0])/cdelta2                          , [255.0/255, 177.0/255, 100.0/255]),
                    ((cref-clevel[0]+(clevel[-1]-cref)*0.3333)/cdelta2, [255.0/255, 102.0/255, 102.0/255]),
                    ((cref-clevel[0]+(clevel[-1]-cref)*0.6666)/cdelta2, [255.0/255,  51.0/255,  51.0/255]),
                    ((clevel[-1]-clevel[0])/cdelta2                      , [153.0/255,   0.0/255,  51.0/255])]
            
    #-->________________________________________________________________________
    elif cname=='seaice':    
        cmap_def = [(0.0                                            , [153.0/255,   0.0/255,  51.0/255]),
                    (((cref-clevel[0])*0.25)/cdelta2                , [204.0/255,   0.0/255,   0.0/255]),
                    (((cref-clevel[0])*0.50)/cdelta2                , [255.0/255, 102.0/255, 102.0/255]),
                    (((cref-clevel[0])*0.75)/cdelta2                , [255.0/255, 153.0/255, 153.0/255]),
                    ((cref-clevel[0])/cdelta2                        , [255.0/255, 255.0/255, 255.0/255]),
                    ((cref-clevel[0]+(clevel[-1]-cref)*0.25)/cdelta2, [153.0/255, 255.0/255, 255.0/255]),
                    ((cref-clevel[0]+(clevel[-1]-cref)*0.50)/cdelta2, [  0.0/255, 153.0/255, 204.0/255]),
                    ((cref-clevel[0]+(clevel[-1]-cref)*0.75)/cdelta2, [  0.0/255,  51.0/255, 204.0/255]),
                    ((clevel[-1]-clevel[0])/cdelta2                    , [  0.0/255,  51.0/255, 153.0/255])]
    
    #-->________________________________________________________________________
    elif cname=='seaice_i':    
        cmap_def = [(0.0                                            , [  0.0/255,  51.0/255, 153.0/255]),
                    (((cref-clevel[0])*0.25)/cdelta2                , [  0.0/255,  51.0/255, 204.0/255]),
                    (((cref-clevel[0])*0.50)/cdelta2                , [  0.0/255, 153.0/255, 204.0/255]),
                    (((cref-clevel[0])*0.75)/cdelta2                , [153.0/255, 255.0/255, 255.0/255]),
                    ((cref-clevel[0])/cdelta2                        , [255.0/255, 255.0/255, 255.0/255]),
                    ((cref-clevel[0]+(clevel[-1]-cref)*0.25)/cdelta2, [255.0/255, 153.0/255, 153.0/255]),
                    ((cref-clevel[0]+(clevel[-1]-cref)*0.50)/cdelta2, [255.0/255, 102.0/255, 102.0/255]),
                    ((cref-clevel[0]+(clevel[-1]-cref)*0.75)/cdelta2, [204.0/255,   0.0/255,   0.0/255]),
                    ((clevel[-1]-clevel[0])/cdelta2                    , [153.0/255,   0.0/255,  51.0/255])]
            
    ##-->________________________________________________________________________
    elif cname=='test':    
        cmap_def = [(0.0                    , [0.0, 0.0, 0.0]),
                    (0.5                    , [0.5, 0.5, 0.5]),
                    (1.0                    , [1.0, 1.0, 1.0])]        
    
    ##-->________________________________________________________________________
    #elif cname.find('cmocean.cm')==0:
        #nmax = 9
        ##cdict = cmocean.tools.get_dict(cmocean.cm.balance, N=nmax)
        #cdict = cmocean.tools.get_dict(eval(cname), N=nmax)
        #aux_icref = [0.0, 
                    #((cref-clevel[0])*0.25)/cdelta2,
                    #((cref-clevel[0])*0.50)/cdelta2,
                    #((cref-clevel[0])*0.75)/cdelta2,
                    #(cref-clevel[0])/cdelta2,
                    #(cref-clevel[0]+(clevel[-1]-cref)*0.25)/cdelta2,
                    #(cref-clevel[0]+(clevel[-1]-cref)*0.50)/cdelta2,
                    #(cref-clevel[0]+(clevel[-1]-cref)*0.75)/cdelta2,
                    #(clevel[-1]-clevel[0])/cdelta2]
                    
        #for ii in range(0,nmax):
            #auxr = list(cdict['red'][ii])
            #auxr[0]=aux_icref[ii]
            #cdict['red'][ii]=tuple(auxr)
            
            #auxr = list(cdict['blue'][ii])
            #auxr[0]=aux_icref[ii]
            #cdict['blue'][ii]=tuple(auxr)
            
            #auxr = list(cdict['green'][ii])
            #auxr[0]=aux_icref[ii]
            #cdict['green'][ii]=tuple(auxr)
            
    #-->________________________________________________________________________
    #elif cname=='jet':    
        #cmap_def = [(0                                                , []),
                    #(((cref-clevel[0])*0.25)/cdelta2                , []),
                    #(((cref-clevel[0])*0.50)/cdelta2                , []),
                    #(((cref-clevel[0])*0.75)/cdelta2                , []),
                    #((cref-clevel[0])/cdelta2                        , []),
                    #((cref-clevel[0]+(clevel[-1]-cref)*0.25)/cdelta2, []),
                    #((cref-clevel[0]+(clevel[-1]-cref)*0.50)/cdelta2, []),
                    #((cref-clevel[0]+(clevel[-1]-cref)*0.75)/cdelta2, []),
                    #((clevel[-1]-clevel[0])/cdelta2                    , [])]
        
    #___________________________________________________________________________
    # make python colormap
    #print(cmap_def)
    #print(cmap_def[0][1])
    #print(cmap_def[1][1])
    
    #cdict = [(0.0,  0.0, 0.0),
             #(0.5,  0.5, 0.5),
             #(1.0,  1.0, 1.0)]

    #if cname.find('cmocean.cm')==0:    
        #cmap = LinearSegmentedColormap(cname, cdict ,N=clevel.size-1)
        #cmap.set_under([cdict['red'][0][1],cdict['green'][0][1],cdict['blue'][0][1]] )
        #cmap.set_over( [cdict['red'][-1][1],cdict['green'][-1][1],cdict['blue'][-1][1]] )
    #else:
        #cmap = LinearSegmentedColormap.from_list(cname, cmap_def, N=clevel.size-1, gamma=1)
        #cmap.set_under(cmap_def[0][1])
        #cmap.set_over(cmap_def[-1][1])
    cmap = LinearSegmentedColormap.from_list(cname, cmap_def, N=clevel.size-1, gamma=1)
    cmap.set_under(cmap_def[0][1])
    cmap.set_over(cmap_def[-1][1])
    
    return(cmap,clevel)
