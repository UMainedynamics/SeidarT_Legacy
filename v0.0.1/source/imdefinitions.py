#!/usr/bin/env python3

import numpy as np

# ---------------------------- Processing functions ---------------------------
def agc(ts, k, agctype):
    # Auto-gain normalization using a running window
    # k is the window length
    # agctype is either mean, rms, median 
    n = len(ts)

    k2 = int(k/2) # This will floor the number if k is even; (k-1)/2 
    if np.mod(k, 2) == 0: # even numbers need will not have a centered window
        k = int( k + 1) 
    
    stat = np.ones([n])
    # Normalize
    if agctype == "std":
        for i in range(0, k2):
            stat[i] = np.std( abs( ts[0:(i+k2)] ) )
            stat[ (n-k2+i) ] = np.std( abs(ts[ (n-2*k2+i):n] ) )
        for i in range(k2,n-k2):
            stat[i] = np.std( abs( ts[ (i-k2):(i+k2) ] ) )
    elif agctype == "mean":
        for i in range(0, k2):
            stat[i] = np.mean( abs( ts[0:(i+k2)] ) )
            stat[ (n-k2+i) ] = np.mean( abs(ts[ (n-2*k2+i):n] ) )
        for i in range(k2,n-k2):
            stat[i] = np.mean( abs( ts[ (i-k2):(i+k2) ] ) )
    else:
        for i in range(0, k2):
            stat[i] = np.std( ts[i:(i+k2)] )
            stat[ (n-k2+i) ] = np.std( ts[ (n-2*k2+i):n] )
        for i in range(k2,n-k2):
            stat[i] = np.std( ts[ (i-k2):(i+k2) ] )
    
    stat[stat == 0] = 1
    ts = ts/stat 
    return ts



# =============================================================================
# =============================================================================

# See rcxdisplay 
# def tsplot(self):
#         m, n = self.timeseries.shape
#         # if the gain is 0, then the window is 1
        
#         # The default is a NoneType so we won't apply gain
#         if self.gain is None:
#             self.gain = m
        
#         if self.gain == 0:
#             self.gain = int(1)
        
#         # The gain can't exceed the length of the time series
#         if self.gain > m:
#             self.gain = m 

#         for ind in range(0, n):
#             self.timeseries[:,ind] = agc(
#                 self.timeseries[:,ind], self.gain, "mean"
#             )
        
#         # Create the values for the y-axis
#         timelocs = np.arange(0,m, int(m/10) ) #10 tick marks
#         timevals = timelocs*self.dt
        
#         # Create the reciever location
#         xlocs = np.arange(0, n, int(n/5) ) #5 tick marks 
        
#         # Create the figure
#         fig = plt.figure(figsize =(n/2,m/2) )
#         ax1 = plt.gca()
#         # ax2 = ax1.twinx() 
#         ax1.imshow(self.timeseries, cmap = 'Greys', aspect = 'auto')
#         ax1.set_xlabel(r'Receiver #')
#         ax1.xaxis.tick_top()
#         ax1.xaxis.set_label_position('top')
#         ax1.set_xticks(xlocs)
#         ax1.set_xticklabels(xlocs)
#         ax1.set_ylabel(r'Two-way Travel Time (s)')
#         ax1.set_yticks(timelocs)
#         ax1.set_yticklabels(timevals)
        
#         # ax2.set_ylabel('Depth (m)')
#         # ax2.set_yticks(timelocs) 
#         # ax2.set_yticklabels(twt)
        
#         # # Label the x-axis
#         # ax.xaxis.tick_top()
#         # ax.xaxis.set_label_position('top')
#         # ax.set_xlabel(r"Receiver #")
        
#         # Label the y-axis
#         if self.channel == 'Vx' or self.channel == 'Vz':
#             mult = 1e2
#         else:
#             mult = 1e6
        
#         time_locations = np.linspace(1, m, 10)
#         time_labels = np.round( time_locations*self.dt*mult, 4)
#         # ax.set_ylabel(r"Two way travel time (s)")
#         plt.yticks(ticks=time_locations, labels = time_labels.astype(str) )
        
#         # # Other figure operations
#         # if self.channel == 'Vx' or self.channel == 'Vz':
#         #     plt.figtext(0.30, 0.07, 'x $10^{-2}$')	
#         # else:
#         #     plt.figtext(0.30, 0.07, 'x $10^{-6}$')
        
#         ax1.set_aspect(aspect = exaggeration	)
#         plt.show()
    
#     def plot_layout(self):
#         figmod, axmod = plt.subplots()
#         img = mpimg.imread(self.modelfile)
#         axes_extent = [0, img.shape[2],img.shape[0], 0 ]
        
#         implot = axmod.imshow(img, extent = extent)
        
#         # add white upside down triangle receivers
#         axmod.scatter(
#             self.receiver_locations[:,0]- self.cpml,
#             self.receiver_locations[:,1] - self.cpml,
#             marker = 'v', s = 30, c = (0.8, 0.8, 0.8, 1),
#             linewidths = 0.5, edgecolor = (0.2, 0.2, 0.2, 1 )
#         )
        
#         # add white star source
#         axmod.scatter(
#             self.source_location[0], self.source_location[1],
#             marker = '*', s = 30, c = (0.8, 0.8, 0.8, 1 ),
#             linewidths = 1, edgecolor = (0.2, 0.2, 0.2, 1 ) 
#         )
        
#         axmod.set_ylabel("z-indice")
#         axmod.set_xlabel("x-indice")
#         plt.show()


# =============================================================================
# =============================================================================