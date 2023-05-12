# -​*- coding: utf-8 -*​-
# @author: Laurent Hébert-Dufresne <lhebertd@uvm.edu>

# Packages
import numpy as np
import matplotlib
import matplotlib.gridspec
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm

from matplotlib.offsetbox import (TextArea, DrawingArea, OffsetImage,
                                  AnnotationBbox)
from matplotlib.cbook import get_sample_data


# Global parameters for the figures.
matplotlib.use('agg')
plt.style.use('seaborn-deep')
plt.rcParams["text.usetex"] = True
plt.rcParams["font.size"] = 24
plt.rcParams["legend.frameon"] = False
plt.rcParams["legend.fancybox"] = False


fig = plt.figure(constrained_layout=True, figsize=(7, 6))

# Add rows or columns here for panels
gs = matplotlib.gridspec.GridSpec(nrows=1,
                                  ncols=1,
                                  figure=fig,
                                  wspace=0.05, hspace=0.05
                                  )
ax1 = fig.add_subplot(gs[0, 0])

#Load simulation results.
y = np.loadtxt('./results/results.dat')
z = np.loadtxt('./matrices/Pnkl.txt')
num_types, num_cols = z.shape

colors = plt.cm.Spectral_r(np.linspace(0,1,6))
# Compute and plot final size
ax1.set_prop_cycle('color', colors)
if num_types>13:
  color = iter(cm.rainbow(np.linspace(0, 1, 13)))
  plots = np.random.choice(np.arange(0,num_types), size=13)
  for m in plots:
      c = next(color)
      ax1.plot(y[:,0],y[:,int(m+1)], c=c, label=r'$n=$ '+str(int(z[m,0]))+', $k=$ '+str(int(z[m,1])) +', $\ell=$ '+str(int(z[m,2])))
else:
  for m in range(0,num_types):
      color = iter(cm.rainbow(np.linspace(0, 1, num_types)))
      c = next(color)
      ax1.plot(y[:,0],y[:,int(m+1)], c=c, label=r'$n=$ '+str(int(z[m,0]))+', $k=$ '+str(int(z[m,1])) +', $\ell=$ '+str(int(z[m,2])))

# Legend.
lines, labels = ax1.get_legend_handles_labels()
ax1.legend(lines, labels, loc=(1.0,0.0), shadow=False, fancybox=False, prop={'size':20}, frameon=False, handlelength=0.75, numpoints=1)

# Labels
#ax1.set_yscale('log')
ax1.set_ylabel(r'Prevalence per layer')
ax1.set_xlabel(r'Transmission $\lambda$')
ax1.set_xlim(0.0,0.6)
ax1.set_ylim(1e-4,1)

# Save to file.
fig.savefig("./figures/Casestudy_formatted.pdf", bbox_inches='tight')
fig.savefig("./figures/Casestudy_formatted.svg", bbox_inches='tight')
fig.savefig("./figures/Casestudy_formatted.png", bbox_inches='tight')
