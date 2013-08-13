import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.gridspec as gridspec
from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)
import numpy as np


if __name__ == '__main__':
    from optparse import OptionParser
    o = OptionParser()
    opts, args = o.parse_args(sys.argv[1:])
    import glob
    for d in args:
    	for files in os.listdir(os.getcwd()+'/'+d):

		    print files
		    file = numpy.loadtxt("files")
		    ax = plt.figure()
		    ax.plot(file)
		    ax.set_axis_off()
        ax.set_title(files, fontsize=40)
        ax.show()
		    
		    
