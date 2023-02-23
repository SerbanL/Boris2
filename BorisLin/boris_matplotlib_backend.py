#custom matplotlib backend so we can redirect matplotlib output to image file
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.backend_bases import FigureManagerBase
import matplotlib.pyplot as plt

class FigureManager(FigureManagerBase):
    fig_num = 0
    def show(self):
        pass

FigureCanvas = FigureCanvasAgg

def show(*args, **kwargs):
    figmanager = plt.get_current_fig_manager()
    plt_fileName = "plt_figure_%d.png" % (figmanager.fig_num)
    print(f"Matplotlib figure directed to image file: {plt_fileName}")
    figmanager.canvas.figure.savefig(plt_fileName)
    figmanager.canvas.figure.clear()
    figmanager.fig_num += 1