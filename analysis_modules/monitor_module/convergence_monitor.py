import json
from abci_monitor import Monitor
import matplotlib.pyplot as plt
from matplotlib.ticker import LinearLocator, FormatStrFormatter
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.animation as animation
from matplotlib import style
abci_mon = Monitor()

class Monitor:
    def __init__(self):
        self.n_list = []
        self.n = 0

    def objective_function(self):
        pass

    def get_quantities(self, code, monitor_list):
        # get data from excel or json saved
        if code.lower() == 'abci':
            data = abci_mon.get_convergence_data(monitor_list)
            return data
        else:
            pass

    def animate(self, i):
        data = self.get_quantities('abci', ['k_loss_long', 'k_loss_trans', 'df_M0D1'])

        ax1.clear()
        for key, val in data.items():
            ax1.plot(val)

        plt.gca().set_prop_cycle(None)

        self.n += 1


if __name__ == '__main__':
    monitor = Monitor()
    style.use('ggplot')

    fig = plt.figure()
    ax1 = fig.add_subplot(1, 1, 1)

    ani = animation.FuncAnimation(fig, monitor.animate, interval=1000)
    plt.show()