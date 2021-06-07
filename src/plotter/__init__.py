from astropy.table import Table, Column
import matplotlib.pyplot as plt
import matplotlib.widgets as wgt

from templates import Spectra2nd


class Plotter:
    def __init__(self):
        self.spectra = Spectra2nd()
        self.lgLTIR = 10.0
        self.lgsSFR = -10.0

        self.fig, self.ax = plt.subplots()
        plt.get_current_fig_manager().set_window_title('Template plotter')
        plt.subplots_adjust(bottom=0.20)

        lgLTIR_ax = plt.axes([0.2, 0.02, 0.56, 0.0275])
        lgLTIR = wgt.Slider(
            lgLTIR_ax,
            'log LTIR',
            9.0,
            12.0,
            valinit=self.lgLTIR,
            valfmt=" %5.2F"
        )
        lgLTIR.on_changed(self.updateLTIR)

        lgsSFR_ax = plt.axes([0.2, 0.06, 0.56, 0.0275])
        lgsSFR = wgt.Slider(
            lgsSFR_ax,
            'log sSFR',
            -11.5,
            -9.0,
            valinit=self.lgsSFR,
            valfmt="%6.2F"
        )
        lgsSFR.on_changed(self.updatesSFR)

        saveFITS_ax = plt.axes([0.86, 0.0625, 0.135, 0.0375])
        saveFITS = wgt.Button(saveFITS_ax, label='Save FITS')
        saveFITS.on_clicked(self.toFITS)

        savePDF_ax = plt.axes([0.86, 0.02, 0.135, 0.0375])
        savePDF = wgt.Button(savePDF_ax, label='Save PDF')
        savePDF.on_clicked(self.toPDF)

        self.update_plot()
        plt.show()

    def updateLTIR(self, val):
        self.lgLTIR = val
        self.update_plot()

    def updatesSFR(self, val):
        self.lgsSFR = val
        self.update_plot()

    def update_plot(self):
        data = self.spectra.estimate(('LTIR', 'sSFR'), self.lgLTIR, self.lgsSFR)
        if len(self.ax.get_lines()) > 0:
            self.plot.set_ydata(data)
        else:
            (self.plot,) = self.ax.plot(self.spectra.wl, data, c='k', lw=1.5)
            self.ax.semilogx()
            self.ax.set_xlim(3000.0, 5e5)
            self.ax.set_ylim(7.0, 12.0)
            self.ax.minorticks_on()
            self.ax.set_xlabel('Wavelength [nm]')
            self.ax.set_ylabel(r'log $\nu L_\nu$ [$L_\odot$]')
            self.fig.canvas.draw()

    def toFITS(self, _):
        wl = Column(self.spectra.wl, name='wavelength', unit='nm')
        spectrum = Column(
            10 ** self.spectra.estimate(("LTIR", "sSFR"), self.lgLTIR, self.lgsSFR
        ),
            name='nuLnu', unit='Lsun')

        table = Table()
        table.add_columns([wl, spectrum])
        table.write('spectrum.fits')

    def toPDF(self, _):
        plt.savefig('spectrum.pdf')


def main():
    _ = Plotter()
