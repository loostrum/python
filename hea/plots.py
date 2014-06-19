from matplotlib import pyplot as plt
from matplotlib import rc
import numpy as np
import math

def MJ(gamma):
    theta=.168637 # =kt/mc^2
    const = 3238.1 # = 1/(theta * K_2(1/theta))
    return const*math.sqrt(1-gamma**-2)*gamma**2*math.exp(-gamma/theta)

        
def create_MJ():
    #create binned MJ 
    #electronbins = binned gamma
    electronbins= np.linspace(1,15,1E3)
    mjdist= [ MJ(i) for i in electronbins ]
    return electronbins,mjdist

def main():    
    electronbins,mjdist = create_MJ()

    d='data'
    e_bins = np.loadtxt(d+'/bins.txt')
    n_sync = np.loadtxt(d+'/sync-input.txt')
    n_planck = np.loadtxt(d+'/planck-input.txt')
    n_non = np.loadtxt(d+'/non-scattered.txt')
    n_comp_sy = np.loadtxt(d+'/sync-output.txt')
    n_comp_pl = np.loadtxt(d+'/planck-output.txt')

    full_input = np.add(n_planck,n_sync)
    comp_output = np.add(n_comp_pl,n_comp_sy)
    full_output = np.add(comp_output,n_non)


    def plotter(y,col,lab,style='-'):
        plt.loglog(e_bins,y,color=col,ls=style,label=lab)


    def plot_finish():
        plt.legend()
        plt.xlim(1E-2,1E4)
        plt.ylim(1E-35,1E-30)
        plt.xlabel(r'E (keV)')
        plt.ylabel(r'$F_\nu$' )

    #plot MJ distribution
    plt.plot(electronbins,mjdist,color='red')
    plt.title('Maxwell-Juttner distribution')
    plt.xscale('log')
    plt.xlabel(r'$\gamma$')
    plt.ylabel(r'$P(\gamma)$')
    plt.xlim(1,5)
    plt.ylim(1E-6,3)
    plt.savefig('mj.png')
    plt.show()


    #plot inputs
    plotter(n_planck,'red','Thermal input','--')
    plotter(n_sync,'blue','Synchrotron input','--')
    plotter(full_input,'black','Total input')
    plt.title('Input spectrum')
    plot_finish()
    plt.savefig('input.png')
    #plt.show()

    plt.close()

    #plot outputs
    plotter(n_non,'green','Non-scattered output','--')
    plotter(n_comp_pl,'red','Comptonized thermal output','--')
    plotter(n_comp_sy,'blue','Comptonized synchrotron output','--')
    plotter(full_output,'black','Total output',)
    plot_finish()
    plt.savefig('output_parts.png')
    #plt.show()

    plt.close()

    #plot total input and output
    plotter(full_output,'black','Total output')
    plotter(full_input,'red','Total input','--')
    plot_finish()
    plt.savefig('output_total.png')

    #plt.show()

if __name__=='__main__':
    main()
