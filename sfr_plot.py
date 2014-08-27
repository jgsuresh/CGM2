import numpy as np
import h5py
import matplotlib
matplotlib.use('PDF')
from matplotlib import rc
# rc('text', usetex=True)
rc('font', family='serif')
import matplotlib.pyplot as plt
import brewer2mpl

from units import AREPO_units
AU = AREPO_units()

class sfr_plot:
    def __init__(self):
        #self.npz_base = "/n/home04/jsuresh/CGM_new/data/npz/"
        self.fig_base = "/n/home04/jsuresh/CGM_new/data/figs/"

        bmap_o = brewer2mpl.get_map('Spectral','Diverging',9, reverse=True)
        bmap = bmap_o.mpl_colormap
        
        mv_base = "/n/hernquistfs1/Illustris/SmallBox/GFM/Production/Cosmo/" 
        # mv_base = mv_base+""

        c0_512_base = mv_base+"Cosmo0_V6/L25n512/output/"
        c0_256_base = mv_base+"Cosmo0_V6/L25n256/output/"
        c2_base = mv_base+"Cosmo2_V6/L25n256/output/"
        c4_512_base = "/n/hernquistfs1/spb/Cosmo/Cosmo4_V6/L25n512/output/"
        c3_512_base = "/n/hernquistfs1/spb/Cosmo/Cosmo3_V6/L25n512/output/"
        c5_256_base = "/n/home04/jsuresh/runs/Cosmo5_V6/output/"
        c0_fw_256_base = mv_base+"Cosmo0_V6_fastWinds/L25n256/output/"
        c0_sw_256_base = mv_base+"Cosmo0_V6_strongWinds/L25n256/output/"

        gam_10_BH_base = "/n/hernquistfs1/jsuresh/Runs/gam_10_BH/output/"
        gam_20_BH_base = "/n/hernquistfs1/jsuresh/Runs/gam_20_BH/output/"
        gam_30_BH_base = "/n/hernquistfs1/jsuresh/Runs/gam_30_BH/output/"
        gam_40_BH_base = "/n/hernquistfs1/jsuresh/Runs/gam_40_BH/output/"
        gam_50_BH_base = "/n/hernquistfs1/jsuresh/Runs/gam_50_BH/output/"

        gam_10_BH_nothermal_base = "/n/hernquistfs1/jsuresh/Runs/gam_10_BH_nothermal/output/"
        gam_20_BH_nothermal_base = "/n/hernquistfs1/jsuresh/Runs/gam_20_BH_nothermal/output/"
        gam_30_BH_nothermal_base = "/n/hernquistfs1/jsuresh/Runs/gam_30_BH_nothermal/output/"
        gam_40_BH_nothermal_base = "/n/hernquistfs1/jsuresh/Runs/gam_40_BH_nothermal/output/"
        gam_50_BH_nothermal_base = "/n/hernquistfs1/jsuresh/Runs/gam_50_BH_nothermal/output/"

        gam_25_BH_base = "/n/hernquistfs1/jsuresh/Runs/gam_25_BH/output/"
        gam_25_noBH_base = "/n/hernquistfs1/jsuresh/Runs/gam_25_noBH/output/"
        gam_50_BH_base = "/n/hernquistfs1/jsuresh/Runs/gam_50_BH/output/"
        gam_50_noBH_base = "/n/hernquistfs1/jsuresh/Runs/gam_50_noBH/output/"
        gam_75_BH_base = "/n/hernquistfs1/jsuresh/Runs/gam_75_BH/output/"
        gam_75_noBH_base = "/n/hernquistfs1/jsuresh/Runs/gam_75_noBH/output/"
        gam_95_BH_base = "/n/hernquistfs1/jsuresh/Runs/gam_95_BH/output/"
        gam_95_noBH_base = "/n/hernquistfs1/jsuresh/Runs/gam_95_noBH/output/"

        gam_25_fixv_base = "/n/hernquistfs1/jsuresh/Runs/gam_25_fixv/output/"
        gam_50_fixv_base = "/n/hernquistfs1/jsuresh/Runs/gam_50_fixv/output/"
        gam_75_fixv_base = "/n/hernquistfs1/jsuresh/Runs/gam_75_fixv/output/"
        gam_95_fixv_base = "/n/hernquistfs1/jsuresh/Runs/gam_95_fixv/output/"

        gam_50_fixv_nothermal_base = "/n/hernquistfs1/jsuresh/Runs/gam_50_fixv_nothermal/output/"

        c0_nometalwinds_base = "/n/hernquistfs1/jsuresh/Runs/c0_nometalwinds/output/"
        c0_fullmetalwinds_base = "/n/hernquistfs1/jsuresh/Runs/c0_fullmetalwinds/output/"

        c0_check_base = "/n/hernquistfs1/jsuresh/Runs/c0_check/output/"
        c2_check_base = "/n/hernquistfs1/jsuresh/Runs/c2_check/output/"
        
        #self.base_list = [c0_512_base,c4_512_base,c3_512_base,c5_256_base,c0_fw_256_base,c0_sw_256_base]
        #self.label_list = ["c0","c4","c3","c5","fw","sw"]
        #self.base_list = [c0_256_base,c0_512_base,gam_25_BH_base,gam_50_noBH_base,gam_75_noBH_base,gam_95_noBH_base]
        #self.label_list = ["c0 256","c0 512","gam 25 BH","gam 50 noBH","gam 75 noBH","gam 95 noBH"]
        #self.base_list = [c0_256_base,c0_512_base,gam_25_BH_base,gam_50_BH_base,gam_75_BH_base,gam_95_BH_base]
        #self.label_list = ["c0 256","c0 512","gam 25 BH","gam 50 BH","gam 75 BH","gam 95 BH"]

        # self.base_list = [c0_256_base,gam_25_BH_base,gam_50_BH_base,gam_75_BH_base,gam_95_BH_base,gam_25_BH_base,gam_50_noBH_base,gam_75_noBH_base,gam_95_noBH_base]
        # self.label_list = ["c0 256","gam 25","gam 50","gam 75","gam 95",None,None,None,None]
        # self.color_list = ["black","blue","yellow","orange","red","blue","yellow","orange","red"]
        # self.ls_list = ["solid","solid","solid","solid","solid","dotted","dotted","dotted","dotted"]

        # self.base_list = [c0_256_base,gam_25_BH_base,gam_50_BH_base,gam_75_BH_base,gam_95_BH_base,gam_25_BH_base,gam_50_noBH_base,gam_75_noBH_base,gam_95_noBH_base]
        # self.label_list = ["c0 256","gam 25","gam 50","gam 75","gam 95",None,None,None,None]
        # self.color_list = ["black","blue","yellow","orange","red","blue","yellow","orange","red"]
        # self.ls_list = ["solid","solid","solid","solid","solid","dotted","dotted","dotted","dotted"]

        # self.base_list = [c0_256_base,gam_25_fixv_base,gam_50_fixv_base,gam_75_fixv_base,gam_95_fixv_base]
        # self.label_list = ["c0 256","gam 25 fixv","gam 50 fixv","gam 75 fixv","gam 95 fixv"]
        # #self.color_list = ["black","blue","yellow","orange","red"]
        # self.color_list = [bmap(0),bmap(0.2),bmap(0.4),bmap(0.6),bmap(0.8)]
        # self.ls_list = ["solid","solid","solid","solid","solid"]

        # self.base_list = [c0_256_base,gam_10_BH_base,gam_20_BH_base,gam_30_BH_base,gam_40_BH_base,gam_50_BH_base]
        # self.label_list = ["c0 256","gam 10","gam 20","gam 30","gam 40","gam 50","gam 60"]
        # self.color_list = [bmap(0),bmap(0.2),bmap(0.4),bmap(0.6),bmap(0.8),bmap(1.0)]
        # self.ls_list = ["solid","solid","solid","solid","solid","solid"]

        # self.base_list = [c0_256_base,gam_10_BH_base,gam_10_BH_nothermal_base,gam_20_BH_base,gam_20_BH_nothermal_base,gam_30_BH_base,gam_30_BH_nothermal_base,gam_40_BH_base,gam_40_BH_nothermal_base,gam_50_BH_base,gam_50_BH_nothermal_base]
        # self.label_list = ["c0 256","gam 10",None,"gam 20",None,"gam 30",None,"gam 40",None,"gam 50",None]
        # self.color_list = [bmap(0),bmap(0.2),bmap(0.2),bmap(0.4),bmap(0.4),bmap(0.6),bmap(0.6),bmap(0.8),bmap(0.8),bmap(1.0),bmap(1.0)]
        # self.ls_list = ["solid","solid","dotted","solid","dotted","solid","dotted","solid","dotted","solid","dotted"]

        # self.base_list = [c0_256_base]
        # self.label_list = ["c0_256"]
        # self.color_list = [bmap(0)]
        # self.ls_list = ["solid"]

        self.base_list = [c0_256_base,c2_base,c0_sw_256_base,gam_50_fixv_nothermal_base,gam_50_BH_base,gam_50_fixv_base,c0_nometalwinds_base,c0_fullmetalwinds_base]
        self.label_list = ["Fiducial","No AGN","Higher Mass-Loading","Faster Winds","Fixed-E Hot Winds","Fixed-v Hot Winds","Pristine Winds","Fully Enriched Winds"]
        self.color_list = ["blue","cyan","chartreuse","gold","saddlebrown","magenta","gray","black"]
        self.ls_list = ["solid","solid","solid","solid","solid","solid","solid","solid"]
        self.extract_sfr_data()


    def extract_sfr_data(self):
        plt.close('all')
        plt.figure(figsize=(3.54331,3.14))
        ax1 = plt.subplot(1,1,1)
        for i in xrange(len(self.base_list)):
            print "Working on {}".format(self.label_list[i])
            base = self.base_list[i]
            fname = base+"sfr.txt"

            dat = np.loadtxt(fname)
            a = dat[:,0]
            z = 1./a - 1.
            print "z: ",z
            sfr = dat[:,2]

            #x = self.ss_dat(np.log10(1+z),n=50)
            #y = self.ss_dat(np.log10(sfr),n=50)
            x = self.ss_dat(np.log10(1+z),n=50)
            y = self.ss_dat(sfr/(25/0.7)**3,n=50)
            # y = self.ss_dat(AU.PhysicalSFRD(sfr,z),n=50) #AU.PhysicalVolume(25**3.,z)
            # y = self.ss_dat(AU.PhysicalMdot(sfr)*0.7**3./25**3.,n=50)
            # y = self.ss_dat(sfr,n=50)
            print "y ",y
            plt.semilogy(x,y,label=self.label_list[i],color=self.color_list[i],linestyle=self.ls_list[i])
        plt.legend(prop={'size':5.7},ncol=2)
        plt.xlim([0.45,1.0])
        plt.ylim([10.**-3.,1.0])
        #plt.ylabel(r' SFRD $\left[M_\odot$ yr$^{-1}$ Mpc$^{-3} \right]$')
        plt.ylabel(r'SFRD $[$'+'M'+'$_\odot$'+' yr'+r'$^{-1}$'+' cMpc'+r'$^{-3}]$')
        plt.xlabel(r'Log$_{10} [1+z]$')
        plt.subplots_adjust(left=0.2,bottom=0.18,top=0.85)

        ax2 = ax1.twiny()
        u = 1/0.55
        new_tick_locations = (1+u)*(np.array([0.47712,0.6021,0.69897,0.8450,1.0]))-u
        # the 0.8181 comes from the scaling of the x-axis to [0.45,1] instead of [0,1]

        def tick_function(X):
            V_foo = (X+u)/(1+u)
            V = 10.**(V_foo) - 1.
            V = np.rint(V)
            return ["{}".format(int(z)) for z in V]

        ax2.set_xticks(new_tick_locations)
        ax2.set_xticklabels(tick_function(new_tick_locations))
        ax2.set_xlabel(r"$z$")
        plt.show()

        # print "saving to {}".format(self.fig_base+"sfr_plot.pdf")
        plt.savefig(self.fig_base+"sfr_plot.pdf")

    def ss_dat(self,dat,n=2):
        s = np.size(dat)
        ind = np.arange(0,s-1,n)
        return dat[ind]


if __name__ == '__main__':
    sfr_plot()




