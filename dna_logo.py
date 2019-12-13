# a class for ploting dna logos
import seaborn
import matplotlib
import matplotlib.pyplot as plt
plt.style.use('seaborn-ticks')

from matplotlib import transforms
import matplotlib.patheffects
import numpy as np

matplotlib.use('AGG')

class Logo:
    COLOR_SCHEME = {'G': 'orange',
                'A': 'red',
                'C': 'blue',
                'T': 'darkgreen'}
    BASES = list(COLOR_SCHEME.keys())

    class Scale(matplotlib.patheffects.RendererBase):
        def __init__(self, sx, sy=None):
            self._sx = sx
            self._sy = sy

        def draw_path(self, renderer, gc, tpath, affine, rgbFace):
            affine = affine.identity().scale(self._sx, self._sy)+affine
            renderer.draw_path(gc, tpath, affine, rgbFace)

    # count_mat : k x 4 numpy matrix, count of each base at each position
    # pfm: k x 4 matrix, frequency, each row sum up to 1
    # pwm: k x 4 matrix, for making logos
    def __init__(self, count_mat=None, pfm=None, pwm=None, out_logo_file='out_logo.png'):
        self.count_mat = count_mat
        self.pfm = pfm
        self.pwm = pwm
        self.out_logo_file = out_logo_file
        # if not background_distribution:
        #     self.background_distribution = np.array([0.25, 0.25, 0.25, 0.25])
        if not (count_mat is None):
            self.pwm = self.count2pwm()
        if not (pfm is None):
            self.pwm = self.pfm2pwm()
        
    def count2pwm(self):
        if not self.pfm:
            self.pfm = self.count2pfm()
        self.pwm = self.pfm2pwm()
        return self.pwm
    
    def count2pfm(self):
        if not (self.count_mat is None):
            counts = self.count_mat + 1e-6
            self.pfm = self._normalize(counts)
            return self.pfm
        else:
            raise ValueError('self.count_mat is None') 
    
    @staticmethod
    def _normalize(in_mat):
        # in_mat: a k x 4 numpy array
        tsum = np.sum(in_mat,1)
        return in_mat/tsum[:,None]

    def pfm2pwm(self):
        if not (self.pfm is None):
            if self.count_mat is None and np.any(self.pfm)==0:
                self.pfm = self._normalize(self.pfm+1e-6)
            # tmat = np.log2( (self.pfm) /self.background_distribution[None,:] )
            # kl_entropy_arr =  np.sum(self.pfm * tmat, 1)
            entropy_arr = np.sum( -1 * np.log2(self.pfm) * self.pfm, 1)
            info_content_arr = 2 - entropy_arr
            self.pwm = info_content_arr[:,None]*self.pfm
            return self.pwm
        else:
            raise ValueError('self.pfm is None') 
    
    def draw_logo(self):
        pwm = self.pwm
        out_file = self.out_logo_file

        fig = plt.figure()
        pwm_len = pwm.shape[0]
        fig.set_size_inches(pwm_len, 2)
        ax = fig.add_subplot(111)

        xshift = 0
        trans_offset = transforms.offset_copy(ax.transAxes,
                                          fig=fig,
                                          x=0,
                                          y=0,
                                          units='points')

        for k in range(pwm_len):
            yshift = 0
            tmparr = pwm[k,]
            tmpidx = np.argsort(tmparr)

            for bi in tmpidx:
                base = self.BASES[bi]
                base_height = tmparr[bi]

                txt = ax.text(0,
                              0,
                              base,
                              transform = trans_offset,
                            fontsize = 80,
                            color = self.COLOR_SCHEME[base],
                            weight = 'bold',
                            ha = 'center',
                            family = 'sans-serif',
                            )
                txt.set_clip_on(False)
                txt.set_path_effects([Logo.Scale(1.0, base_height)])
                fig.canvas.draw()
                window_ext = txt.get_window_extent(txt._renderer)
                yshift = window_ext.height*base_height*0.7  # times 0.7 to make letter more compact
                trans_offset = transforms.offset_copy(txt._transform, fig=fig, y=yshift, units='points')
            xshift += window_ext.width
            trans_offset = transforms.offset_copy(ax.transAxes, fig=fig, x=xshift, units='points')

        ax.set_yticks(range(0,3))
        ax.set_ylabel("bits")

        plt.tick_params(bottom = False, labelbottom = False)
        seaborn.despine(ax=ax, offset=30, trim=False, bottom=True)

        ax.set_yticklabels(np.arange(0,3,1))

        plt.savefig(out_file, dpi=200, bbox_inches='tight')
        plt.close()


if __name__=="__main__":

    counts = np.random.rand(6,4)*10
    logo = Logo(count_mat=counts)
    logo.draw_logo()
    
