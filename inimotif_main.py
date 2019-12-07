#!/usr/bin/env python3
import os    
import pickle
import matplotlib.pyplot as plt
from inimotif_core import KmerCounter, MotifManager
from yattag import Doc,indent

class FileProcessor:
    def __init__(self, file_name=None, file_type="fasta", out_dir=".",
              kmer_len=0, unique_kmer_in_seq_mode=True, revcom_flag=True, 
              consensus_seq=None, n_max_mutation=2, kmer_dict=None):
        assert os.path.exists(file_name), f"input file {file_name} does not exist"
        
        # store input parameters
        self.out_dir = out_dir

        self.file_name = file_name
        self.file_type = file_type

        self.kmer_len = kmer_len
        self.unique_kmer_in_seq_mode = unique_kmer_in_seq_mode
        self.revcom_flag = revcom_flag

        self.consensus_seq = consensus_seq
        self.n_max_mutation = n_max_mutation
        self.kmer_dict = kmer_dict

        # make output directory
        # preproc results, figures are stored in this directory
        self.mkdir(out_dir)

        # file names to be saved
        self.preproc_res_file = 'preproc.pickle'
        self.logo_forward_file = 'logo.forward.png'
        self.logo_revcom_file = 'logo.revcom.png'
        self.motif_posdis_file = 'posdis.png'
        self.kmer_hamdis_file = 'hamdis.png'
        self.motif_cooccur_dis_file = 'cooccurdis.png'

        # kmer counter and motif manager to be generated
        self.kmer_counter = None
        self.motif_manager = None
    
    def run(self):
        # output general information
        print(f'Start processing {self.file_name}, kmer_len={self.kmer_len}')

        # create kmer counts and motif manager
        self.kmer_counter = KmerCounter(self.kmer_len, unique_kmer_in_seq_mode=self.unique_kmer_in_seq_mode, revcom_flag=self.revcom_flag)
        self.kmer_counter.scan_file(self.file_name, file_type=self.file_type)
        print('kmer counter has scaned input file')
        
        self.motif_manager =  MotifManager(self.kmer_counter,self.consensus_seq, n_max_mutation=self.n_max_mutation, kmer_dict=self.kmer_dict, revcom_flag=self.revcom_flag)
        self.motif_manager.scan_file(self.file_name)
        print('motif manager has scaned input file')
        
        # make plots and save results
        with open( self.gen_absolute_path(self.preproc_res_file), 'wb') as f:
            pickle.dump(self, f)   # not sure if self could be pickled, TODO
        
        self.mk_plots()

    def mk_plots(self):
        kc = self.kmer_counter 
        mm = self.motif_manager

        kc.mk_kmer_dis_plot()
        self.save_figure( self.gen_absolute_path(self.kmer_hamdis_file) )
       
        mm.mk_logo_plot(mm.forward_motif_mat)
        self.save_figure( self.gen_absolute_path(self.logo_forward_file) )
       
        mm.mk_logo_plot(mm.revcom_motif_mat)
        self.save_figure( self.gen_absolute_path(self.logo_revcom_file) )
       
        mm.mk_motif_posdis_plot()
        self.save_figure( self.gen_absolute_path(self.motif_posdis_file) )
       
        mm.mk_bubble_plot()
        self.save_figure( self.gen_absolute_path(self.motif_cooccur_dis_file) )

    # generate html file string for displaying figures etc.
    def gen_html_str(self, img_dir):
        assert len(img_dir)>0, "image directory {img_dir} must be non-empty."
        if img_dir[-1]==os.sep:
            img_dir = img_dir[:-1]

        mm = self.motif_manager
        kc = self.kmer_counter
        img_files = [self.logo_forward_file, self.logo_revcom_file,
                     self.kmer_hamdis_file, self.motif_posdis_file, self.motif_cooccur_dis_file]
        doc, tag, text = Doc().tagtext()
        with tag('h2'):
            text(f'K={kc.k}')
        if mm.is_palindrome:
            tmpstr = 'palindrome'
        else:
            tmpstr = 'non-palindrome'
        with tag('p'):
            text(f'Consensus (forward) [{tmpstr}]: {mm.consensus_seq}')
        with tag('p'):
            text(f'Consensus (revcom) [{tmpstr}]: {kc.revcom(mm.consensus_seq)}')
        with tag('p'):
            text(f'Number of maximum allowed mutations: {mm.n_max_mutation}')
        with tag('p'):
            text(f'Total number of input sequences: {mm.n_seq}')
        with tag('p'):
            tmp_prec = round(mm.n_tfbs_forward_seq/mm.n_seq*100,2)
            text(f'Number of Forward motif Sequences: {mm.n_tfbs_forward_seq} ({tmp_prec}%) ')
        with tag('p'):
            tmp_prec = round(mm.n_tfbs_revcom_seq/mm.n_seq*100,2)
            text(f'Number of Revcom motif Sequences: {mm.n_tfbs_revcom_seq} ({tmp_prec}%)')
        with tag('p'):
            tmp_prec = round(mm.n_tfbs_seq/mm.n_seq*100,2)
            text(f'Number of motif (forward/revcom) Sequences: {mm.n_tfbs_seq} ({tmp_prec}%)')
        with tag('div'):
            for imgf in img_files:
                if imgf==self.kmer_hamdis_file:
                    doc.stag('img', klass="hamdis", src=img_dir+'/'+imgf, alt=imgf,  onclick=f"window.open('{img_dir}/{imgf}', '_blank');")
                else:
                    doc.stag('img', src=img_dir+'/'+imgf, alt=imgf,  onclick=f"window.open('{img_dir}/{imgf}', '_blank');")
        html_str = indent(doc.getvalue(), indent_text = True) # will also indent the text directly contained between <tag> and </tag>
        return html_str

    # make output directory if "outdir" does not exist
    def mkdir(self, outdir):
        if not os.path.exists(outdir):
            os.makedirs(outdir)
    
    def gen_absolute_path(self, filename):
        return os.path.join(self.out_dir,filename)
    
    # save figure, to do 
    def save_figure(self, file_name):
        plt.savefig(file_name)
        plt.close()


class ChipSeqProcessor:
    def __init__(self, file_name=None, file_type="fasta", identifier='out', out_dir=".",
              min_kmer_len=0, max_kmer_len=0, unique_kmer_in_seq_mode=True, revcom_flag=True, 
              consensus_seq=None, n_max_mutation=2, kmer_dict=None):
        assert len(out_dir)>0, "output directory must be non-empty string"
        if out_dir[-1]==os.sep:
            out_dir=out_dir[:-1]

        # store input parameters
        self.out_dir = out_dir
        self.identifier = identifier

        self.file_name = file_name
        self.file_type = file_type

        self.min_kmer_len = min_kmer_len
        self.max_kmer_len = max_kmer_len

        self.unique_kmer_in_seq_mode = unique_kmer_in_seq_mode
        self.revcom_flag = revcom_flag

        self.consensus_seq = consensus_seq
        self.n_max_mutation = n_max_mutation
        self.kmer_dict = kmer_dict

        # make output directory
        if not os.path.exists(out_dir):
            os.makedirs(out_dir)

    def run(self):
        html_div_list = []
        # run for different kmers
        for kmer_len in range(self.min_kmer_len, self.max_kmer_len+1):
            stem_dir = f'k{kmer_len}'
            out_dir = self.out_dir + os.sep + stem_dir
            fp = FileProcessor(file_name=self.file_name, file_type=self.file_type, out_dir=out_dir,
              kmer_len=kmer_len, unique_kmer_in_seq_mode=self.unique_kmer_in_seq_mode, revcom_flag=self.revcom_flag, 
              consensus_seq=self.consensus_seq, n_max_mutation=self.n_max_mutation, kmer_dict=self.kmer_dict)
            fp.run()
            html_div_list.append(fp.gen_html_str('./'+stem_dir))
        
        # generate html file
        style_str = self._get_style_str()
        doc, tag, text = Doc().tagtext()
        doc.asis('<!DOCTYPE html>')
        with tag('html',lang="en"):
            with tag('head'):
                with tag('title'):
                    text('Chip-Seq Result')
                doc.stag('meta', charset="utf-8")
                doc.stag('meta', name="viewport", content="width=device-width, initial-scale=1")
                with tag('style'):
                    text(style_str) # specify style
            with tag('body'):
                with tag('h1'):
                    text('IniMotif: Chip-seq Results')
                with tag('h2'):
                    text(f'identifier={self.identifier}')
                    doc.stag('br')
                    text(f'minimum kmer length: {self.min_kmer_len}')
                    doc.stag('br')
                    text(f'maximum kmer length: {self.max_kmer_len}')
                    doc.stag('br')
                for tmpstr in html_div_list:
                    doc.stag('hr')
                    doc.asis(tmpstr)
        
        # output html string
        html_str = indent(doc.getvalue(), indent_text = True) # will also indent the text directly contained between <tag> and </tag>
        outfile = self.out_dir + os.sep + self.identifier + '.html'
        with open(outfile,'w') as out_fh:
            out_fh.write(html_str)
        
    def _get_style_str(self):
        style_str = """
        body{
            font-family: sans-serif, Arial, Helvetica;
        }
        h2 {
            text-align: center;
        }
        div {
            text-align: justify;
        }
        .hamdis{
            display: block;
            margin-left: auto;
            margin-right: auto;
            width: 80%;
        }
        div img {
            display: inline-block;
            width: 48%;
        }
        div:after {
            content: '';
            display: inline-block;
            width: 100%;
        }
        """
        return style_str



class SelexSeqProcessor:
    pass

if __name__=="__main__":
    in_file = "/Users/lcheng/Documents/github/IniMotif-py/exampledata/NF1-1"
    csp = ChipSeqProcessor(file_name=in_file,identifier='NF',min_kmer_len=6, max_kmer_len=7,out_dir='/Users/lcheng/Documents/github/IniMotif/NF')
    csp.run()