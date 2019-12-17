#!/usr/bin/env python3
import os
import pickle
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from inimotif_core import KmerCounter, MotifManager, save_figure
from yattag import Doc,indent
import numpy as np

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
        #self.kmer_dict = {k: v for k, v in sorted(self.kmer_dict.items(), key=lambda item: item[1], reverse=True)}

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
            pickle.dump(self, f)   # current FileProcessor be pickled

        self.mk_plots()

    def mk_plots(self):
        kc = self.kmer_counter
        mm = self.motif_manager

        kc.mk_kmer_dis_plot(outfile=self.gen_absolute_path(self.kmer_hamdis_file))

        mm.mk_logo_plot(mm.forward_motif_mat, outfile=self.gen_absolute_path(self.logo_forward_file))

        mm.mk_logo_plot(mm.revcom_motif_mat, outfile=self.gen_absolute_path(self.logo_revcom_file))

        mm.mk_motif_posdis_plot(outfile=self.gen_absolute_path(self.motif_posdis_file))

        mm.mk_bubble_plot(outfile=self.gen_absolute_path(self.motif_cooccur_dis_file))

    # generate html file string for displaying figures etc.
    def gen_html_str(self, img_dir, title=None):
        assert len(img_dir)>0, "image directory {img_dir} must be non-empty."
        if img_dir[-1]==os.sep:
            img_dir = img_dir[:-1]

        mm = self.motif_manager
        kc = self.kmer_counter
        img_files = [self.logo_forward_file, self.logo_revcom_file,
                     self.kmer_hamdis_file, self.motif_posdis_file, self.motif_cooccur_dis_file]
        doc, tag, text = Doc().tagtext()
        with tag('h2'):
            if title:
                text(title)
            else:
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
        with tag('p'):
            text(f'Forward-Forward motif co-occurence index: {mm.ff_co_occur_index}')
        with tag('p'):
            text(f'Forward-RevCom motif co-occurence index: {mm.fr_co_occur_index}')

        with tag('div'):
            for imgf in img_files:
                if imgf==self.kmer_hamdis_file:
                    doc.stag('img', klass="hamdis", src=img_dir+'/'+imgf, alt=imgf,  onclick=f"window.open('{img_dir}/{imgf}', '_blank');")
                else:
                    doc.stag('img', src=img_dir+'/'+imgf, alt=imgf,  onclick=f"window.open('{img_dir}/{imgf}', '_blank');")
        html_str = indent(doc.getvalue(), indent_text = True) # will also indent the text directly contained between <tag> and </tag>
        return html_str

    # make output directory if "outdir" does not exist
    @staticmethod
    def mkdir(outdir):
        if not os.path.exists(outdir):
            os.makedirs(outdir)

    def gen_absolute_path(self, filename):
        return os.path.join(self.out_dir,filename)

    # load preprocessed data
    @staticmethod
    def load_pickle(in_file):
        with open(in_file,'rb') as file:
            fp = pickle.load(file)   # an FileProcessor object
            return fp

    @staticmethod
    def get_style_str():
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

    def gen_html(self, html_div_list):
        # generate html file
        style_str = FileProcessor.get_style_str()
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
        return html_str

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

        html_str = self.gen_html(html_div_list)
        outfile = self.out_dir + os.sep + self.identifier + '.html'
        with open(outfile,'w') as out_fh:
            out_fh.write(html_str)

class SelexSeqProcessor:
    def __init__(self, file_name_arr=None, file_type="fasta", identifier='out', out_dir=".",
              min_kmer_len=0, max_kmer_len=0, min_selex_round=0, max_selex_round=0,
              unique_kmer_in_seq_mode=True, revcom_flag=True, consensus_seq=None, n_max_mutation=2, kmer_dict=None):
        assert len(out_dir)>0, "output directory must be non-empty string"
        if out_dir[-1]==os.sep:
            out_dir=out_dir[:-1]

        # store input parameters
        self.out_dir = out_dir
        self.identifier = identifier

        self.file_name_arr = file_name_arr
        self.file_type = file_type

        self.min_kmer_len = min_kmer_len
        self.max_kmer_len = max_kmer_len

        self.min_selex_round = min_selex_round
        self.max_selex_round = max_selex_round

        self.unique_kmer_in_seq_mode = unique_kmer_in_seq_mode
        self.revcom_flag = revcom_flag

        self.consensus_seq = consensus_seq
        self.n_max_mutation = n_max_mutation
        self.kmer_dict = kmer_dict

        self.trend_figure_dir = 'trend_figure'

        # make output directory
        FileProcessor.mkdir(out_dir)

        # make trend figure directory
        FileProcessor.mkdir(self.out_dir + os.sep + self.trend_figure_dir)

    def run(self):
        html_div_k_list = [[] for _ in range(self.max_kmer_len+1)]
        html_div_r_list = [[] for _ in range(self.max_selex_round+1)]
        # run for different kmers
        for kmer_len in range(self.min_kmer_len, self.max_kmer_len+1):
            selex_res = []
            for i_round,file_name in zip(range(self.min_selex_round, self.max_selex_round+1),self.file_name_arr):
                stem_dir = f'r{i_round}k{kmer_len}'
                out_dir = self.out_dir + os.sep + stem_dir
                fp = FileProcessor(file_name=file_name, file_type=self.file_type, out_dir=out_dir,
                    kmer_len=kmer_len, unique_kmer_in_seq_mode=self.unique_kmer_in_seq_mode, revcom_flag=self.revcom_flag,
                    consensus_seq=self.consensus_seq, n_max_mutation=self.n_max_mutation, kmer_dict=self.kmer_dict)
                fp.run()
                html_div_k_list[kmer_len].append(fp.gen_html_str('./'+stem_dir, title=f'Round={i_round} K={kmer_len}'))
                selex_res.append(fp)
                
            # generate kmer trend figures
            trend_fig_file = self.out_dir + os.sep + self.trend_figure_dir + os.sep + f'k{kmer_len}.png'
            self.mk_kmer_trend_fig(selex_res, trend_fig_file)
            
        for kmer_len in range(self.min_kmer_len, self.max_kmer_len+1):
            k_list = html_div_k_list[kmer_len]
            for i_round,div in zip(range(self.min_selex_round, self.max_selex_round+1), k_list):
                html_div_r_list[i_round].append(div)

        # generate html for each round
        for i_round in range(self.min_selex_round, self.max_selex_round+1):
            html_str = self.gen_html_round(html_div_r_list[i_round], i_round)
            outfile = self.out_dir + os.sep + self.identifier + f'_round_{i_round}.html'
            with open(outfile,'w') as out_fh:
                out_fh.write(html_str)

        # generate html for each kmer_len
        for kmer_len in range(self.min_kmer_len, self.max_kmer_len+1):
            html_str = self.gen_html_k(html_div_k_list[kmer_len], kmer_len, f'./{self.trend_figure_dir}', f'k{kmer_len}.png')
            outfile = self.out_dir + os.sep + self.identifier + f'_k_{kmer_len}.html'
            with open(outfile,'w') as out_fh:
                out_fh.write(html_str)

    # make kmer trend figure
    def mk_kmer_trend_fig(self, selex_round_res_list, outfile="selex_trend.png"):
        n_round = self.max_selex_round - self.min_selex_round + 1
        kc = selex_round_res_list[n_round-1].kmer_counter
        kmer_len = kc.k

        # random sample kmers to be displayed, top kmers are always included
        # draw large amount of lines is slow 
        top_kh_arr = [x for x in kc.top_kmers_list[0]]
        all_kh_arr = list(kc.kmer_dict.keys())
        sub_kh_arr = top_kh_arr.copy()
        n_disp_sample = 800  # number of kmers to be displayed
        if n_disp_sample<len(all_kh_arr):
            tmparr = np.random.choice(all_kh_arr,n_disp_sample)
            sub_kh_arr += [kh for kh in tmparr]  # tolist change element type to int
        else:
            sub_kh_arr += all_kh_arr
        
        # only keep one kmer of a pair (forward / revcom)
        # top kmers are in the front, so will be kept
        def filter_kmer(kc, kh_arr):
            res = []
            for kh in kh_arr:
                if kh in res or kc.revcom_hash(kh) in res:
                    continue
                else:
                    res.append(kh)
            return res

        sub_kh_arr = filter_kmer(kc, sub_kh_arr)
        n_disp_sample = len(sub_kh_arr)

        sub_kh_cnt_mat = np.zeros((n_disp_sample,n_round),dtype="float")
        for r,res in enumerate(selex_round_res_list):
            for i,kh in enumerate(sub_kh_arr):
                sub_kh_cnt_mat[i,r] = res.kmer_counter.get_pair_cnt(kh)
        
        n_total_kmer_arr = np.array([res.kmer_counter.n_total_kmer for res in selex_round_res_list])
        sub_kh_freq_mat = sub_kh_cnt_mat/n_total_kmer_arr[None, :]
        sub_kh_log_freq_mat = np.log10( (sub_kh_freq_mat+1e-9)/(1-sub_kh_freq_mat+1e-9) )

        x_round = np.arange(self.min_selex_round, self.max_selex_round+1)

        fig = plt.figure(figsize=(10,10))
        grid = plt.GridSpec(2, 3, wspace=0.4, hspace=0.3)

        top = fig.add_subplot(grid[:-1,:])
        top.set_xlabel("SELEX round")
        top.set_ylabel("log10(f/(1-f))")
        top.set_title(f"log10 {kmer_len}-mer frequency trend")
        top.set_xlim([self.min_selex_round-1, self.max_selex_round+2])
        top.set_xticks(np.linspace( self.min_selex_round-1, self.max_selex_round, num=n_round+2, endpoint=True))
        top.spines['right'].set_visible(False)
        top.spines['top'].set_visible(False)

        bottom = fig.add_subplot(grid[-1,:-1])
        bottom.set_xlabel("SELEX round")
        bottom.set_ylabel("f = #kmer/#total_kmer")
        bottom.set_title(f"{kmer_len}-mer frequency trend")
        bottom.set_xlim([self.min_selex_round-1, self.max_selex_round+2])
        bottom.set_xticks(np.linspace(self.min_selex_round-1, self.max_selex_round, num=n_round+2, endpoint=True))
        bottom.spines['right'].set_visible(False)
        bottom.spines['top'].set_visible(False)

        bar = fig.add_subplot(grid[-1,-1:])
        bar.set_xlabel("SELEX round")
        bar.set_ylabel("Total kmers")
        bar.set_title(f"#total {kmer_len}-mers")
        bar.set_xticks(x_round)
        bar.set_xlim(self.min_selex_round-1, self.max_selex_round+1)

        colourslist = ['C0', 'C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'C7']

        # plot randomly sampled kmers
        for i in range(n_disp_sample):
            top.plot(x_round, sub_kh_log_freq_mat[i,], color = '0.75', linestyle='--', linewidth=0.5, marker="x", markevery=None, alpha=0.5, zorder=0)
            bottom.plot(x_round, sub_kh_freq_mat[i,], color = '0.75', linestyle='--', linewidth=0.5, marker="x", markevery=None, alpha=0.5, zorder=0)
        
        # plot lines for top kmers
        n_top_kmer = len(top_kh_arr)
        for i in range(n_top_kmer):
            ind = n_top_kmer-1-i
            top.plot(x_round, sub_kh_log_freq_mat[ind,],
                    color=colourslist[ind], linewidth=2, marker="s", markevery=None, zorder=(i+1)*2)
            
            bottom.plot(x_round, sub_kh_freq_mat[ind,],
                    color=colourslist[ind], linewidth=2, marker="s", markevery=None, zorder=(i+1)*2)
        
        # plot annotations
        ymint, ymaxt = top.get_ylim()
        ypost = np.linspace(ymint, ymaxt, num=20, endpoint=True)
        yminb, ymaxb = bottom.get_ylim()
        yposb = np.linspace(yminb, ymaxb, num=20, endpoint=True)
        for i in range(n_top_kmer):
            ind = n_top_kmer-1-i
            kh = top_kh_arr[ind]
            tmptext = f'{ind}. {kc.hash2kmer(kh)} / {kc.hash2kmer( kc.revcom_hash(kh) )}'
            top.annotate(tmptext, 
                         (x_round[-1],sub_kh_log_freq_mat[ind,][-1]),
                         (x_round[-1]+0.2, ypost[-(ind+2)] ),
                         size=10, fontname='monospace', weight='bold', 
                         arrowprops=dict(color=colourslist[ind], shrink=0.05, width=0.05, headwidth=0.4), 
                         color=colourslist[ind]
                          )
            bottom.annotate(tmptext, 
                         (x_round[-1],sub_kh_freq_mat[ind,][-1]),
                         (x_round[-1]+0.2, yposb[-(ind+2)] ),
                         size=10, fontname='monospace', weight='bold', 
                         arrowprops=dict(color=colourslist[ind], shrink=0.05, width=0.05, headwidth=0.4), 
                         color=colourslist[ind]
                          )
        
        # plot kmer counts in different rounds    
        bar.bar(x_round, n_total_kmer_arr)

        save_figure(outfile)
        

    # generate html for kmer_len=k
    def gen_html_k(self, html_div_list, kmer_len, trend_fig_dir, trend_fig_name):
        # generate html file
        style_str = FileProcessor.get_style_str()
        doc, tag, text = Doc().tagtext()
        doc.asis('<!DOCTYPE html>')
        with tag('html',lang="en"):
            with tag('head'):
                with tag('title'):
                    text('SELEX-Seq Result')
                doc.stag('meta', charset="utf-8")
                doc.stag('meta', name="viewport", content="width=device-width, initial-scale=1")
                with tag('style'):
                    text(style_str) # specify style
            with tag('body'):
                with tag('h1'):
                    text(f'IniMotif: SELEX-seq Results. kmer_len={kmer_len}')
                with tag('h2'):
                    text(f'identifier={self.identifier}')
                    doc.stag('br')
                    text(f'minimum round number: {self.min_selex_round}')
                    doc.stag('br')
                    text(f'maximum round number: {self.max_selex_round}')
                    doc.stag('br')
                for tmpstr in html_div_list:
                    doc.stag('hr')
                    doc.asis(tmpstr)
                # add trend figure
                doc.stag('hr')
                with tag('h2'):
                    text('SELEX kmer trend figure')
                with tag('div'):
                    doc.stag('img', klass="hamdis", src=trend_fig_dir+'/'+trend_fig_name, alt=trend_fig_name,  onclick=f"window.open('{trend_fig_dir}/{trend_fig_name}', '_blank');")

        # output html string
        html_str = indent(doc.getvalue(), indent_text = True) # will also indent the text directly contained between <tag> and </tag>
        return html_str

    def gen_html_round(self, html_div_list, i_round):
        # generate html file
        style_str = FileProcessor.get_style_str()
        doc, tag, text = Doc().tagtext()
        doc.asis('<!DOCTYPE html>')
        with tag('html',lang="en"):
            with tag('head'):
                with tag('title'):
                    text(f'SELEX-Seq Result, Round {i_round}')
                doc.stag('meta', charset="utf-8")
                doc.stag('meta', name="viewport", content="width=device-width, initial-scale=1")
                with tag('style'):
                    text(style_str) # specify style
            with tag('body'):
                with tag('h1'):
                    text(f'IniMotif: SELEX-Seq Result, Round {i_round}')
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
        return html_str


if __name__=="__main__":
    # in_file = "/Users/lcheng/Documents/github/IniMotif-py/exampledata/NF1-1"
    # csp = ChipSeqProcessor(file_name=in_file,identifier='NF',min_kmer_len=6, max_kmer_len=7,out_dir='/Users/lcheng/Documents/github/IniMotif/NF')
    # csp.run()

    # load preprocessed result
    #in_file = "/Users/lcheng/Documents/github/IniMotif/NF/k6/preproc.pickle"
    #fp = FileProcessor.load_pickle(in_file)
    #fp.kmer_counter.disp_kmer_info()

    # in_dir = "/home/alex/Desktop/IniMotif/Data/"
    # file_list = ['NF1-1t','NF1-2t','NF1-3t','NF1-4t']
    # file_name_arr = [in_dir+f for f in file_list]
    # ssp = SelexSeqProcessor(file_name_arr=file_name_arr,identifier='NF',min_kmer_len=6, max_kmer_len=7,min_selex_round=1,max_selex_round=4,out_dir='/home/alex/Desktop/IniMotif/IniMotif-master/NF1')
    # ssp.run()

    in_dir = "/Users/lcheng/Documents/github/IniMotif-py/exampledata/"
    file_list = ['NF1-1','NF1-2','NF1-3','NF1-4']
    file_name_arr = [in_dir+f for f in file_list]
    out_dir='./test'
    ssp = SelexSeqProcessor(file_name_arr=file_name_arr,identifier='NF',min_kmer_len=18, max_kmer_len=20,min_selex_round=1,max_selex_round=4,out_dir=out_dir)
    ssp.run()
