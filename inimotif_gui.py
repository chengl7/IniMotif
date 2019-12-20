#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  2 22:35:43 2019

@author: lcheng
"""
from tkinter import filedialog
from tkinter import messagebox
from tkinter import Tk, END, HORIZONTAL, PhotoImage, N, W, E, S, Text
from tkinter import StringVar, IntVar, BooleanVar
from tkinter.ttk import Button, Label, Entry, Frame, Notebook, Progressbar, Style, Combobox, Radiobutton, Checkbutton, Scrollbar
import os
from inimotif_main import FileProcessor,ChipSeqProcessor,SelexSeqProcessor
from inimotif_util import Masker, MotifScanner
#from tkinter import *

# import time
import threading
import pickle
import sys

class Command:
    def __init__(self, callback, *args, **kwargs):
        self.callback = callback
        self.args = args
        self.kwargs = kwargs

    def __call__(self):
        return self.callback(*self.args, **self.kwargs)
    
class NumGenerator:
    def __init__(self, num=0):
        self.num = num-1
    
    def get(self):
        self.num += 1
        return self.num

def set_style():
    app_style = Style()
    app_style.configure('main.TLabel', font=('Helvetica', 25))
    app_style.configure('main.TButton', font=('Helvetica', 25, 'bold'), padding=20)
    app_style.configure('grey.TEntry', fg="light grey")

def grid_widgets_line(master, widgets_list, irow, icol=None):
    if icol is None:
        icol=0
    for i,w in enumerate(widgets_list):
        w.grid(row=irow, column=(icol+i) )

def gen_multi_labels(master, label_text_list, grid_on=False, irow=None, icol=None):
    label_list = []
    for ltext in label_text_list:
        label_list.append( Label(master, text=ltext) )
    if grid_on:
        grid_widgets_line(master, label_list, irow, icol)
    return label_list

def gen_label_entry(master, label_text, entry_text, grid_on=False, irow=None, icol=None):
    label = Label(master, text=label_text)
    entry = Entry(master, style='grey.TEntry')
    entry.insert(END, entry_text)
    if grid_on:
        grid_widgets_line(master, [label, entry], irow, icol)
    return label, entry

def gen_file_entry(master, label_text, entry_text, button_text, grid_on=False, irow=None, icol=None):
    def enter_filename():
        file = filedialog.askopenfilename(initialdir='.',title = "Select file", filetypes = (("fasta files",".fasta .fa .gz"),("all files","*.*")))
        entry.delete(0,END)
        entry.insert(0,file)
    label = Label(master, text=label_text)
    entry = Entry(master, style='grey.TEntry')
    entry.insert(END, entry_text)
    
    button = Button(master, text=button_text, command=enter_filename)
    
    if grid_on:
        grid_widgets_line(master, [label, entry, button], irow, icol)
        
    return label, entry, button

def gen_pickel_file_entry(master, label_text, button_text, grid_on=False, irow=None, icol=None):
    def enter_filename():
        file = filedialog.askopenfilename(initialdir='.',title = "Select file", filetypes = (("preproc files",".pickle"),("all files","*.*")))
        entry.delete(0,END)
        entry.insert(0,file)
    label = Label(master, text=label_text)
    entry = Entry(master, style='grey.TEntry')
    entry.insert(END, 'output_directory/preproc.pickle')
    button = Button(master, text=button_text, command=enter_filename)
    if grid_on:
        grid_widgets_line(master, [label, entry, button], irow, icol)
    return label, entry, button

def gen_directory_entry(master, label_text, entry_text, button_text, grid_on=False, irow=None, icol=None):
    def enter_outdir():
        outdir = filedialog.askdirectory(initialdir='.')
        entry.delete(0,END)
        entry.insert(0,outdir)
    label = Label(master, text=label_text)
    entry = Entry(master, style='grey.TEntry')
    entry.insert(END, entry_text)
    
    button = Button(master, text=button_text, command=enter_outdir)
    
    if grid_on:
        grid_widgets_line(master, [label, entry, button], irow, icol)
    return label, entry, button

def gen_mask_pattern(master, label_text, type_value, revcom_value, grid_on=False, irow=None, icol=None):
    # type_value: StringVar, type of pattern, 'Repeat' or 'Motif'
    # revcom_value: IntVar, if revcom of pattern should be used
    option_list = ['Repeat', 'Motif']
    label = Label(master, text=label_text)
    combobox = Combobox(master, textvariable = type_value, values=option_list)
    check = Checkbutton(master, text='', variable=revcom_value, onvalue=1, offvalue=0)
    entry1 = Entry(master, style='grey.TEntry')
    # entry1.insert(END, entry1_text)
    entry2 = Entry(master, style='grey.TEntry')
    # entry2.insert(END, entry2_text)
    
    if grid_on:
        grid_widgets_line(master, [label, combobox, check, entry1, entry2], irow, icol)
    
    return label, combobox, check, entry1, entry2

def gen_motif_row(master, label_text, revcom_value, grid_on=False, irow=None, icol=None):
    # type_value: StringVar, type of pattern, 'Repeat' or 'Motif'
    # revcom_value: IntVar, if revcom of pattern should be used
    label = Label(master, text=label_text)
    check = Checkbutton(master, text='', variable=revcom_value, onvalue=1, offvalue=0)
    
    entry1 = Entry(master, style='grey.TEntry')
    # entry1.insert(END, entry1_text)
    
    entry2 = Entry(master, style='grey.TEntry')
    # entry2.insert(END, entry2_text)
    
    if grid_on:
        grid_widgets_line(master, [label, check, entry1, entry2], irow, icol)
    
    return label, check, entry1, entry2

class IORedirector(object):
    '''A general class for redirecting I/O to this Text widget.'''
    def __init__(self,text_area):
        self.text_area = text_area

class StdoutRedirector(IORedirector):
    '''A class for redirecting stdout to this Text widget.'''
    def write(self,str):
        self.text_area.write(str,False)

class Application(Frame):
    def __init__(self, master=None):
        super().__init__(master, borderwidth=5, width=0.8*screen_width, height=0.7*screen_height)
        self.master = master
        self.proc = None
        
        self.main_widgets = []
        
        script_dir = os.path.dirname(os.path.realpath(__file__))
        sep = os.sep
        self.logoimage = PhotoImage(file=f'{script_dir}{sep}GUIgraphics{sep}logo.png')
        self.maskerimage = PhotoImage(file=f'{script_dir}{sep}GUIgraphics{sep}masker.png')

        self.run_func = None
        self.run_button = None
        self.progress = None
        self.progress_row = None
        
        set_style()

        self.grid(row=0, column=0,sticky=N+S+W+E)
        # self.grid_propagate(False)

        self.select_analysis()

        self.chip_para_dict = {}
        self.selex_para_dict = {}
        self.mask_para = []
        self.motif_scan_para = []

        self.selex_gui_dict = {}
        
        self.master.grid_columnconfigure(0, weight=1)
        self.master.grid_rowconfigure(0, weight=1)

    # clear all widgets in the main panel
    def clear(self):
        # destroy widgets, use pack_forget() or grid_forget() to hide
        for w in self.main_widgets:
            w.destroy()

    def select_analysis(self):
        row_gen = NumGenerator()
        
        logo = Label(self, image=self.logoimage)
        logo.grid(row=row_gen.get())
        
        # l1 = Label(self, text="Choose Analysis Type", style="main.TLabel")
        l1 = Label(self, text="", style="main.TLabel")
        l1.grid(row=row_gen.get())

        chip_btn = Button(self, text="Chip-seq", style="main.TButton", command=self.init_chipseq_gui)
        chip_btn.grid(row=row_gen.get(), pady=10)

        selex_btn = Button(self, text="Selex-seq", style="main.TButton", command=self.init_selexseq_gui)
        selex_btn.grid(row=row_gen.get(), pady=10)
        
        masker_btn = Button(self, text="Masker", style="main.TButton", command=self.init_masker_gui)
        masker_btn.grid(row=row_gen.get(), pady=10)

        l2 = Label(self, text="")
        l2.grid(row=row_gen.get(), padx=40, pady=10)

        quit_btn = Button(self, text="QUIT", command=self.master.destroy)
        quit_btn.grid(row=row_gen.get(), pady=10, ipadx=10, ipady=5)
        
        # https://stackoverflow.com/questions/45847313/what-does-weight-do-in-tkinter
        self.master.grid_columnconfigure(0, weight=1)    
        # self.columnconfigure(1, weight=2)
        # self.columnconfigure(2, weight=1)
        # self.rowconfigure(0, weight=1)
        self.main_widgets += [logo, l1, chip_btn, selex_btn, masker_btn, l2, quit_btn]
    
    def init_masker_gui(self):
        self.clear()
        row_gen = NumGenerator()

        logo = Label(self, image=self.maskerimage)
        # logo.grid_columnconfigure(0, weight=1)
        logo.grid(row=row_gen.get(),columnspan=6)
        
        irow = row_gen.get()
        infile_label, infile_entry, infile_button = gen_file_entry(self, 'Input file', './inputfile.fasta', 'Open File')
        infile_label.grid(row = irow, column=0, columnspan=1, sticky='w')
        infile_entry.grid(row = irow, column=1, columnspan=1, sticky='w')
        infile_button.grid(row = irow, column=2, columnspan=1, sticky='w')
        
        irow = row_gen.get()
        outfile_label, outfile_entry, outfile_button = gen_file_entry(self, 'Output file', './output_file.fasta', 'Open File')
        outfile_label.grid(row = irow, column=0, columnspan=1, sticky='w')
        outfile_entry.grid(row = irow, column=1, columnspan=1, sticky='w')
        outfile_button.grid(row = irow, column=2, columnspan=1, sticky='w')
        
        irow = row_gen.get()
        text_list = ['Repeat/\nMotif', 'Type', 'Include\nRevCom\nPattern', 'Sequence', 'n_min_repeats/\nn_max_mutation']
        for i,text in enumerate(text_list):
            tmplabel = Label(self, text=text)
            tmplabel.grid(row=irow, column=i, pady=20, sticky='w')
                
        def add_pattern(master, label_text, type_value, revcom_value, irow, seq_entry_arr, num_entry_arr):
            if len(seq_entry_arr)>=n_max_pat:
                return
            # label, combobox, checkbox, entry1, entry2
            widget_list = gen_mask_pattern(master, label_text, type_value, revcom_value)
            for i,w in enumerate(widget_list):
                w.grid(row=irow, column=i, sticky='w')
            
            revcom_checkbox, seq_entry, num_entry = widget_list[2:]
                
            revcom_checkbox.invoke()
            # print(f'irow={irow}')
            seq_entry_arr.append(seq_entry)
            num_entry_arr.append(num_entry)
            
            shift_buttons(irow)
            
        def shift_buttons(irow):
            self.progress_row += 1
            add_entry_button.grid_forget()
            add_entry_button['command']=callback_fun_arr[func_ind_gen.get()]
            add_entry_button.grid(row=irow, column=5, sticky='w')
            self.run_button.grid_forget()
            self.run_button.grid(row=irow+2, column=2, pady=50, sticky='w')
            quit_button.grid_forget()
            quit_button.grid(row=irow+2, column=4, pady=50, sticky='w')
        
        def run_analysis():
            self.mask_para = []
            
            n_input_pat = len(seq_entry_arr)
            masker = Masker()
            for i in range(n_input_pat):
                seq = seq_entry_arr[i].get()
                num = int(num_entry_arr[i].get())
                revcom_flag = revcom_arr[i].get()
                if pat_type_arr[i].get()=='Repeat':
                    masker.add_reppat(seq, num, revcom_flag)
                else:
                    masker.add_motif(seq, num, revcom_flag)
            
            infile = infile_entry.get()
            outfile = outfile_entry.get()
            
            self.mask_para.append(masker)
            self.mask_para.append(infile)
            self.mask_para.append(outfile)
            
            self.run_mask()
        
        n_max_pat = 20
        pat_type_arr = [StringVar() for _ in range(n_max_pat)]
        revcom_arr = [BooleanVar() for _ in range(n_max_pat)]
        seq_entry_arr = []
        num_entry_arr = []
        func_ind_gen = NumGenerator()
        callback_fun_arr = [Command(add_pattern, self, f'Pattern {i+1}', pat_type_arr[i], revcom_arr[i], irow+i+1, seq_entry_arr, num_entry_arr)
                            for i in range(n_max_pat) ]
        # for cmd in callback_fun_arr:
        #     print(cmd.args)
        add_entry_button = Button(self, text="Add Pattern", command=callback_fun_arr[func_ind_gen.get()])
        add_entry_button.grid(row=irow, column=5, sticky='w')

        self.progress_row = row_gen.get()
        irow = row_gen.get()
        self.run_button = Button(self, text="RUN", command=run_analysis)
        self.run_button.grid(row=irow, column=2, pady=50, sticky='w')
        
        quit_button = Button(self, text="QUIT", command=self.master.destroy)
        quit_button.grid(row=irow, column=4, pady=50, sticky='w')

        self.progress = Progressbar(self, orient=HORIZONTAL, length=300,  mode='indeterminate')
        
        self.columnconfigure(0, weight=1)
        self.columnconfigure(1, weight=1)
        self.columnconfigure(2, weight=1)
        self.columnconfigure(3, weight=1)
        self.columnconfigure(4, weight=1)
        self.columnconfigure(5, weight=1)
    
    def run_mask(self):
        def run_mask_task(*args):
            masker = args[0]
            infile = args[1]
            outfile = args[2]
            self.progress.grid(row=self.progress_row ,column=0, columnspan=5, pady=30)
            self.progress.start()
            masker.mask_file(infile, outfile)
            self.progress.stop()
            self.progress.grid_forget()

            messagebox.showinfo('Info', "Process completed!")
            self.run_button['state']='normal'

        self.run_button['state']='disabled'

        threading.Thread(target=run_mask_task, args=self.mask_para).start()

    def init_chipseq_gui(self):
        self.clear()
        
        logo = Label(self, image=self.logoimage)
        logo.grid(row=0, column=0)

        note = Notebook(self)
        main = Frame(note)
        motif_scan = Frame(note)
        topkmer = Frame(note)
        kmer_query = Frame(note)
        
        note.add(main, text = "Chip-Seq Main")
        note.add(motif_scan, text = "Motif Scan")
        note.add(topkmer, text='Top kmer')
        note.add(kmer_query, text='Kmer query')

        self.init_chipseq_gui_tabmain(main)
        self.init_motifscan_guitab(motif_scan)
        self.init_topkmer_guitab(topkmer)
        self.init_kmerquery_guitab(kmer_query)

        note.grid(row=1,column=0)
        
        self.master.grid_columnconfigure(0, weight=1)
        self.master.grid_rowconfigure(0, weight=1)
        self.master.grid_rowconfigure(1, weight=5)        
    
    # display top kmer information
    def init_topkmer_guitab(self, master):
        row_gen = NumGenerator()
        
        irow = row_gen.get()
        infile_label, infile_entry, infile_button = gen_pickel_file_entry(master, 'Result File', 'Open file')
        infile_label.grid(row=irow, column=0, sticky='w')
        infile_entry.grid(row=irow, column=1, sticky='w')
        infile_button.grid(row=irow, column=2, sticky='w')
        
        irow = row_gen.get()
        n_topkmer_label, n_topkmer_entry = gen_label_entry(master, 'Number of top kmer', '6')
        n_topkmer_label.grid(row=irow, column=0, sticky='w')
        n_topkmer_entry.grid(row=irow, column=1, sticky='w')
        
        # create a Text widget
        irow = row_gen.get()
        txt = Text(master,width=32, height=10)
        txt.grid(row=irow, column=0, columnspan=3, sticky="nsew", padx=2, pady=2)

        # create a Scrollbar and associate it with txt
        scrollb = Scrollbar(master, command=txt.yview)
        scrollb.grid(row=irow, column=3, sticky='nsew')
        txt['yscrollcommand'] = scrollb.set
        
        def run_analysis():
            file = infile_entry.get()
            with open(file, 'rb') as f:
                fp = pickle.load(f)
            kc = fp.kmer_counter
            res = kc.get_top_kmers( int(n_topkmer_entry.get()) )
            str_list = kc.disp_kmer_info(kmer_list=res[0])
            txt.delete('1.0', END)
            txt.insert(END, "\n".join(str_list) )
        
        irow = row_gen.get()
        self.run_button = Button(master, text="RUN", command=run_analysis)
        self.run_button.grid(row=irow, column=1, pady=20, sticky='w')
        
        quit_button = Button(master, text="QUIT", command=self.master.destroy)
        quit_button.grid(row=irow, column=2, pady=20, sticky='w')
        
        master.columnconfigure(0, weight=5)
        master.columnconfigure(1, weight=10)
        master.columnconfigure(2, weight=5)
        master.columnconfigure(3, weight=1)
        master.rowconfigure(0, weight=1)
        master.rowconfigure(1, weight=1)
        master.rowconfigure(2, weight=7)
        master.rowconfigure(3, weight=1) 
        
    # query kmer information in a list
    def init_kmerquery_guitab(self, master):
        row_gen = NumGenerator()
        
        irow = row_gen.get()
        infile_label, infile_entry, infile_button = gen_pickel_file_entry(master, 'Result File', 'Open file')
        infile_label.grid(row=irow, column=0, sticky='w')
        infile_entry.grid(row=irow, column=1, sticky='w')
        infile_button.grid(row=irow, column=2, sticky='w')
        
        irow = row_gen.get()
        kmer_info_label = Label(master, text='Input kmer list')
        kmer_info_label.grid(row=irow, column=0, sticky='w')
        
        irow = row_gen.get()
        kmer_input_txt = Text(master, width=32, height=10)
        kmer_input_txt.grid(row=irow, column=0, columnspan=3, sticky="nsew", padx=2, pady=2)
        kmer_input_txt.insert(END, 'Input kmers here (one kmer per line)')
        
        # create a Scrollbar and associate it with txt
        kmer_scrollb = Scrollbar(master, command=kmer_input_txt.yview)
        kmer_scrollb.grid(row=irow, column=3, sticky='nsew')
        kmer_input_txt['yscrollcommand'] = kmer_scrollb.set
        
        irow = row_gen.get()
        output_info_label = Label(master, text='Output')
        output_info_label.grid(row=irow, column=0, sticky='w')
        
        # create a Text widget
        irow = row_gen.get()
        txt = Text(master, width=32, height=10)
        txt.grid(row=irow, column=0, columnspan=3, sticky="nsew", padx=2, pady=2)

        # create a Scrollbar and associate it with txt
        scrollb = Scrollbar(master, command=txt.yview)
        scrollb.grid(row=irow, column=3, sticky='nsew')
        txt['yscrollcommand'] = scrollb.set
        
        def run_analysis():
            file = infile_entry.get()
            with open(file, 'rb') as f:
                fp = pickle.load(f)
            kc = fp.kmer_counter
            
            tmpstr = kmer_input_txt.get("1.0", "end-1c")
            kmer_list = tmpstr.split("\n")
            str_list = []
            for i,kmer in enumerate(kmer_list):
                tmpstr = f'line {i}, '
                if len(kmer)!=kc.k:
                    tmpstr += f'{kmer}, kmer length must be {kc.k}. Return!'
                    txt.delete('1.0', END)
                    txt.insert(END, tmpstr )
                    return
            
            str_list = kc.disp_kmer_info(kmer_list=kmer_list)
            txt.delete('1.0', END)
            txt.insert(END, "\n".join(str_list) )
        
        irow = row_gen.get()
        self.run_button = Button(master, text="RUN", command=run_analysis)
        self.run_button.grid(row=irow, column=1, pady=20, sticky='w')
        
        quit_button = Button(master, text="QUIT", command=self.master.destroy)
        quit_button.grid(row=irow, column=2, pady=20, sticky='w')
        
        master.columnconfigure(0, weight=5)
        master.columnconfigure(1, weight=10)
        master.columnconfigure(2, weight=5)
        master.columnconfigure(3, weight=1)
        master.rowconfigure(0, weight=1)
        master.rowconfigure(1, weight=1)
        master.rowconfigure(2, weight=3)
        master.rowconfigure(3, weight=1)
        master.rowconfigure(4, weight=3)
        master.rowconfigure(5, weight=1)
    
    def run_scan(self):
        def run_scan_task(*args):
            scanner = args[0]
            infile = args[1]
            outfile = args[2]
            self.progress.grid(row=self.progress_row ,column=0, columnspan=3, pady=30)
            self.progress.start()
            scanner.scan_file(infile, outfile)
            self.progress.stop()
            self.progress.grid_forget()

            messagebox.showinfo('Info', "Process completed!")
            self.run_button['state']='normal'

        self.run_button['state']='disabled'

        threading.Thread(target=run_scan_task, args=self.motif_scan_para).start()
    
    def init_motifscan_guitab(self, master):
        row_gen = NumGenerator()
        
        irow = row_gen.get()
        infile_label, infile_entry, infile_button = gen_file_entry(master, 'Input file', './inputfile.fasta', 'Open File')
        infile_label.grid(row = irow, column=0, columnspan=1, sticky='w')
        infile_entry.grid(row = irow, column=1, columnspan=1, sticky='w')
        infile_button.grid(row = irow, column=2, columnspan=1, sticky='w')
        
        irow = row_gen.get()
        outfile_label, outfile_entry, outfile_button = gen_file_entry(master, 'Output file', './output_file.html', 'Open File')
        outfile_label.grid(row = irow, column=0, columnspan=1, sticky='w')
        outfile_entry.grid(row = irow, column=1, columnspan=1, sticky='w')
        outfile_button.grid(row = irow, column=2, columnspan=1, sticky='w')
        
        irow = row_gen.get()
        text_list = ['Motif', 'Include\nRevCom', 'Consensus Sequence', 'n_max_mutation']
        for i,text in enumerate(text_list):
            tmplabel = Label(master, text=text)
            tmplabel.grid(row=irow, column=i, pady=20, sticky='w')
        
        def add_pattern(master, label_text, revcom_value, irow, seq_entry_arr, num_entry_arr):
            # allow no more than 2 motifs
            if len(seq_entry_arr)>=2:
                return
            # label, checkbox, entry1, entry2
            widget_list = gen_motif_row(master, label_text, revcom_value)
            for i,w in enumerate(widget_list):
                w.grid(row=irow, column=i, sticky='w')
            
            revcom_checkbox, seq_entry, num_entry = widget_list[1:]
                
            revcom_checkbox.invoke()
            seq_entry_arr.append(seq_entry)
            num_entry_arr.append(num_entry)
            
            shift_buttons(irow)
            
        def shift_buttons(irow):
            self.progress_row += 1
            add_entry_button.grid_forget()
            add_entry_button['command']=callback_fun_arr[func_ind_gen.get()]
            add_entry_button.grid(row=irow, column=4, sticky='w')
            run_button.grid_forget()
            run_button.grid(row=irow+2, column=2, pady=50, sticky='w')
            quit_button.grid_forget()
            quit_button.grid(row=irow+2, column=4, pady=50, sticky='w')
        
        def run_analysis():
            self.motif_scan_para = []
            
            n_input_pat = len(seq_entry_arr)
            
            scanner = MotifScanner()
            for i in range(n_input_pat):
                seq = seq_entry_arr[i].get()
                num = int(num_entry_arr[i].get())
                revcom_flag = revcom_arr[i].get()
                scanner.add_motif(seq, num, revcom_flag)
            
            infile = infile_entry.get()
            outfile = outfile_entry.get()
            
            self.motif_scan_para.append(scanner)
            self.motif_scan_para.append(infile)
            self.motif_scan_para.append(outfile)
            
            self.run_scan()
            
        n_max_pat = 10
        revcom_arr = [BooleanVar() for _ in range(n_max_pat)]
        seq_entry_arr = []
        num_entry_arr = []
        func_ind_gen = NumGenerator()
        callback_fun_arr = [Command(add_pattern, master, f'Motif {i+1}', revcom_arr[i], irow+i+1, seq_entry_arr, num_entry_arr)
                            for i in range(n_max_pat) ]
        # for cmd in callback_fun_arr:
        #     print(cmd.args)
        add_entry_button = Button(master, text="Add Motif", command=callback_fun_arr[func_ind_gen.get()])
        add_entry_button.grid(row=irow, column=4, sticky='w')

        self.progress_row = row_gen.get()
        irow = row_gen.get()
        run_button = Button(master, text="RUN", command=run_analysis)
        run_button.grid(row=irow, column=2, pady=50, sticky='w')
        
        quit_button = Button(master, text="QUIT", command=self.master.destroy)
        quit_button.grid(row=irow, column=4, pady=50, sticky='w')

        self.progress = Progressbar(master, orient=HORIZONTAL, length=300,  mode='indeterminate')

    def init_chipseq_gui_tabmain(self,master):

        # def enter_filename():
        #     file = filedialog.askopenfilename(initialdir='.',title = "Select file", filetypes = (("fasta files",".fasta .fa .gz"),("all files","*.*")))
        #     infile_entry.delete(0,END)
        #     infile_entry.insert(0,file)

        # def enter_outdir():
        #     outdir = filedialog.askdirectory(initialdir='.')
        #     outdir_entry.delete(0,END)
        #     outdir_entry.insert(0,outdir)

        def run_analysis():
            self.chip_para_dict = {"identifier":identifier_entry.get(),
                                   "file_name":infile_entry.get(),
                                   "out_dir":outdir_entry.get(),
                                   "min_kmer_len":int(min_kmer_len_entry.get()),
                                   "max_kmer_len":int(max_kmer_len_entry.get())}
            self.run_chipseq()

        # def validate_identifier(in_str):
        #     if in_str:
        #         return True
        #     else:
        #         return False
        # vcmd1 = (master.register(validate_identifier), '%P')

        # def identifiercolour(event):
        #     identifier_entry.config({"background": "light green"})
        #     #identifier_entry.config({"background": "tomoato"})

        # def infilecolour(event):
        #     if os.path.exists(str(infile_entry.get())):
        #         infile_entry.config({"background": "light green"})
        #     else:
        #         infile_entry.config({"background": "tomato"})

        # def outdircolour(event):
        #     if os.access(str(outdir_entry.get()), os.W_OK):
        #         outdir_entry.config({"background": "light green"})
        #     else:
        #         outdir_entry.config({"background": "tomato"})

        # def minkmercolour(event):
        #     if int(min_kmer_len_entry.get()) >= 0:
        #         min_kmer_len_entry.config({"background": "light green"})
        #     else:
        #         min_kmer_len_entry.config({"background": "tomato"})

        # def maxkmercolour(event):
        #     if int(max_kmer_len_entry.get()) >= int(min_kmer_len_entry.get()):
        #         max_kmer_len_entry.config({"background": "light green"})
        #     else:
        #         max_kmer_len_entry.config({"background": "tomato"})
        
        row_gen = NumGenerator()
        
        irow = row_gen.get()
        identifier_label, identifier_entry = gen_label_entry(master, "Identifier", "TF_name", grid_on=True, irow=irow)
        # identifier_label = Label(master, text="Identifier")
        # identifier_entry = Entry(master, validate="focusout", validatecommand=vcmd1)
        # identifier_entry.insert(END, 'TF_name')
        # identifier_entry['validate']="focusout"
        # identifier_entry['validatecommand']=vcmd1
        # identifier_label.grid(row=irow, column=0)
        # identifier_entry.grid(row=irow, column=1)
        # identifier_entry.bind("<FocusOut>", identifiercolour)

        # def validate_inputfile(in_str):
        #     return os.path.exists(in_str)
        # vcmd2 = (master.register(validate_inputfile), '%P')

        irow = row_gen.get()
        infile_label, infile_entry, infile_button = gen_file_entry(master, "Input file", 'path_to_input_fasta_file', 'Open File', 
                                                                   grid_on=True, irow=irow )
        
        # infile_label = Label(master, text="Input file")
        # infile_entry = Entry(master)
        # # infile_entry = Entry(master, validate="focusout", validatecommand=vcmd2)
        # infile_entry.insert(END, 'path_to_input_fasta_file')
        # infile_button = Button(master, text="Open File", command=enter_filename)
        # infile_label.grid(row=irow, column=0)
        # infile_entry.grid(row=irow, column=1)
        # infile_button.grid(row=irow, column=2)
        # infile_entry.bind("<FocusOut>", infilecolour)

        # irow = row_gen.get()
        # outdir_label = Label(master, text="Output Directory")
        # outdir_entry = Entry(master)
        # # outdir_entry.insert(END, 'path_to_output_directory')
        # outdir_button = Button(master, text="Open Directory", command=enter_outdir)
        
        irow = row_gen.get()
        outdir_label, outdir_entry, outdir_button = gen_directory_entry(master, "Output Directory", "", "Open Directory",
                                                                        grid_on=True, irow=irow)
        # outdir_label.grid(row=irow, column=0)
        # outdir_entry.grid(row=irow, column=1)
        # outdir_button.grid(row=irow, column=2)
        # outdir_entry.bind("<FocusOut>", outdircolour)

        # def validate_kmer_len(in_str):
        #     val = int(in_str)
        #     if val<1:
        #         return False
        #     else:
        #         return True
        # vcmd4 = (master.register(validate_kmer_len), '%P')

        irow = row_gen.get()
        min_kmer_len_label, min_kmer_len_entry = gen_label_entry(master, "Minimum Kmer Length", "", grid_on=True, irow=irow)
        # min_kmer_len_label = Label(master, text="Minimum Kmer Length")
        # min_kmer_len_entry = Entry(master)
        # min_kmer_len_entry = Entry(master, validate="key", validatecommand=vcmd4)
        # min_kmer_len_label.grid(row=irow, column=0)
        # min_kmer_len_entry.grid(row=irow, column=1)
        # min_kmer_len_entry.bind("<FocusOut>", minkmercolour)

        irow = row_gen.get()
        max_kmer_len_label, max_kmer_len_entry = gen_label_entry(master, "Maximum Kmer Length", "", grid_on=True, irow=irow)
        # max_kmer_len_label = Label(master, text="Maximum Kmer Length")
        # max_kmer_len_entry = Entry(master)
        # max_kmer_len_entry = Entry(master,validate="key", validatecommand=vcmd4)
        # max_kmer_len_label.grid(row=irow, column=0)
        # max_kmer_len_entry.grid(row=irow, column=1)
        # max_kmer_len_entry.bind("<FocusOut>", maxkmercolour)

        self.progress_row = row_gen.get()

        irow = row_gen.get()
        self.run_button = Button(master, text="Run", command=run_analysis)
        self.run_button.grid(row=irow,column=1, pady=10)

        quit_button = Button(master, text="QUIT", command=self.master.destroy)
        quit_button.grid(row=irow, column=2, pady=10)

        self.progress = Progressbar(master, orient=HORIZONTAL, length=300,  mode='indeterminate')
        
# see
# https://stackoverflow.com/questions/33768577/tkinter-gui-with-progress-bar
# https://www.youtube.com/watch?v=o_Ct13fHeck
    def run_chipseq(self):
        def chipseq(**kwargs):
            self.progress.grid(row=self.progress_row,column=0,columnspan=3)
            self.progress.start()
            csp = ChipSeqProcessor(**kwargs)
            csp.run()
            self.progress.stop()
            self.progress.grid_forget()

            messagebox.showinfo('Info', "Process completed!")
            self.run_button['state']='normal'

        # print(self.chip_para_dict)

        self.run_button['state']='disabled'

        threading.Thread(target=chipseq, kwargs=self.chip_para_dict).start()

    def init_selexseq_gui(self):
        self.clear()

        row_gen = NumGenerator()
        
        logo = Label(self, image=self.logoimage)
        logo.grid_columnconfigure(0, weight=1)
        logo.grid(row=row_gen.get())

        note = Notebook(self)
        main = Frame(note)
        motif_scan = Frame(note)
        topkmer = Frame(note)
        kmer_query = Frame(note)
        
        note.add(main, text = "SELEX-SEQ Main")
        note.add(motif_scan, text = "Motif Scan")
        note.add(topkmer, text='Top kmer')
        note.add(kmer_query, text='Kmer query')

        self.init_selexseq_gui_tabmain_step0(main)
        self.init_motifscan_guitab(motif_scan)
        self.init_topkmer_guitab(topkmer)
        self.init_kmerquery_guitab(kmer_query)

        note.grid(row=row_gen.get(), column=0)
        # self.columnconfigure(0,weight=1)
        
        self.master.grid_columnconfigure(0, weight=1)
        self.master.grid_rowconfigure(0, weight=1)
        self.master.grid_rowconfigure(1, weight=5)
        
    def init_selexseq_gui_tabmain_step0(self, master):
        # def validate_identifier(in_str):
        #     if in_str:
        #         return True
        #     else:
        #         return False
        # vcmd1 = (master.register(validate_identifier), '%P')

        def call_next():
            self.selex_gui_dict = {"identifier_label":identifier_label, "identifier_entry":identifier_entry,
             "min_round_label":min_round_label, "min_round_entry":min_round_entry,
             "max_round_label":max_round_label, "max_round_entry":max_round_entry, "next_button":next_button}
            self.selex_para_dict = {
                "curr_row":irow,
                "identifier":identifier_entry.get(),
                "min_selex_round":int(min_round_entry.get()),
                "max_selex_round":int(max_round_entry.get())
            }
            self.init_selexseq_gui_tabmain_step1(master)

        # def identifiercolour(event):
        #     identifier_entry.config({"background": "light green"})
        #     #identifier_entry.config({"background": "tomoato"})
        # def minroundcolour(event):
        #     if int(min_round_entry.get()) >= 0:
        #         min_round_entry.config({"background": "light green"})
        #     else:
        #         min_round_entry.config({"background": "tomato"})

        # def maxroundcolour(event):
        #     if int(max_round_entry.get()) >= int(min_round_entry.get()):
        #         max_round_entry.config({"background": "light green"})
        #     else:
        #         max_round_entry.config({"background": "tomato"})

        
        row_gen = NumGenerator()
        irow = row_gen.get()
        identifier_label, identifier_entry = gen_label_entry(master, "Identifier", "TF_name", grid_on=True, irow=irow)
        
        # identifier_label = Label(master, text="Identifier")
        # identifier_entry = Entry(master, validate="focusout", validatecommand=vcmd1)
        # identifier_entry.insert(END, 'TF_name')
        # identifier_label.grid(row=0, column=0)
        # identifier_entry.grid(row=0, column=1)
        # identifier_entry.bind("<FocusOut>", identifiercolour)

        irow = row_gen.get()
        min_round_label, min_round_entry = gen_label_entry(master, "Minimum SELEX round number", "", grid_on=True, irow=irow)
        
        # irow = 1
        # min_round_label = Label(master, text="Minimum SELEX round number")
        # min_round_entry = Entry(master)
        # # min_round_entry = Entry(master, validate="key", validatecommand=vcmd4)
        # min_round_label.grid(row=irow, column=0)
        # min_round_entry.grid(row=irow, column=1)
        # min_round_entry.bind("<FocusOut>", minroundcolour)

        irow = row_gen.get()
        max_round_label, max_round_entry = gen_label_entry(master, "Maximum SELEX round numer", "", grid_on=True, irow=irow)
        
        # irow = 2
        # max_round_label = Label(master, text="Maximum SELEX round numer")
        # max_round_entry = Entry(master)
        # # max_round_entry = Entry(master, validate="key", validatecommand=vcmd4)
        # max_round_label.grid(row=irow, column=0)
        # max_round_entry.grid(row=irow, column=1)
        # max_round_entry.bind("<FocusOut>", maxroundcolour)

        # irow = 3
        irow = row_gen.get()
        next_button = Button(master, text="Next", command=call_next)
        next_button.grid(row=irow,column=1, pady=10)

    def init_selexseq_gui_tabmain_step1(self, master):
        def call_prev():
            self.selex_para_dict['curr_row'] = curr_row
            for obj in selex_infile_entry_arr + selex_infile_label_arr + selex_infile_button_arr:
                obj.destroy()
            back_button.destroy()
            for key in self.selex_gui_dict:
                self.selex_gui_dict[key].destroy()
            self.init_selexseq_gui_tabmain_step0(master)

        def call_next():
            back_button.destroy()
            self.selex_gui_dict['next_button'].destroy()
            del self.selex_gui_dict['next_button']

            self.selex_para_dict['curr_row'] = curr_row
            self.selex_gui_dict['selex_infile_entry_arr'] = selex_infile_entry_arr

            file_name_arr = []
            for obj in selex_infile_entry_arr:
                file_name_arr.append(obj.get())
            self.selex_para_dict['file_name_arr'] = file_name_arr

            self.init_selexseq_gui_tabmain_step2(master)

        # def enter_filename(i):
        #     file = filedialog.askopenfilename(initialdir='.',title = "Select file", filetypes = (("fasta files",".fasta .fa .gz"),("all files","*.*")))
        #     selex_infile_entry_arr[i].delete(0,END)
        #     selex_infile_entry_arr[i].insert(0,file)

        self.selex_gui_dict['next_button'].grid_forget()

        # def roundfilescolour(event, relround):
        #     if os.path.exists(str(selex_infile_entry_arr[relround].get())):
        #         selex_infile_entry_arr[relround].config({"background": "light green"})
        #     else:
        #         selex_infile_entry_arr[relround].config({"background": "tomato"})

        curr_row = self.selex_para_dict['curr_row']
        selex_infile_label_arr = []
        selex_infile_entry_arr = []
        selex_infile_button_arr = []
        for i_row, i_round in enumerate(range(self.selex_para_dict["min_selex_round"],self.selex_para_dict["max_selex_round"]+1)):
            tmplabel, tmpentry, tmpbutton = gen_file_entry(master, f'Selex Round {i_round}', "", "Open File", grid_on=True, irow=curr_row+i_row)
            # tmplabel = Label(master, text=f'Selex Round {i_round}')
            # tmpentry = Entry(master)
            # #tmpentry.bind("<FocusOut>", lambda event : roundfilescolour(event, i_row))
            # tmpbutton = Button(master, text="Open File", command=Command(enter_filename, i_row) )
            # tmplabel.grid(row=curr_row+i_row, column=0)
            # tmpentry.grid(row=curr_row+i_row, column=1)
            # tmpbutton.grid(row=curr_row+i_row, column=2)
            selex_infile_label_arr.append(tmplabel)
            selex_infile_entry_arr.append(tmpentry)
            # selex_infile_entry_arr[i_row].bind("<FocusOut>", lambda event, relround=(i_row): roundfilescolour(event, relround))
            selex_infile_button_arr.append(tmpbutton)

        curr_row += i_row

        back_button =  Button(master, text="Back", command=call_prev)
        back_button.grid(row=curr_row+1,column=0, pady=10)

        self.selex_gui_dict['next_button']['command']=call_next
        self.selex_gui_dict['next_button'].grid(row=curr_row+1, column=2, pady=10)

    def init_selexseq_gui_tabmain_step2(self, master):
        def run_analysis():
            self.selex_para_dict["out_dir"]=outdir_entry.get()
            self.selex_para_dict["min_kmer_len"]=int(min_kmer_len_entry.get())
            self.selex_para_dict["max_kmer_len"]=int(max_kmer_len_entry.get())
            self.run_selexseq()

        # def enter_outdir():
        #     outdir = filedialog.askdirectory(initialdir='.')
        #     outdir_entry.delete(0,END)
        #     outdir_entry.insert(0,outdir)

        curr_row = self.selex_para_dict['curr_row']+1

        # def outdircolour(event):
        #     if os.access(str(outdir_entry.get()), os.W_OK):
        #         outdir_entry.config({"background": "light green"})
        #     else:
        #         outdir_entry.config({"background": "tomato"})

        # def minkmercolour(event):
        #     if int(min_kmer_len_entry.get()) >= 0:
        #         min_kmer_len_entry.config({"background": "light green"})
        #     else:
        #         min_kmer_len_entry.config({"background": "tomato"})

        # def maxkmercolour(event):
        #     if int(max_kmer_len_entry.get()) >= int(min_kmer_len_entry.get()):
        #         max_kmer_len_entry.config({"background": "light green"})
        #     else:
        #         max_kmer_len_entry.config({"background": "tomato"})

        # outdir_label = Label(master, text="Output Directory")
        # outdir_entry = Entry(master)
        # outdir_entry.insert(END, 'path_to_output_directory')
        # outdir_button = Button(master, text="Open Directory", command=enter_outdir)
        # outdir_label.grid(row=curr_row, column=0)
        # outdir_entry.grid(row=curr_row, column=1)
        # outdir_button.grid(row=curr_row, column=2)
        # outdir_entry.bind("<FocusOut>", outdircolour)
        
        row_gen = NumGenerator()
        irow = row_gen.get()
        outdir_label, outdir_entry, outdir_button = gen_directory_entry(master, "Output Directory", 'path_to_output_directory',"Open Directory", 
                                                                        grid_on=True, irow=curr_row+irow)


        irow = row_gen.get()
        min_kmer_len_label, min_kmer_len_entry = gen_label_entry(master, "Minimum Kmer Length", "", grid_on=True, irow=curr_row+irow)
        # min_kmer_len_label = Label(master, text="Minimum Kmer Length")
        # min_kmer_len_entry = Entry(master)
        # min_kmer_len_label.grid(row=curr_row+irow, column=0)
        # min_kmer_len_entry.grid(row=curr_row+irow, column=1)
        # min_kmer_len_entry.bind("<FocusOut>", minkmercolour)

        irow = row_gen.get()
        max_kmer_len_label, max_kmer_len_entry = gen_label_entry(master, "Maximum Kmer Length", "", grid_on=True, irow=curr_row+irow)
        
        # irow = 2
        # max_kmer_len_label = Label(master, text="Maximum Kmer Length")
        # max_kmer_len_entry = Entry(master)
        # max_kmer_len_label.grid(row=curr_row+irow, column=0)
        # max_kmer_len_entry.grid(row=curr_row+irow, column=1)
        # max_kmer_len_entry.bind("<FocusOut>", maxkmercolour)

        self.progress_row = row_gen.get() + curr_row # reserve for progress bar
        
        # irow = 3
        irow = row_gen.get()
        self.run_button = Button(master, text="Run", command=run_analysis)
        self.run_button.grid(row=curr_row+irow,column=1, pady=10)

        # irow = 3
        quit_button = Button(master, text="QUIT", command=self.master.destroy)
        quit_button.grid(row=curr_row+irow, column=2, pady=10)

        self.progress = Progressbar(master, orient=HORIZONTAL, length=300,  mode='indeterminate')

    def run_selexseq(self):
        def selexseq(**kwargs):
            self.progress.grid(row=self.progress_row,column=0,columnspan=3)
            self.progress.start()
            if 'curr_row' in kwargs:
                del kwargs['curr_row']
            ssp = SelexSeqProcessor(**kwargs)
            ssp.run()
            self.progress.stop()
            self.progress.grid_forget()

            messagebox.showinfo('Info', "Process completed!")
            self.run_button['state']='normal'

        self.run_button['state']='disabled'

        threading.Thread(target=selexseq, kwargs=self.selex_para_dict).start()

if __name__=="__main__":
    root = Tk()
    root.title("IniMotif")
    #root.geometry("800x600")

    # https://stackoverflow.com/questions/51973653/tkinter-grid-fill-empty-space
    screen_width = root.winfo_screenwidth()
    screen_height = root.winfo_screenheight()

    app = Application(master=root)
    app.mainloop()
