#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  2 22:35:43 2019

@author: lcheng
"""
from tkinter import filedialog
from tkinter import messagebox
from tkinter import *
from tkinter.ttk import *

import os
from inimotif_main import FileProcessor,ChipSeqProcessor,SelexSeqProcessor

import time
import threading

class Application(Frame):
    def __init__(self, master=None):
#        screen_width = master.winfo_screenwidth()
#        screen_height = master.winfo_screenheight()
#        print(screen_width, screen_height)
        super().__init__(master, borderwidth=5, width=0.5*screen_width, height=0.5*screen_height)
        self.master = master
        self.proc = None
        
        self.run_func = None
        self.run_button = None
        self.progress = None
        
        self.grid(row=0, column=0,sticky='nswe')
        self.grid_propagate(False)
        self.select_analysis()

        self.chip_para_dict = {}
        self.selex_para_dict = {}
        
    def select_analysis(self):
        self.l1 = Label(self, text="Choose Analysis Type")
        self.l1.grid(row=0,pady=10)
        
        self.chip = Button(self, text="Chip-seq", command=self.init_chipseq_gui)
        self.chip.grid(row=1,pady=10)
        
        self.selex = Button(self, text="Selex-seq", command=self.start_selex_seq)
        self.selex.grid(row=2,pady=10)
        
        self.l2 = Label(self, text="")
        self.l2.grid(row=3, padx=50, pady=10)
        
        self.quit = Button(self, text="QUIT", command=self.master.destroy)
        self.quit.grid(row=4, pady=10, ipadx=10, ipady=5)
    
    def init_chipseq_gui(self):
        self.clear()
        
        note = Notebook()
        self.notebook  = note
        main = Frame(note)
        masker = Frame(note)
        win_extract = Frame(note)
        note.add(main, text = "Chip-Seq Main")
        note.add(masker, text = "Masker")
        note.add(win_extract, text = "Chip-Seq Window Extract")
        
        self.init_chipseq_gui_tabmain(main)
        
        note.grid_rowconfigure(0, weight=1)
        note.grid_columnconfigure(0, weight=1)
        note.grid(row=0,column=0)
    
    def init_chipseq_gui_tabmain(self,master):
        
        def enter_filename():
            file = filedialog.askopenfilename(initialdir='.',title = "Select file", filetypes = (("fasta files",".fasta .fa .gz"),("all files","*.*")))
            infile_entry.delete(0,END)
            infile_entry.insert(0,file)
        
        def enter_outdir():
            outdir = filedialog.askdirectory(initialdir='.')
            outdir_entry.delete(0,END)
            outdir_entry.insert(0,outdir)
        
        def run_analysis():
            self.chip_para_dict = {"identifier":identifier_entry.get(),
                                   "file_name":infile_entry.get(),
                                   "out_dir":outdir_entry.get(),
                                   "min_kmer_len":int(min_kmer_len_entry.get()),
                                   "max_kmer_len":int(max_kmer_len_entry.get())}
            self.self.run_chipseq()

        def validate_identifier(in_str):
            if in_str:
                return True
            else:
                return False
        vcmd1 = (master.register(validate_identifier), '%P')

        identifier_label = Label(master, text="Identifier")
        identifier_entry = Entry(master, validate="focusout", validatecommand=vcmd1)
        identifier_entry.insert(END, 'TF_name')
        identifier_label.grid(row=0, column=0)
        identifier_entry.grid(row=0, column=1)
        
        # def validate_inputfile(in_str):
        #     return os.path.exists(in_str)
        # vcmd2 = (master.register(validate_inputfile), '%P')

        irow = 1
        infile_label = Label(master, text="Input file")
        infile_entry = Entry(master)
        # infile_entry = Entry(master, validate="focusout", validatecommand=vcmd2)
        infile_entry.insert(END, 'path_to_input_fasta_file')
        infile_button = Button(master, text="Open File", command=enter_filename)
        infile_label.grid(row=irow, column=0)
        infile_entry.grid(row=irow, column=1)
        infile_button.grid(row=irow, column=2)
        
        irow = 2
        outdir_label = Label(master, text="Output Directory")
        outdir_entry = Entry(master)
        outdir_entry.insert(END, 'path_to_output_directory')
        outdir_button = Button(master, text="Open Directory", command=enter_outdir)
        outdir_label.grid(row=irow, column=0)
        outdir_entry.grid(row=irow, column=1)
        outdir_button.grid(row=irow, column=2)
        
        # def validate_kmer_len(in_str):
        #     val = int(in_str)
        #     if val<1:
        #         return False
        #     else:
        #         return True
        # vcmd4 = (master.register(validate_kmer_len), '%P')

        irow = 3
        min_kmer_len_label = Label(master, text="Minimum Kmer Length")
        min_kmer_len_entry = Entry(master)
        # min_kmer_len_entry = Entry(master, validate="key", validatecommand=vcmd4)
        min_kmer_len_label.grid(row=irow, column=0)
        min_kmer_len_entry.grid(row=irow, column=1)
        
        irow = 4
        max_kmer_len_label = Label(master, text="Maximum Kmer Length")
        max_kmer_len_entry = Entry(master)
        # max_kmer_len_entry = Entry(master,validate="key", validatecommand=vcmd4)
        max_kmer_len_label.grid(row=irow, column=0)
        max_kmer_len_entry.grid(row=irow, column=1)
        
        irow = 6
        self.run_button = Button(master, text="Run", command=run_analysis)
        self.run_button.grid(row=irow,column=1, pady=10)
        
        irow = 6
        quit_button = Button(master, text="QUIT", command=self.master.destroy)
        quit_button.grid(row=irow, column=2, pady=10)
        
        self.progress = Progressbar(master, orient=HORIZONTAL, length=300,  mode='indeterminate')
# see    
# https://stackoverflow.com/questions/33768577/tkinter-gui-with-progress-bar
# https://www.youtube.com/watch?v=o_Ct13fHeck
    def run_chipseq(self):
        def chipseq(**kwargs):
            self.progress.grid(row=5,column=0,columnspan=3)
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
    
    def start_selex_seq(self):
        self.state = 'selex-seq'
        print("start selex-seq gui")
        pass
    
    # clear all widgets in the main panel
    def clear(self):
        self.l1.destroy()
        self.l2.destroy()
        self.chip.destroy()
        self.selex.destroy()
        self.quit.destroy()
#        self.chip.pack_forget()
#        self.selex.pack_forget()
#        self.quit.pack_forget()

if __name__=="__main__":    
    root = Tk()
    root.title("IniMotif")
    #root.geometry("800x600")

    # https://stackoverflow.com/questions/51973653/tkinter-grid-fill-empty-space
    screen_width = root.winfo_screenwidth()
    screen_height = root.winfo_screenheight()

    app = Application(master=root)
    app.mainloop()
