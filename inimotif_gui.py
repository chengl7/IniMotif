#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  2 22:35:43 2019

@author: lcheng
"""
from tkinter import filedialog
from tkinter import messagebox
from tkinter import Tk, END, HORIZONTAL, PhotoImage, N, W, E, S
from tkinter import StringVar, IntVar, BooleanVar
from tkinter.ttk import Button, Label, Entry, Frame, Notebook, Progressbar, Style, Combobox, Radiobutton, Checkbutton
import os
from inimotif_main import FileProcessor,ChipSeqProcessor,SelexSeqProcessor
from inimotif_util import Masker
#from tkinter import *

import time
import threading

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
    app_style.configure('grey.TEntry', fg="grey")

def gen_label_entry(master, label_text, entry_text):
    label = Label(master, text=label_text)
    entry = Entry(master, style='grey.TEntry')
    entry.insert(END, entry_text)
    return label, entry

def gen_file_entry(master, label_text, entry_text, button_text):
    def enter_filename():
        file = filedialog.askopenfilename(initialdir='.',title = "Select file", filetypes = (("fasta files",".fasta .fa .gz"),("all files","*.*")))
        entry.delete(0,END)
        entry.insert(0,file)
    label = Label(master, text=label_text)
    entry = Entry(master, style='grey.TEntry')
    entry.insert(END, entry_text)
    
    button = Button(master, text=button_text, command=enter_filename)
    return label, entry, button

def gen_directory_entry(master, label_text, entry_text, button_text):
    def enter_outdir():
        outdir = filedialog.askdirectory(initialdir='.')
        entry.delete(0,END)
        entry.insert(0,outdir)
    label = Label(master, text=label_text)
    entry = Entry(master, style='grey.TEntry')
    entry.insert(END, entry_text)
    
    button = Button(master, text=button_text, command=enter_outdir)
    return label, entry, button

def gen_mask_pattern(master, label_text, type_value, revcom_value):
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
    return label, combobox, check, entry1, entry2

class Application(Frame):
    def __init__(self, master=None):
        super().__init__(master, borderwidth=5, width=0.8*screen_width, height=0.7*screen_height)
        self.master = master
        self.proc = None
        
        self.main_widgets = []
        self.logoimage = PhotoImage(file="GUIgraphics/logo.png")
        self.maskerimage = PhotoImage(file="GUIgraphics/masker.png")

        self.run_func = None
        self.run_button = None
        self.progress = None
        self.progress_row = None
        
        set_style()

        self.grid(row=0, column=0,sticky=N+S+W+E)
        self.grid_propagate(False)

        self.select_analysis()

        self.chip_para_dict = {}
        self.selex_para_dict = {}
        self.mask_para = []

        self.selex_gui_dict = {}
        
        self.master.columnconfigure(0, weight=1)
        self.master.rowconfigure(0, weight=1)

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
        self.columnconfigure(0, weight=1)    
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
        add_entry_button.grid(row=irow, column=5)
        add_entry_button.grid(row=irow, column=5)

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
        row_gen = NumGenerator()

        logo = Label(self, image=self.logoimage)
        logo.grid_columnconfigure(0, weight=1)
        logo.grid(row=row_gen.get())
        # self.logo.place(relx=0.5, rely=0.1, anchor="n")

        note = Notebook(self)
        self.notebook  = note
        main = Frame(note)
        win_extract = Frame(note)
        note.add(main, text = "Chip-Seq Main")
        note.add(win_extract, text = "Chip-Seq Window Extract")

        self.init_chipseq_gui_tabmain(main)

        # note.grid_rowconfigure(0, weight=1)
        note.grid_columnconfigure(0, weight=1)
        # note.grid(row=0,column=0, sticky=N, pady=(250, 10))
        note.grid(row=row_gen.get(),column=0)


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
            self.run_chipseq()

        def validate_identifier(in_str):
            if in_str:
                return True
            else:
                return False
        vcmd1 = (master.register(validate_identifier), '%P')

        def identifiercolour(event):
            identifier_entry.config({"background": "light green"})
            #identifier_entry.config({"background": "tomoato"})

        def infilecolour(event):
            if os.path.exists(str(infile_entry.get())):
                infile_entry.config({"background": "light green"})
            else:
                infile_entry.config({"background": "tomato"})

        def outdircolour(event):
            if os.access(str(outdir_entry.get()), os.W_OK):
                outdir_entry.config({"background": "light green"})
            else:
                outdir_entry.config({"background": "tomato"})

        def minkmercolour(event):
            if int(min_kmer_len_entry.get()) >= 0:
                min_kmer_len_entry.config({"background": "light green"})
            else:
                min_kmer_len_entry.config({"background": "tomato"})

        def maxkmercolour(event):
            if int(max_kmer_len_entry.get()) >= int(min_kmer_len_entry.get()):
                max_kmer_len_entry.config({"background": "light green"})
            else:
                max_kmer_len_entry.config({"background": "tomato"})
        
        identifier_label, identifier_entry = gen_label_entry(master, "Identifier", "TF_name")
        # identifier_label = Label(master, text="Identifier")
        # identifier_entry = Entry(master, validate="focusout", validatecommand=vcmd1)
        # identifier_entry.insert(END, 'TF_name')
        identifier_entry['validate']="focusout"
        identifier_entry['validatecommand']=vcmd1
        identifier_label.grid(row=0, column=0)
        identifier_entry.grid(row=0, column=1)
        identifier_entry.bind("<FocusOut>", identifiercolour)

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
        infile_entry.bind("<FocusOut>", infilecolour)

        irow = 2
        outdir_label = Label(master, text="Output Directory")
        outdir_entry = Entry(master)
        # outdir_entry.insert(END, 'path_to_output_directory')
        outdir_button = Button(master, text="Open Directory", command=enter_outdir)
        outdir_label.grid(row=irow, column=0)
        outdir_entry.grid(row=irow, column=1)
        outdir_button.grid(row=irow, column=2)
        outdir_entry.bind("<FocusOut>", outdircolour)

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
        min_kmer_len_entry.bind("<FocusOut>", minkmercolour)

        irow = 4
        max_kmer_len_label = Label(master, text="Maximum Kmer Length")
        max_kmer_len_entry = Entry(master)
        # max_kmer_len_entry = Entry(master,validate="key", validatecommand=vcmd4)
        max_kmer_len_label.grid(row=irow, column=0)
        max_kmer_len_entry.grid(row=irow, column=1)
        max_kmer_len_entry.bind("<FocusOut>", maxkmercolour)


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

    def init_selexseq_gui(self):
        self.clear()

        row_gen = NumGenerator()
        
        logo = Label(self, image=self.logoimage)
        logo.grid_columnconfigure(0, weight=1)
        logo.grid(row=row_gen.get())

        note = Notebook(self)
        # self.notebook  = note
        main = Frame(note)
        note.add(main, text = "SELEX-SEQ Main")

        self.init_selexseq_gui_tabmain_step0(main)

        # note.grid_rowconfigure(0, weight=1)
        note.grid_columnconfigure(0, weight=1)
        # note.grid(row=0,column=0, sticky=N, pady=(250, 10))
        note.grid(row=row_gen.get(), column=0)

    def init_selexseq_gui_tabmain_step0(self, master):
        def validate_identifier(in_str):
            if in_str:
                return True
            else:
                return False
        vcmd1 = (master.register(validate_identifier), '%P')

        def call_next():
            self.selex_gui_dict = {"identifier_label":identifier_label, "identifier_entry":identifier_entry,
             "min_round_label":min_round_label, "min_round_entry":min_round_entry,
             "max_round_label":max_round_label, "max_round_entry":max_round_entry, "next_button":next_button}
            self.selex_para_dict = {
                "curr_row":3,
                "identifier":identifier_entry.get(),
                "min_selex_round":int(min_round_entry.get()),
                "max_selex_round":int(max_round_entry.get())
            }
            self.init_selexseq_gui_tabmain_step1(master)

        def identifiercolour(event):
            identifier_entry.config({"background": "light green"})
            #identifier_entry.config({"background": "tomoato"})
        def minroundcolour(event):
            if int(min_round_entry.get()) >= 0:
                min_round_entry.config({"background": "light green"})
            else:
                min_round_entry.config({"background": "tomato"})

        def maxroundcolour(event):
            if int(max_round_entry.get()) >= int(min_round_entry.get()):
                max_round_entry.config({"background": "light green"})
            else:
                max_round_entry.config({"background": "tomato"})


        identifier_label = Label(master, text="Identifier")
        identifier_entry = Entry(master, validate="focusout", validatecommand=vcmd1)
        identifier_entry.insert(END, 'TF_name')
        identifier_label.grid(row=0, column=0)
        identifier_entry.grid(row=0, column=1)
        identifier_entry.bind("<FocusOut>", identifiercolour)

        irow = 1
        min_round_label = Label(master, text="Minimum SELEX round number")
        min_round_entry = Entry(master)
        # min_round_entry = Entry(master, validate="key", validatecommand=vcmd4)
        min_round_label.grid(row=irow, column=0)
        min_round_entry.grid(row=irow, column=1)
        min_round_entry.bind("<FocusOut>", minroundcolour)

        irow = 2
        max_round_label = Label(master, text="Maximum SELEX round numer")
        max_round_entry = Entry(master)
        # max_round_entry = Entry(master, validate="key", validatecommand=vcmd4)
        max_round_label.grid(row=irow, column=0)
        max_round_entry.grid(row=irow, column=1)
        max_round_entry.bind("<FocusOut>", maxroundcolour)

        irow = 3
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

        def enter_filename(i):
            file = filedialog.askopenfilename(initialdir='.',title = "Select file", filetypes = (("fasta files",".fasta .fa .gz"),("all files","*.*")))
            selex_infile_entry_arr[i].delete(0,END)
            selex_infile_entry_arr[i].insert(0,file)

        self.selex_gui_dict['next_button'].grid_forget()

        def roundfilescolour(event, relround):
            if os.path.exists(str(selex_infile_entry_arr[relround].get())):
                selex_infile_entry_arr[relround].config({"background": "light green"})
            else:
                selex_infile_entry_arr[relround].config({"background": "tomato"})

        curr_row = self.selex_para_dict['curr_row']
        selex_infile_label_arr = []
        selex_infile_entry_arr = []
        selex_infile_button_arr = []
        for i_row, i_round in enumerate(range(self.selex_para_dict["min_selex_round"],self.selex_para_dict["max_selex_round"]+1)):
            tmplabel = Label(master, text=f'Selex Round {i_round}')
            tmpentry = Entry(master)
            #tmpentry.bind("<FocusOut>", lambda event : roundfilescolour(event, i_row))
            tmpbutton = Button(master, text="Open File", command=Command(enter_filename, i_row) )
            tmplabel.grid(row=curr_row+i_row, column=0)
            tmpentry.grid(row=curr_row+i_row, column=1)
            tmpbutton.grid(row=curr_row+i_row, column=2)
            selex_infile_label_arr.append(tmplabel)
            selex_infile_entry_arr.append(tmpentry)
            selex_infile_entry_arr[i_row].bind("<FocusOut>", lambda event, relround=(i_row): roundfilescolour(event, relround))
            selex_infile_button_arr.append(tmpbutton)

        curr_row += i_row

        back_button =  Button(master, text="Back", command=call_prev)
        back_button.grid(row=curr_row+1,column=0, pady=10)

        self.selex_gui_dict['next_button']['command']=call_next
        self.selex_gui_dict['next_button'].grid(row=curr_row+1, column=2, pady=10)

    def init_selexseq_gui_tabmain_step2(self, master):
        def run_analysis():
            self.selex_para_dict["curr_row"] = curr_row
            self.selex_para_dict["out_dir"]=outdir_entry.get()
            self.selex_para_dict["min_kmer_len"]=int(min_kmer_len_entry.get())
            self.selex_para_dict["max_kmer_len"]=int(max_kmer_len_entry.get())
            self.run_selexseq()

        def enter_outdir():
            outdir = filedialog.askdirectory(initialdir='.')
            outdir_entry.delete(0,END)
            outdir_entry.insert(0,outdir)

        curr_row = self.selex_para_dict['curr_row']+1

        def outdircolour(event):
            if os.access(str(outdir_entry.get()), os.W_OK):
                outdir_entry.config({"background": "light green"})
            else:
                outdir_entry.config({"background": "tomato"})

        def minkmercolour(event):
            if int(min_kmer_len_entry.get()) >= 0:
                min_kmer_len_entry.config({"background": "light green"})
            else:
                min_kmer_len_entry.config({"background": "tomato"})

        def maxkmercolour(event):
            if int(max_kmer_len_entry.get()) >= int(min_kmer_len_entry.get()):
                max_kmer_len_entry.config({"background": "light green"})
            else:
                max_kmer_len_entry.config({"background": "tomato"})

        outdir_label = Label(master, text="Output Directory")
        outdir_entry = Entry(master)
        outdir_entry.insert(END, 'path_to_output_directory')
        outdir_button = Button(master, text="Open Directory", command=enter_outdir)
        outdir_label.grid(row=curr_row, column=0)
        outdir_entry.grid(row=curr_row, column=1)
        outdir_button.grid(row=curr_row, column=2)
        outdir_entry.bind("<FocusOut>", outdircolour)


        irow = 1
        min_kmer_len_label = Label(master, text="Minimum Kmer Length")
        min_kmer_len_entry = Entry(master)
        min_kmer_len_label.grid(row=curr_row+irow, column=0)
        min_kmer_len_entry.grid(row=curr_row+irow, column=1)
        min_kmer_len_entry.bind("<FocusOut>", minkmercolour)

        irow = 2
        max_kmer_len_label = Label(master, text="Maximum Kmer Length")
        max_kmer_len_entry = Entry(master)
        max_kmer_len_label.grid(row=curr_row+irow, column=0)
        max_kmer_len_entry.grid(row=curr_row+irow, column=1)
        max_kmer_len_entry.bind("<FocusOut>", maxkmercolour)

        irow = 3
        self.run_button = Button(master, text="Run", command=run_analysis)
        self.run_button.grid(row=curr_row+irow,column=1, pady=10)

        irow = 3
        quit_button = Button(master, text="QUIT", command=self.master.destroy)
        quit_button.grid(row=curr_row+irow, column=2, pady=10)

        curr_row += 4  # reserve for progress bar

        self.progress = Progressbar(master, orient=HORIZONTAL, length=300,  mode='indeterminate')

    def run_selexseq(self):
        def selexseq(**kwargs):
            self.progress.grid(row=self.selex_para_dict["curr_row"],column=0,columnspan=3)
            self.progress.start()
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
