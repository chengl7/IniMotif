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

import sys
#import subprocess as sub
from inimotif_main import FileProcessor,ChipSeqProcessor

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
        
#        self.rowconfigure(0, weight=1)
#        self.columnconfigure(0, weight=1)
        
        self.run_func = None
        self.run_button = None
        self.progress = None
        
#        self.pack()
        
#        self.width = width
#        self.height = height
        
        self.grid(row=0, column=0,sticky='nswe')
        self.grid_propagate(False)
        self.select_analysis()
#        self.create_widgets()
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
        
    
#    def run_chip_analysis(self):
##        p = sub.Popen([sys.executable, os.path.dirname(os.path.realpath(__file__))+os.path.sep+"inimotif.py",
##                     sys.argv[1]]]) 
##        p = sub.Popen('./script',stdout=sub.PIPE,stderr=sub.PIPE)
##        output, errors = p.communicate()
#        
#        self.proc = Process(target=FileProcessor, args=(), kwargs=self.chip_para_dict)
#        self.proc.start()
#        #self.proc.join()   # will block untill process finished
        
        
#https://stackoverflow.com/questions/665566/redirect-command-line-results-to-a-tkinter-gui
#from tkinter import *
#import subprocess as sub
#p = sub.Popen('./script',stdout=sub.PIPE,stderr=sub.PIPE)
#output, errors = p.communicate()
#
#root = Tk()
#text = Text(root)
#text.pack()
#text.insert(END, output)
#root.mainloop()
    
#    def start_chip_seq(self):
#        self.state = "chip-seq"
#        print("start chip-seq gui")
#        self.clear()
#        self.create_widgets()
#        
#        self.btn.grid(row=0,column=0)
#        self.progress = Progressbar(self, orient=HORIZONTAL,length=100,  mode='indeterminate')
#        
#        pass
    
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
            file = filedialog.askopenfilename()
            infile_entry.delete(0,END)
            infile_entry.insert(0,file)
        
        self.chip_para_dict = {"file_name":"/Users/lcheng/Documents/github/IniMotif-py/exampledata/NF1-1",
                               "out_dir":"./test",
                               "identifier":"NF1",
                               "min_kmer_len":10,
                               "max_kmer_len":12}
        
        identifier_label = Label(master, text="Identifier")
        identifier_entry = Entry(master)
        identifier_label.grid(row=0, column=0)
        identifier_entry.grid(row=0, column=1)
        
        irow = 1
        infile_label = Label(master, text="Input file")
        infile_entry = Entry(master)
        infile_button = Button(master, text="Open File", command=enter_filename)
        infile_label.grid(row=irow, column=0)
        infile_entry.grid(row=irow, column=1)
        infile_button.grid(row=irow, column=2)
        
        irow = 2
        outdir_label = Label(master, text="Output Directory")
        outdir_entry = Entry(master)
        outdir_label.grid(row=irow, column=0)
        outdir_entry.grid(row=irow, column=1)
        
        irow = 3
        kmer_len_label = Label(master, text="Kmer Length")
        kmer_len_entry = Entry(master)
        kmer_len_label.grid(row=irow, column=0)
        kmer_len_entry.grid(row=irow, column=1)
        
        
        
        # collect input parameters, TODO
        self.chip_para_dict = {"file_name":"/Users/lcheng/Documents/github/IniMotif-py/exampledata/NF1-1",
                               "out_dir":"./test",
                               "identifier":"NF1",
                               "min_kmer_len":10,
                               "max_kmer_len":12}
        
        
        irow = 5
        self.run_func = self.run_chipseq
        self.run_button = Button(master, text="Run",command=self.run_func)
        self.run_button.grid(row=irow,column=1, pady=10)
        
        irow = 5
        quit_button = Button(master, text="QUIT", command=self.master.destroy)
        quit_button.grid(row=irow, column=2, pady=10)
        
        
        self.progress = Progressbar(master, orient=HORIZONTAL, length=300,  mode='indeterminate')
# see    
# https://stackoverflow.com/questions/33768577/tkinter-gui-with-progress-bar
# https://www.youtube.com/watch?v=o_Ct13fHeck
    def run_chipseq(self):
        def chipseq(**kwargs):
            self.progress.grid(row=4,column=0,columnspan=3)
            self.progress.start()
            csp = ChipSeqProcessor(**kwargs)
            csp.run()
            self.progress.stop()
            self.progress.grid_forget()
            
            messagebox.showinfo('Info', "Process completed!")
            self.run_button['state']='normal'

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
    
#    def create_widgets(self):
##        self.hi_there = Button(self)
##        self.hi_there["text"] = "Hello World\n(click me)"
##        self.hi_there["command"] = self.say_hi
##        self.hi_there.pack(side="top")
##
##        self.quit = Button(self, text="QUIT",
##                              command=self.master.destroy)
##        self.quit.pack(side="bottom")
#        
#        note = Notebook()
#        self.notebook  = note
#        tab1 = Frame(note)
#        tab2 = Frame(note)
#        tab3 = Frame(note)
#        note.add(tab1, text = "Tracing", compound=TOP)
#        note.add(tab2, text = "Network Details")
#        note.add(tab3, text = "Tab Three")
#        note.grid()
#        
#        self.hi_there = Button(tab1, text="identifier",command=self.say_hi)
#        self.hi_there.grid(row=0)
#        
#        text = Text(tab1)
#        text.grid(row=1)
#        text.insert(END, sys.stdout)
#        
#        self.run_chip = Button(tab1, text="run chip",command=self.run_chip_analysis)
#        self.run_chip.grid(row=1)
#        
#        self.quit = Button(tab2, text="QUIT", command=self.master.destroy)
#        self.quit.grid(row=0)
#
#    def say_hi(self):
#        print("hi there, everyone!")

#style = Style()
#style.configure("LargeButton", borderwidth=5, relief="sunken", width=600, height=200)

root = Tk()
root.title("IniMotif")
#root.geometry("800x600")

# https://stackoverflow.com/questions/51973653/tkinter-grid-fill-empty-space
screen_width = root.winfo_screenwidth()
screen_height = root.winfo_screenheight()

#frame = Frame(root, borderwidth=5, relief="sunken", width=0.5*screen_width, height=0.5*screen_height)

#frame1.rowconfigure(3, minsize=30)
#frame1.columnconfigure(3, minsize=30)
#content.grid(column=0, row=0)
#frame.grid(row=0,column=0)
#root.maxsize(1000, 400)
#app = Application(master=frame)

app = Application(master=root)
app.mainloop()

#root.mainloop()
