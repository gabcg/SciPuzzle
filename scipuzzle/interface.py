#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import tkinter

main_window = tkinter.Tk()

input = tkinter.filedialog.askdirectory()

output = tkinter.filedialog.asksaveasfilename(initialdir="/",
                                              title="Select file",
                                              filetypes=(("PDB files", "*.pdb")))

input.pack()
output.pack()

main_window.mainloop()
