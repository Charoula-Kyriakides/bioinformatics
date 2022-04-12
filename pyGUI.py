import PySimpleGUI as sg
from genomics_class import *

# The names of algorithms that can be called 

funcs = {"Highest GC content":"gc_content", 
        "Nucleotide Frequency":"frequency", 
        "Common Ancestor":"common_ancestor"}

# A button for each algorithm
alg_buttons = [sg.Button(i, size=(19,3)) for i in funcs.keys()]

sg.theme('DarkBlue7')

layout = [
    [sg.T("")],
    [sg.Text("Chose a FASTA file: ")],
    [sg.Input(size=(60,5)), sg.FileBrowse(key="-FASTA-", file_types=(("Text Files","*.txt"),))],
    [sg.Button("Upload")],
    [sg.T("")],
    [sg.Text("Results will appear here.", key='-output-', size=(45, 10), font='25', justification='centre')],
    [sg.T("")],
    alg_buttons
    ]


def pick_event(event):
    """Picks an algorithm depending on the button clicked"""
    
    out = "The aglorithm was unsecessful. \n\nPlease check your file." # Print out this message in case algorithm fails
    for key in funcs.keys():
        if key == event:
            fn = "genome." + funcs[event] + "()"
            out = eval(fn)
            
    # Present the results
    window['-output-'].update(value=out)

    # Change the color of the buttons so that only one appears selected
    rems = list(funcs.keys())
    rems.remove(event)
    for a in rems:
        window[a].Update(button_color=('white','dark blue'))
    
    window[event].Update(button_color=('black', 'yellow'))


# create window
window = sg.Window("Genomics Calculations", layout)
genome = None
# create an event loop
while True:
    event, values = window.read()
    # End program if user closes window
    if event==sg.WIN_CLOSED:
        window.close()
        break

    elif event == "Upload":
        if values["-FASTA-"]:   # Checks if file is selected
            FASTA = values["-FASTA-"]
            genome = Read_Genome(FASTA)
            window['-output-'].update(value="File succesfully uploaded.\n\n Now you can select an algorithm.")

    elif (event in funcs.keys()) and genome:
    # Checks if any of the algorithm buttons are clicked and runs the algorithm
        pick_event(event)
