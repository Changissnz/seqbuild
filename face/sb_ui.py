import tkinter as tk
from tkinter import filedialog,font
from .comm_lang import * 

"""
the main Tkinter application class for seqbuild user 
interface (SB UI). 
"""
class SBApplication(tk.Frame):

    def __init__(self, master=None):
        tk.Frame.__init__(self, master)

        self.clp = CommLangParser(None) 
        self.grid()
        self.createWidgets()

    def createWidgets(self):
        self.quitButton = tk.Button(self, text='Quit',
            command=self.destroy)
        self.quitButton.grid()

        # Create a Text widget to display the content
        bold_font = font.Font(family="Arial", size=12, weight="bold")
        self.text_widget = tk.Text(self, wrap="word", font=bold_font,bg="red",fg="white",width=80, height=15)
        self.text_widget.grid()

        # Create a button to open the file
        self.open_button = tk.Button(self, text="OpEn CMD FiLe", command=self.open_file__UI)
        self.open_button.grid(row=0,column=0)

        self.baseview_button = tk.Button(self,text="bASe ViEw",command=self.reset_to_baseview)
        self.baseview_button.grid(row=0,column=10)


        # Create a button to clear the SHOW window. 
        def tw_del(): 
            return self.text_widget.delete(1.0, tk.END)
        self.clear_button = tk.Button(self,text="CLeAr ScReeN",command=tw_del)
        self.clear_button.grid(row=14,column=0) 

        self.fullreset_button = tk.Button(self,text="fuLL REseT",command=self.fullreset_clp)
        self.fullreset_button.grid(row=14,column=10)

        # Create another Text widget to write info. 
        self.text_widget2 = tk.Text(self, wrap="word", width=80, height=15)
        self.text_widget2.grid()

        # Create a button to send written commands 
        self.send_button = tk.Button(self, text="SeND cMd", command=self.send_cmd__UI)
        self.send_button.grid()


    def open_file__UI(self): 
        file_path = filedialog.askopenfilename(
            title="Select a Text File", \
            filetypes=[("Text files", "*.txt")])
        if file_path:
            self.clp.reload_file(file_path) 
            self.clp.process_file()
            self.reset_to_baseview() 
            return

    def reset_to_baseview(self): 
        self.text_widget.delete(1.0, tk.END)
        self.text_widget.insert(tk.END, str(self.clp))  

    def send_cmd__UI(self): 
        cmd = self.text_widget2.get("1.0", tk.END)
        self.text_widget2.delete("1.0", tk.END)
        cmd_lines = cmd.split("\n")
        cmd_lines = [cline.strip() for cline in cmd_lines]
        self.clp.load_cmd_lines(cmd_lines)
        self.clp.process_cmdlines()

        self.reset_to_baseview()
        self.process_SHOW_cmd() 
        return cmd 

    def process_SHOW_cmd(self): 
        if len(self.clp.show_commands) == 0: 
            return 
        
        self.text_widget.delete(1.0,tk.END) 
        while len(self.clp.show_commands) > 0:
            q = self.clp.show_commands.pop(0) 
            splitstr_cmd = q.split(" ")
            splitstr_cmd = [s.strip(".") for s in splitstr_cmd] 
            splitstr_cmd = [c_ for c_ in splitstr_cmd if len(c_) > 0]
            q = splitstr_cmd

            assert q[0] == "show" 
            assert len(q) == 2 

            assert q[1] in self.clp.vartable 

            x = self.clp.vartable[q[1]]
            self.text_widget.insert(tk.END, "SHOW {}\n{}\n".format(q[1],str(x)))

    def fullreset_clp(self): 
        self.clp = CommLangParser(None) 
        self.reset_to_baseview()

def run_sb_app(): 
    app = SBApplication()
    app.master.title('SeQbUiLd WiNdOw')  
    app.mainloop()