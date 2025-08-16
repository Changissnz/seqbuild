import tkinter as tk
from tkinter import filedialog,font
from .comm_lang import * 
from .cl_guide_parser import * 


"""
the main Tkinter application class for seqbuild user 
interface (SB UI). 
"""
class SBApplication(tk.Frame):

    def __init__(self, master=None):
        tk.Frame.__init__(self, master)

        self.clp = CommLangParser(None) 
        self.clgp = CLGuideParser() 
        self.clgp.process() 

        self.grid()
        self.create_widgets()

        # bool 
        self.is_helpwindow = False 

    def create_widgets(self):
        self.init_primary_window_detils() 
        self.set_primary_window_details()
        self.text_widget3 = None 

    def init_primary_window_detils(self):
        
        bold_font = font.Font(family="Arial", size=12, weight="bold")

        self.open_button = tk.Button(self, text="OpEn CMD FiLe", command=self.open_file__UI)
        self.text_widget = tk.Text(self, wrap="word", font=bold_font,bg="red",fg="white",width=80, height=15)
        self.baseview_button = tk.Button(self,text="bASe ViEw",command=self.reset_to_baseview)

        # Create a button to clear the SHOW window. 
        def tw_del(): 
            return self.text_widget.delete(1.0, tk.END)
        self.clear_button = tk.Button(self,text="CLeAr ScReeN",command=tw_del)

        self.fullreset_button = tk.Button(self,text="fuLL REseT",command=self.fullreset_clp)
        self.send_button = tk.Button(self, text="SeND cMd", command=self.send_cmd__UI)
        self.help_button = tk.Button(self, text="heLP", command=self.switch_HELP_display)

        self.text_widget2 = tk.Text(self, wrap="word", width=80, height=15)

    def set_primary_window_details(self): 

        self.open_button.grid(row=0,column=0)

        # Create a Text widget to display the content
        self.text_widget.grid()
        self.baseview_button.grid(row=0,column=10)
        self.clear_button.grid(row=14,column=0) 
        self.fullreset_button.grid(row=14,column=10)
        self.text_widget2.grid()

        # Create a button to send written commands 
        self.send_button.grid()

        # Create a help button for user-friendliness  
        self.help_button.grid(row=24,column=10)

    def switch_HELP_display(self): 
        if self.is_helpwindow:
            self.text_widget3.grid_forget() 
            self.goback_button.grid_forget()
            self.set_primary_window_details()
            self.is_helpwindow = not self.is_helpwindow
            return

        # hide all the primary display's buttons
        self.text_widget.grid_forget() 
        self.text_widget2.grid_forget() 

        self.open_button.grid_forget() 
        self.baseview_button.grid_forget() 
        self.clear_button.grid_forget() 
        self.fullreset_button.grid_forget() 
        self.send_button.grid_forget() 
        self.help_button.grid_forget() 

        self.is_helpwindow = not self.is_helpwindow

        # open up text widget 3 
        if type(self.text_widget3) == type(None): 
            self.init_help_display()

        self.text_widget3.grid() 
        self.goback_button.grid(row=30,column=0)
        self.reload_entire_help_display()

        return

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
            if type(x) == list: x = np.array(x) 
            
            self.text_widget.insert(tk.END, "SHOW {}\n{}\n".format(q[1],str(x)))

    def fullreset_clp(self): 
        self.clp = CommLangParser(None) 
        self.reset_to_baseview()

    def init_help_display(self): 
        bold_font = font.Font(family="Arial", size=12, weight="bold")
        self.text_widget3 = tk.Text(self, wrap="word", font=bold_font,bg="purple",\
            fg="yellow",width=86, height=31.5) 

        self.goback_button = tk.Button(self, text="GO bacK", command=self.switch_HELP_display)

    def reload_entire_help_display(self): 
        s = str(self.clgp) 
        self.text_widget3.insert(tk.END,s)
        return

def run_sb_app(): 
    app = SBApplication()
    app.master.title('SeQbUiLd WiNdOw')  
    app.mainloop()