import tkinter as tk
from tkinter import filedialog,font
from .comm_lang import * 

"""
the main Tkinter application class for seqbuild user 
interface (SB UI). 
"""
class Application(tk.Frame):

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

        # Create a button to clear the SHOW window. 
        def tw_del(): 
            return self.text_widget.delete(1.0, tk.END)
        self.clear_button = tk.Button(self,text="CLeAr ScReeN",command=tw_del)
        self.clear_button.grid(row=14,column=0) 

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
            self.text_widget.delete(1.0, tk.END)
            self.text_widget.insert(tk.END, str(self.clp))  
            return

    def send_cmd__UI(self): 
        cmd = self.text_widget2.get("1.0", tk.END)
        self.text_widget2.delete("1.0", tk.END)
        cmd_lines = cmd.split("\n")
        cmd_lines = [cline.strip() for cline in cmd_lines]
        self.clp.load_cmd_lines(cmd_lines)
        self.clp.process_cmdlines()

        self.text_widget.delete(1.0, tk.END)  # Clear previous content
        self.text_widget.insert(tk.END, str(self.clp))  
        return cmd 

app = Application()
app.master.title('SeQbUiLd WiNdOw')  
app.mainloop()