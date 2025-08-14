import tkinter as tk
from tkinter import filedialog

"""
the main Tkinter application class for seqbuild user 
interface (SB UI). 
"""
class Application(tk.Frame):

    def __init__(self, master=None):
        tk.Frame.__init__(self, master)
        self.grid()
        self.createWidgets()

    def createWidgets(self):
        self.quitButton = tk.Button(self, text='Quit',
            command=self.destroy)
        self.quitButton.grid()

        # Create a Text widget to display the content
        self.text_widget = tk.Text(self, wrap="word", width=80, height=15)
        self.text_widget.grid()

        # Create a button to open the file
        self.open_button = tk.Button(self, text="Open File", command=self.open_file__UI)
        self.open_button.grid()

        # Create another Text widget to write info. 
        self.text_widget2 = tk.Text(self, wrap="word", width=80, height=15)
        self.text_widget2.grid()

        # Create a button to send written commands 
        self.open_button2 = tk.Button(self, text="Send Cmd", command=self.send_cmd__UI)
        self.open_button2.grid()


    def open_file__UI(self): 
        file_path = filedialog.askopenfilename(
            title="Select a Text File", \
            filetypes=[("Text files", "*.txt")])
        if file_path:
            with open(file_path, 'r') as file:
                content = file.read()
                self.text_widget.delete(1.0, tk.END)  # Clear previous content
                self.text_widget.insert(tk.END, content)

    # TODO: route this output to a <CommLangParser> 
    def send_cmd__UI(self): 
        cmd = self.text_widget2.get("1.0", tk.END)
        self.text_widget2.delete("1.0", tk.END)
        return cmd 

app = Application()
app.master.title('SeQbUiLd WiNdOw')  
app.mainloop()