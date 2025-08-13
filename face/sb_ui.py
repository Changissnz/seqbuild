import tkinter as tk
from tkinter import filedialog

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
        self.text_widget = tk.Text(self, wrap="word", width=40, height=10)
        self.text_widget.grid()

        # Create a button to open the file
        self.open_button = tk.Button(self, text="Open File", command=self.open_file__UI)
        self.open_button.grid()

    def open_file__UI(self): 
        file_path = filedialog.askopenfilename(
            title="Select a Text File", \
            filetypes=[("Text files", "*.txt")])
        if file_path:
            with open(file_path, 'r') as file:
                content = file.read()
                self.text_widget.delete(1.0, tk.END)  # Clear previous content
                self.text_widget.insert(tk.END, content)

app = Application()
app.master.title('SeQbUiLd WiNdOw')  
app.mainloop()