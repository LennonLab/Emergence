import Tkinter as tk
from Tkinter import *
import time


class Alien(object):
    def __init__(self, canvas, *args, **kwargs):
        self.canvas = canvas
        self.id = canvas.create_oval(*args, **kwargs)
        self.vx = 5
        self.vy = 0
    def move(self):
        x1, y1, x2, y2 = self.canvas.bbox(self.id)
        if x2 > 400: self.vx = -5
        if x1 < 0: self.vx = 5
        self.canvas.move(self.id, self.vx, self.vy)

class App(object):
     def __init__(self, master, **kwargs):
        self.master = master
        self.canvas = tk.Canvas(self.master, width = 400, height = 400)
        self.canvas.pack()
        self.aliens = [
            Alien(self.canvas, 20, 260, 120, 360, outline = 'white', fill = 'blue'),
            Alien(self.canvas, 2, 2, 40, 40, outline = 'white', fill = 'red'),
            ]
        self.canvas.pack()
        self.master.after(0, self.animation)

     def animation(self):
         for alien in self.aliens:
             alien.move()
         self.master.after(12, self.animation)

root = tk.Tk()
app = App(root)
root.mainloop()
