import sys
import numpy as np
import tkinter as tk
from tkinter import messagebox, ttk 
from PIL import Image, ImageTk
import matplotlib.pyplot as plt
from scipy.signal import sawtooth

def schematic_display(image_path, parent):
    original_image = Image.open(image_path)

    max_width = 600

    frame = tk.Frame(parent, bg="#852323", pady=20)
    frame.pack()

    aspect_ratio = original_image.height / original_image.width
    initial_width = max_width
    initial_height = int(initial_width * aspect_ratio)
    resized = original_image.resize((initial_width, initial_height), Image.Resampling.LANCZOS)
    photo = ImageTk.PhotoImage(resized)

    label = tk.Label(frame, image=photo, bg="#852323")
    label.image = photo
    label.pack()

    #Skalowanie
    def resize_image(event):
        new_width = min(event.width, max_width)
        new_height = int(new_width * aspect_ratio)
        resized = original_image.resize((new_width, new_height), Image.Resampling.LANCZOS)
        photo = ImageTk.PhotoImage(resized)
        label.config(image=photo)
        label.image = photo

    frame.bind("<Configure>", resize_image)


def aktualizacja_zmiennych():
    try:
        J1 = float(J1_entry.get())
        n1 = float(n1_entry.get())
        n2 = float(n2_entry.get())
        J2 = float(J2_entry.get())
        b = float(b_entry.get())
        k = float(k_entry.get())
    except ValueError:
        messagebox.showerror(title='Błędne dane',message='Błędne dane - wprowadź poprawne wartości')

# simulation functions:
# derrivatives:
def x1_derrivative(x2):
    x1=x2
    return x1

def x2_derrivative(J1,Tm,k,n,x1,x3):
    x2=(1/J1)*(Tm-k*(x1-x3/n))
    return x2

def x3_derrivative(x4):
    x3=x4
    return x3

def x4_derrivative(J2,k,b,n,x1,x3,x4):
    x4=(1/J2)*(k*n*(x1-x3/n)-b*x4)
    return x4


def Euler(J1,J2,n,k,b,sygnal,T,czas_values):
    x0=np.zeros(4)
    tn=0
    h=T/10000
    xn=x0
    x_values=np.zeros((4,len(czas_values)))
    for i in range(0, len(czas_values)):
        Tm=sygnal[i]
        #print(czas_values[i],xn)
        tn_plus_one=tn+h
        x_values[0,i]=xn[0]+h*x1_derrivative(xn[1])                        # x1 theta1
        x_values[1,i]=xn[1]+h*x2_derrivative(J1,Tm,k,n,xn[0],xn[2])        # x2 omega1
        x_values[2,i]=xn[2]+h*x3_derrivative(xn[3])                        # x3 theta2
        x_values[3,i]=xn[3]+h*x4_derrivative(J2,k,b,n,xn[0],xn[2],xn[3])   # x4 omega2
        tn=tn_plus_one
        xn=np.transpose(x_values)[i]
    #print(czas_values,x_values) 
    return x_values
    
def RK4(J1,J2,n,k,b,sygnal,T,czas_values): # runge-kutta
    x0 = np.zeros(4)
    h = T / 10000
    xn = x0
    x_values = np.zeros((4, len(czas_values)))
    for i in range(len(czas_values)):
        Tm = sygnal[i]

        def f(x, Tm):
            dx1 = x1_derrivative(x[1])
            dx2 = x2_derrivative(J1,Tm,k,n,x[0],x[2])
            dx3 = x3_derrivative(x[3])
            dx4 = x4_derrivative(J2,k,b,n,x[0],x[2],x[3])
            return np.array([dx1,dx2,dx3,dx4])

        k1 = f(xn, Tm)
        k2 = f(xn + 0.5 * h * k1, Tm)
        k3 = f(xn + 0.5 * h * k2, Tm)
        k4 = f(xn + h * k3, Tm)
        xn = xn + (h / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)
        x_values[:, i] = xn
    return x_values

def generowanie_odpowiedzi():
    try:
        J1 = float(J1_entry.get())
        n1 = float(n1_entry.get())
        n2 = float(n2_entry.get())
        J2 = float(J2_entry.get())
        b = float(b_entry.get())
        k = float(k_entry.get())
        n = float(n1 / n2)

        opcja = x.get()
        T = float(czas_entry.get())
        amplituda = float(amplituda_entry.get())

        fig, axs = plt.subplots(2, 2, figsize=(10, 8))
        czas_values = np.linspace(0, T, 10000)
        if opcja == 1:
            czestotliwosc = float(czestotliwosc_entry.get())
            faza = float(faza_entry.get())
            sygnal = amplituda * np.sin(2 * np.pi * czestotliwosc * czas_values + np.deg2rad(faza))
        elif opcja == 2:
            sygnal = amplituda * np.heaviside(czas_values, 1)
        elif opcja == 3:
            czestotliwosc = float(czestotliwosc_entry.get())
            sygnal = amplituda * sawtooth(2 * np.pi * czestotliwosc * czas_values, 0.5)
        elif opcja == 4:
            czas_impulsu = float(czas_impulsu_entry.get())
            sygnal = amplituda * (czas_values < czas_impulsu)

        # obliczanie odpowiedzi oboma metodami
        x_euler = Euler(J1, J2, n, k, b, sygnal, T, czas_values)
        x_rk4 = RK4(J1, J2, n, k, b, sygnal, T, czas_values)

        try:
            # theta1
            ax = axs[0][0]
            ax.grid(True, linestyle='-.')
            ax.tick_params(labelcolor='r', labelsize='medium', width=3)
            ax.plot(czas_values, x_euler[0], label='Euler')
            ax.plot(czas_values, x_rk4[0], label='RK4', linestyle='dashed')
            ax.set_title('\u03B8\u2081 [rad]')
            ax.set_xlabel("Czas [s]")
            ax.set_ylabel("Wartość")
            ax.legend()

            # theta2
            ax = axs[0][1]
            ax.grid(True, linestyle='-.')
            ax.tick_params(labelcolor='r', labelsize='medium', width=3)
            ax.plot(czas_values, x_euler[2], label='Euler')
            ax.plot(czas_values, x_rk4[2], label='RK4', linestyle='dashed')
            ax.set_title('\u03B8\u2082 [rad]')
            ax.set_xlabel("Czas [s]")
            ax.set_ylabel("Wartość")
            ax.legend()

            # omega1
            ax = axs[1][0]
            ax.grid(True, linestyle='-.')
            ax.tick_params(labelcolor='r', labelsize='medium', width=3)
            ax.plot(czas_values, x_euler[1], label='Euler')
            ax.plot(czas_values, x_rk4[1], label='RK4', linestyle='dashed')
            ax.set_title('\u03C9\u2081 [rad]')
            ax.set_xlabel("Czas [s]")
            ax.set_ylabel("Wartość")
            ax.legend()

            # omega2
            ax = axs[1][1]
            ax.grid(True, linestyle='-.')
            ax.tick_params(labelcolor='r', labelsize='medium', width=3)
            ax.plot(czas_values, x_euler[3], label='Euler')
            ax.plot(czas_values, x_rk4[3], label='RK4', linestyle='dashed')
            ax.set_title('\u03C9\u2082 [rad]')
            ax.set_xlabel("Czas [s]")
            ax.set_ylabel("Wartość")
            ax.legend()

            fig.align_labels()
            fig.align_titles()
            fig.suptitle("Porównanie metod: Euler vs RK4")
            plt.show()
        except ValueError:
            messagebox.showerror("Błąd symulacji", "Błąd symulacji - błąd podczas generowania odpowiedzi układu.")

    except ValueError:
        messagebox.showerror("Błąd symulacji", "Błąd symulacji - wprowadź poprawnie wszystkie liczby")


def main():
    okno = tk.Tk()
    okno.title("Projekt MMM - model układu mechanicznego - 197697 , 197763")
    width=int(okno.winfo_screenwidth()/2)
    height=int(okno.winfo_screenheight()/2)
    w_offset=int(width/2)
    h_offset=int(height/3)
    okno.geometry('{}x{}+{}+{}'.format(width,height,w_offset,h_offset)) # wymiary okna
    okno.call('wm','iconphoto',okno._w,tk.PhotoImage(file='icon.png')) # ikona na oknie
    okno.config(background="#852323")
 
    etykieta = tk.Label(okno, text="Symulator układu mechanicznego",bg='#852323',fg='white',font=("bold"))
    etykieta.pack()

    schematic_display("./scheme.png", okno) # zdjęcie schematu

    global J1_entry, n1_entry, n2_entry, J2_entry, b_entry, k_entry # zmienne do pobierania danych

    canvas = tk.Canvas(okno)
    scrollbar = ttk.Scrollbar(okno, orient="vertical", command=canvas.yview)
    scrollable_frame = ttk.Frame(canvas)

    canvas_frame = canvas.create_window((0, 0), window=scrollable_frame, anchor="nw") # utworzenie okna w canvasie

    def on_frame_configure(event):
        canvas.configure(scrollregion=canvas.bbox("all"))  # automatyczne dopasowanie scrolla

    scrollable_frame.bind("<Configure>", on_frame_configure)

    canvas.configure(yscrollcommand=scrollbar.set)
    canvas.pack(side="left", fill="both", expand=True)
    scrollbar.pack(side="right", fill="y")

    def _on_mousewheel(event): 
        canvas.yview_scroll(int(-1*(event.delta/120)), "units")  # umożliwia używanie scrolla do przewijania
    canvas.bind_all("<MouseWheel>", _on_mousewheel)

    input_frame = tk.Frame(scrollable_frame)
    input_frame.grid(row=0, column=0, padx=20, pady=20, sticky='nw')

    tk.Label(input_frame, text="J1:", font=("Fira Code", 14)).grid(row=0, column=0, sticky='e') # Pola do wprowadzania zmiennych
    J1_entry = tk.Entry(input_frame)
    J1_entry.grid(row=0, column=1, sticky='w')

    tk.Label(input_frame, text="n1:", font=("Fira Code", 14)).grid(row=1, column=0, sticky='e')
    n1_entry = tk.Entry(input_frame)
    n1_entry.grid(row=1, column=1, sticky='w')

    tk.Label(input_frame, text="n2:", font=("Fira Code", 14)).grid(row=2, column=0, sticky='e')
    n2_entry = tk.Entry(input_frame)
    n2_entry.grid(row=2, column=1, sticky='w')

    tk.Label(input_frame, text="J2:", font=("Fira Code", 14)).grid(row=3, column=0, sticky='e')
    J2_entry = tk.Entry(input_frame)
    J2_entry.grid(row=3, column=1, sticky='w')

    tk.Label(input_frame, text="b:", font=("Fira Code", 14)).grid(row=4, column=0, sticky='e')
    b_entry = tk.Entry(input_frame)
    b_entry.grid(row=4, column=1, sticky='w')

    tk.Label(input_frame, text="k:", font=("Fira Code", 14)).grid(row=5, column=0, sticky='e')
    k_entry = tk.Entry(input_frame)
    k_entry.grid(row=5, column=1, sticky='w')

    update_button = tk.Button(input_frame, text="Zatwierdź", command=aktualizacja_zmiennych)
    update_button.grid(row=6, column=2)

    global czas_entry, amplituda_entry, czestotliwosc_entry, faza_entry, x, czas_impulsu_entry
    czas_label = tk.Label(input_frame, text="Czas pobudzenia [s]:", font=("Fira Code", 10))
    czas_label.grid(row=7, column=1, sticky='e', padx=10, pady=10)
    czas_entry = tk.Entry(input_frame)
    czas_entry.grid(row=7, column=2, sticky='w', padx=10, pady=10)

    tk.Label(input_frame, text="Rodzaj pobudzenia:", font=("Fira Code", 10)).grid(row=8, column=1, sticky='e', padx=10, pady=10)
    x = tk.IntVar(value=2)
    tk.Radiobutton(input_frame, text="Sinus", variable=x, value=1, font=("Fira Code", 10)).grid(row=8, column=2, sticky='w', padx=10, pady=5)
    tk.Radiobutton(input_frame, text="Skok jednostkowy", variable=x, value=2, font=("Fira Code", 10)).grid(row=9, column=2, sticky='w', padx=10, pady=5)
    tk.Radiobutton(input_frame, text="Trojkatne", variable=x, value=3, font=("Fira Code", 10)).grid(row=10, column=2, sticky='w', padx=10, pady=5)
    tk.Radiobutton(input_frame, text="Prostokątne", variable=x, value=4, font=("Fira Code", 10)).grid(row=11, column=2, sticky='w', padx=10, pady=5)

    amplituda_label = tk.Label(input_frame, text="Amplituda:", font=("Fira Code", 10))
    amplituda_label.grid(row=12, column=1, sticky='e', padx=10, pady=10)
    amplituda_entry = tk.Entry(input_frame)
    amplituda_entry.grid(row=12, column=2, sticky='w', padx=10, pady=10)

    czestotliwosc_label = tk.Label(input_frame, text="Częstotliwość [Hz]:", font=("Fira Code", 10))
    czestotliwosc_label.grid(row=13, column=1, sticky='e', padx=10, pady=10)
    czestotliwosc_entry = tk.Entry(input_frame)
    czestotliwosc_entry.grid(row=13, column=2, sticky='w', padx=10, pady=10)
    czestotliwosc_label = tk.Label(input_frame, text="(tylko sinusoida i trojkatny)", font=("Fira Code", 10))
    czestotliwosc_label.grid(row=13, column=3, sticky='e', padx=5, pady=10)

    faza_label = tk.Label(input_frame, text="Faza w stopniach:", font=("Fira Code", 10))
    faza_label.grid(row=14, column=1, sticky='e', padx=10, pady=10)
    faza_entry = tk.Entry(input_frame)
    faza_entry.grid(row=14, column=2, sticky='w', padx=10, pady=10)
    faza_label = tk.Label(input_frame, text="(tylko sinusoida)", font=("Fira Code", 10))
    faza_label.grid(row=14, column=3, sticky='e', padx=5, pady=10)

    czas_impulsu_label = tk.Label(input_frame, text="Czas trwania impulsu [s]:", font=("Fira Code", 10))
    czas_impulsu_label.grid(row=15, column=1, sticky='e', padx=10, pady=10)
    czas_impulsu_entry = tk.Entry(input_frame)
    czas_impulsu_entry.grid(row=15, column=2, sticky='w', padx=10, pady=10)

    plot_button = tk.Button(input_frame, text="Generuj odpowiedź", command=generowanie_odpowiedzi)
    plot_button.grid(row=16, column=2, sticky='w', padx=10, pady=10)

    okno.mainloop()

if __name__ == "__main__":
    main()
sys.exit(0)