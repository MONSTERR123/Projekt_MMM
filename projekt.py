import sys
import traceback
import os
import numpy as np
import tkinter as tk
import scipy as scp
from tkinter import messagebox, ttk 
from PIL import Image, ImageTk
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg,
NavigationToolbar2Tk)

from scipy.signal import sawtooth

old_figs = []
figure_canvases = []
selected_exc_plot = None
exc_fig_czas = None
exc_fig_time = None
exc_canvas = None

def on_closing(okno):
    for fig in plt.get_fignums():
        plt.close(fig)
    okno.quit()
    okno.destroy()
    sys.exit(0)

def resource_path(relative_path):
    """ Get absolute path to resource, works for dev and for PyInstaller """
    try:
        base_path = sys._MEIPASS  # PyInstaller sets this in the temp folder
    except Exception:
        base_path = os.path.abspath("D:/Kuba/sem4/MMM_projekt")

    return os.path.join(base_path, relative_path)

def schematic_display(image_path, parent):
    original_image = Image.open(image_path)
    max_width = 600

    frame = tk.Frame(parent, bg="#852323", pady=20)
    frame.pack()

    aspect_ratio = original_image.height / original_image.width
    initial_width = max_width
    initial_height = int(initial_width * aspect_ratio)
    resized = original_image.resize((initial_width, initial_height), Image.Resampling.LANCZOS)

   
    persistent = {}
    persistent['photo'] = ImageTk.PhotoImage(resized)

    label = tk.Label(frame, image=persistent['photo'], bg="#852323")
    label.image = persistent['photo']  
    label.pack()

    # Resize function to adjust image when the window is resized
    def resize_image(event):
        new_width = min(event.width, max_width)
        new_height = int(new_width * aspect_ratio)
        resized_img = original_image.resize((new_width, new_height), Image.Resampling.LANCZOS)
        persistent['photo'] = ImageTk.PhotoImage(resized_img)
        label.config(image=persistent['photo'])
        label.image = persistent['photo']

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


def Euler(J1,J2,n,k,b,sygnal,T,czas_values,time_step):
    x0=np.zeros(4)
    tn=0
    h=T/time_step
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
    
def RK4(J1, J2, n, k, b, sygnal, T, czas_values, time_step): # runge kutta
    x0 = np.zeros(4)
    h = T / time_step
    xn = x0
    x_values = np.zeros((4, len(czas_values)))
    for i in range(len(czas_values)):
        Tm=sygnal[i]
        def f(x, Tm):
            dx1 = x1_derrivative(x[1])
            dx2 = x2_derrivative(J1, Tm, k, n, x[0], x[2])
            dx3 = x3_derrivative(x[3])
            dx4 = x4_derrivative(J2, k, b, n, x[0], x[2], x[3])
            return np.array([dx1, dx2, dx3, dx4])

        k1 = f(xn, Tm)
        k2 = f(xn + 0.5 * h * k1, Tm)
        k3 = f(xn + 0.5 * h * k2, Tm)
        k4 = f(xn + h * k3, Tm)
        xn = xn + (h / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)
        x_values[:, i] = xn
    return x_values

def pobudzenie(opcja,amplituda,czas_values):
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
    elif opcja == 5:
        sygnal = scp.signal.unit_impulse(len(czas_values))
        margin=0.01
    return sygnal

def update_exc_plot():
        global exc_canvas
        if exc_canvas:
            exc_canvas.get_tk_widget().destroy()
        if selected_exc_plot.get() == 1:
            fig = exc_fig_czas
        else:
            fig = exc_fig_time
        if fig:
            exc_canvas = FigureCanvasTkAgg(fig, master=frame1)
            exc_canvas.draw()
            exc_canvas.get_tk_widget().pack(fill="both", expand=True)
            figure_canvases.append(exc_canvas)

def generowanie_odpowiedzi():
    try:
        global old_figs
        for fig in old_figs:
            plt.close(fig)
        old_figs = []

        J1 = float(J1_entry.get())
        n1 = float(n1_entry.get())
        n2 = float(n2_entry.get())
        J2 = float(J2_entry.get())
        b = float(b_entry.get())
        k = float(k_entry.get())

        time_mode = time_mode_var.get()

        opcja = x.get()
        T = float(czas_entry.get())

        if step_entry.get() <= "0":
            messagebox.showerror("Błąd symulacji", "Ilość próbek musi być większa od 0.")
            return

        # wybór rozmiaru kroku lub liczby próbek 
        if time_mode == 1:
            time_step = int(step_entry.get())
            step_size = float(T / time_step)
        else:
            step_size = float(step_entry.get())
            time_step = int(T / step_size)
        
        # sprawdzenie poprawności danych wejściowych
        if step_entry.get() == "":
            messagebox.showerror("Błąd symulacji", "Podaj ilość próbek lub rozmiar kroku.")
            return
        
        # Sprawdzenie czy liczba próbek nie jest zbyt duża 
        if time_mode == 1:
            if time_step > 1e5:
                messagebox.showerror("Błąd symulacji", "Ilość próbek jest zbyt duża (max 100000).")
                return
        else:
            if step_size < 1e-5:
                messagebox.showerror("Błąd symulacji", "Rozmiar kroku jest zbyt mały.")
                return

        if T <= 0:
            messagebox.showerror("Błąd symulacji", "Czas pobudzenia musi być dodatni.")
            return
        
        if step_size <= 0:
            messagebox.showerror("Błąd symulacji", "Rozmiar kroku musi być dodatni.")
            return
        elif step_size > T:
            messagebox.showerror("Błąd symulacji", "Rozmiar kroku musi być mniejszy niż czas pobudzenia.")
            return

        if J1 <= 0 or J2 <= 0:
            messagebox.showerror("Błąd symulacji", "J1 i J2 nie mogą być zerowe.")
            return
        
        if n1 <= 0 or n2 <= 0:
            messagebox.showerror("Błąd symulacji", "n1 i n2 nie mogą być zerowe.")
            return
        elif n1.is_integer() == False or n2.is_integer() == False:
            messagebox.showerror("Błąd symulacji", "n1 i n2 muszą być liczbami całkowitymi.")
            return
        
        # Sprawdzenie amplitudy 
        if opcja != 5:
            if amplituda_entry.get() == "":
                messagebox.showerror("Błąd symulacji", "Podaj wartość amplitudy.")
                return
            try:
                _ = float(amplituda_entry.get())
            except Exception:
                messagebox.showerror("Błąd symulacji", "Amplituda musi być liczbą.")
                return

        # Sprawdzenie częstotliwości i fazy 
        if opcja in [1, 3]:
            if czestotliwosc_entry.get() == "":
                messagebox.showerror("Błąd symulacji", "Podaj częstotliwość.")
                return
            try:
                _ = float(czestotliwosc_entry.get())
            except Exception:
                messagebox.showerror("Błąd symulacji", "Częstotliwość musi być liczbą.")
                return
            if step_size >= 0.5 / float(czestotliwosc_entry.get()):
                messagebox.showerror("Błąd symulacji", "Rozmiar kroku musi być mniejszy niż okres sygnału.")
                return

        
        n = float(n1 / n2)
        
        margin = 0
        if opcja == 5:
            amplituda = 1.0  # domyślna amplituda dla Diraca
        else:
            amplituda = float(amplituda_entry.get())


        # dobieranie ilości próbek 
        min_samples = 2000
        samples = min_samples
        if opcja in [1, 3]:  # Sinus lub trójkątny
            try:
                freq = float(czestotliwosc_entry.get())
                if freq > 0:
                    samples = max(int(20 * freq * T), min_samples)
            except Exception:
                samples = min_samples
        czas_values = np.linspace(0, T, samples)
        # ----------------------------------------------------
        if time_step <= 0:
            messagebox.showerror("Błąd symulacji", "Ilość próbek musi być dodatnia.")
            return
        # pobudzenie sygnału
        
        u_t=pobudzenie(opcja,amplituda,czas_values)

        # pobudzenie
        fig1, ax1 = plt.subplots(figsize=(6, 4))
        ax1.plot(czas_values, u_t)
        ax1.grid(True, linestyle='-.')
        ax1.tick_params(labelcolor='r', labelsize='medium', width=3)
        ax1.set_xlabel("Czas [s]")
        ax1.set_ylabel("T\u2098 [Nm]")
        ax1.margins(x=margin)
        fig1.suptitle("Pobudzenie T\u2098(t)")

        time_values = np.linspace(0, T, time_step+1)
        sygnal=pobudzenie(opcja,amplituda,time_values)
        
        # obliczanie odpowiedzi oboma metodami
        x_euler = Euler(J1, J2, n, k, b, sygnal, T, time_values, time_step)
        x_rk4 = RK4(J1, J2, n, k, b, sygnal, T, time_values, time_step)

        # theta & omega plots
        fig2, axs2 = plt.subplots(2, 2, figsize=(10, 8))  # Euler
        fig3, axs3 = plt.subplots(2, 2, figsize=(10, 8))  # RK4
        fig4, axs4 = plt.subplots(2, 2, figsize=(10, 8))  # both methods

        fig_list = [fig1, fig2, fig3, fig4]
        axs_list = [None, axs2, axs3, axs4]
        frame_list = [frame1, frame2, frame3, frame4]

        def plot_subplot(ax, data1, data2, label, title):
            if j == 1:  # both methods
                ax.plot(time_values, data1, label='Euler')
                ax.plot(time_values, data2, label='RK4', linestyle='dashed')
                ax.legend()
            elif j == 2: # Euler
                ax.plot(time_values, data1)
            elif j == 3: # RK4
                ax.plot(time_values, data2)
            ax.set_title(title)
            ax.set_xlabel("Czas [s]")
            ax.set_ylabel("Wartość")
            ax.grid(True, linestyle='-.')
            ax.margins(x=0)

        for j in range(4):
            if j == 0:
                for widget in frame_list[j].winfo_children():
                    if isinstance(widget, FigureCanvasTkAgg):
                        widget.get_tk_widget().destroy()
                continue
            else:
                for widget in frame_list[j].winfo_children():
                    widget.destroy()

            axs = axs_list[j]
            method = ['porównanie metod', 'metoda Eulera', 'metoda RK4'][j - 1]
            fig_list[j].suptitle(f"{method}")
            fig_list[j].subplots_adjust(hspace=0.5)

            plot_subplot(axs[0][0], x_euler[0], x_rk4[0], '\u03B8\u2081', '\u03B8\u2081 [rad]')
            plot_subplot(axs[0][1], x_euler[2], x_rk4[2], '\u03B8\u2082', '\u03B8\u2082 [rad]')
            plot_subplot(axs[1][0], x_euler[1], x_rk4[1], '\u03C9\u2081', '\u03C9\u2081 [rad/s]')
            plot_subplot(axs[1][1], x_euler[3], x_rk4[3], '\u03C9\u2082', '\u03C9\u2082 [rad/s]')

            # New plot
            canvas = FigureCanvasTkAgg(fig_list[j], master=frame_list[j])
            canvas.draw()
            canvas.get_tk_widget().pack(fill="both", expand=True)
            figure_canvases.append(canvas)
            

        # Wykres pobudzenia dla czas_values i time_values 
        global exc_fig_czas, exc_fig_time, exc_canvas

        # Wykres pobudzenia dla czas_values
        exc_fig_czas, ax_exc_czas = plt.subplots(figsize=(6, 4))
        ax_exc_czas.plot(czas_values, u_t)
        ax_exc_czas.grid(True, linestyle='-.')
        ax_exc_czas.tick_params(labelcolor='r', labelsize='medium', width=3)
        ax_exc_czas.set_xlabel("Czas [s]")
        ax_exc_czas.set_ylabel("T\u2098 [Nm]")
        ax_exc_czas.margins(x=margin)
        exc_fig_czas.suptitle("Pobudzenie T\u2098(t) rozdzielczość optymalna")

        # Wykres pobudzenia dla time_values
        exc_fig_time, ax_exc_time = plt.subplots(figsize=(6, 4))
        ax_exc_time.plot(time_values, sygnal)
        ax_exc_time.grid(True, linestyle='-.')
        ax_exc_time.tick_params(labelcolor='r', labelsize='medium', width=3)
        ax_exc_time.set_xlabel("Czas [s]")
        ax_exc_time.set_ylabel("T\u2098 [Nm]")
        ax_exc_time.margins(x=margin)
        exc_fig_time.suptitle("Pobudzenie T\u2098(t) rozdzielczość zgodna z symulacją")

        old_figs.extend([fig1, fig2, fig3, fig4, exc_fig_czas, exc_fig_time])

        update_exc_plot()

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
    icon_image = tk.PhotoImage(file=resource_path("icon.png"))  
    okno.call('wm', 'iconphoto', okno._w, icon_image) # ikona na oknie
    notebook=ttk.Notebook(okno)
    notebook.pack(expand=True,fill="both")
    global frame0, frame1, frame2, frame3, frame4
    frame0=tk.Frame(notebook)
    frame1=tk.Frame(notebook)
    frame2=tk.Frame(notebook)   
    frame3=tk.Frame(notebook)
    frame4=tk.Frame(notebook)
    frame0.config(background="#852323")
 
    etykieta = tk.Label(frame0, text="Symulator układu mechanicznego",bg='#852323',fg='white',font=("bold"))
    etykieta.pack()

    schematic_display(resource_path("scheme.png"), frame0)# zdjęcie schematu

    global J1_entry, n1_entry, n2_entry, J2_entry, b_entry, k_entry # zmienne do pobierania danych

    canvas = tk.Canvas(frame0)
    scrollbar = ttk.Scrollbar(frame0, orient="vertical", command=canvas.yview)
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
    update_button.grid(row=6, column=2,pady=15)

    global step_entry, czas_entry, amplituda_entry, czestotliwosc_entry, faza_entry, x, czas_impulsu_entry
    czas_label = tk.Label(input_frame, text="Czas pobudzenia [s]:", font=("Fira Code", 10))
    czas_label.grid(row=7, column=1, sticky='e', padx=10, pady=10)
    czas_entry = tk.Entry(input_frame)
    czas_entry.grid(row=7, column=2, sticky='w', padx=10, pady=10)

    # Wybór sposobu podania czasu dyskretyzacji 
    global time_mode_var, step_label, step_entry
    time_mode_var = tk.IntVar(value=1)  # 1 - liczba próbek, 2 - rozmiar kroku

    def update_step_label():
        if time_mode_var.get() == 1:
            step_label.config(text="Ilość próbek:")
        else:
            step_label.config(text="Rozmiar kroku [s]:")

    tk.Radiobutton(input_frame, text="Liczba próbek", variable=time_mode_var, value=1, font=("Fira Code", 10), command=update_step_label).grid(row=8, column=4, sticky='w')
    tk.Radiobutton(input_frame, text="Rozmiar kroku [s]", variable=time_mode_var, value=2, font=("Fira Code", 10), command=update_step_label).grid(row=8, column=5, sticky='w')

    step_label = tk.Label(input_frame, text="Ilość próbek:", font=("Fira Code", 10))
    step_label.grid(row=7, column=3, sticky='e', padx=10, pady=10)
    step_entry = tk.Entry(input_frame)
    step_entry.grid(row=7, column=4, sticky='w', padx=10, pady=10)

    tk.Label(input_frame, text="Rodzaj pobudzenia:", font=("Fira Code", 10)).grid(row=8, column=1, sticky='e', padx=10, pady=10)
    x = tk.IntVar(value=2)
    tk.Radiobutton(input_frame, text="Sinus", variable=x, value=1, font=("Fira Code", 10)).grid(row=8, column=2, sticky='w', padx=10, pady=5)
    tk.Radiobutton(input_frame, text="Skok jednostkowy", variable=x, value=2, font=("Fira Code", 10)).grid(row=9, column=2, sticky='w', padx=10, pady=5)
    tk.Radiobutton(input_frame, text="Trojkatne", variable=x, value=3, font=("Fira Code", 10)).grid(row=10, column=2, sticky='w', padx=10, pady=5)
    tk.Radiobutton(input_frame, text="Prostokątne", variable=x, value=4, font=("Fira Code", 10)).grid(row=11, column=2, sticky='w', padx=10, pady=5)
    tk.Radiobutton(input_frame, text="Impuls Diraca", variable=x, value=5, font=("Fira Code", 10)).grid(row=12, column=2, sticky='w', padx=10, pady=5)

    amplituda_label = tk.Label(input_frame, text="Amplituda:", font=("Fira Code", 10))
    amplituda_label.grid(row=13, column=1, sticky='e', padx=10, pady=10)
    amplituda_entry = tk.Entry(input_frame)
    amplituda_entry.grid(row=13, column=2, sticky='w', padx=10, pady=10)

    czestotliwosc_label = tk.Label(input_frame, text="Częstotliwość [Hz]:", font=("Fira Code", 10))
    czestotliwosc_label.grid(row=14, column=1, sticky='e', padx=10, pady=10)
    czestotliwosc_entry = tk.Entry(input_frame)
    czestotliwosc_entry.grid(row=14, column=2, sticky='w', padx=10, pady=10)
    czestotliwosc_label = tk.Label(input_frame, text="(tylko sinusoida i trojkatny)", font=("Fira Code", 10))
    czestotliwosc_label.grid(row=14, column=3, sticky='e', padx=5, pady=10)

    faza_label = tk.Label(input_frame, text="Faza w stopniach:", font=("Fira Code", 10))
    faza_label.grid(row=15, column=1, sticky='e', padx=10, pady=10)
    faza_entry = tk.Entry(input_frame)
    faza_entry.grid(row=15, column=2, sticky='w', padx=10, pady=10)
    faza_label = tk.Label(input_frame, text="(tylko sinusoida)", font=("Fira Code", 10))
    faza_label.grid(row=15, column=3, sticky='e', padx=5, pady=10)

    czas_impulsu_label = tk.Label(input_frame, text="Czas trwania impulsu [s]:", font=("Fira Code", 10))
    czas_impulsu_label.grid(row=16, column=1, sticky='e', padx=10, pady=10)
    czas_impulsu_entry = tk.Entry(input_frame)
    czas_impulsu_entry.grid(row=16, column=2, sticky='w', padx=10, pady=10)
    faza_label = tk.Label(input_frame, text="(tylko prostokątny)", font=("Fira Code", 10))
    faza_label.grid(row=16, column=3, sticky='e', padx=5, pady=10)

    plot_button = tk.Button(input_frame, text="Generuj odpowiedź", command=generowanie_odpowiedzi)
    plot_button.grid(row=17, column=2, sticky='w', padx=10, pady=10)

    global selected_exc_plot, exc_canvas
    selected_exc_plot = tk.IntVar(value=1)  # 1 - czas_values, 2 - time_values

    tk.Label(frame1, text="Wybierz rozdzielczość wykresu pobudzenia:", font=("Fira Code", 10)).pack(pady=5)
    tk.Radiobutton(frame1, text="optymalna rozdzielczość", variable=selected_exc_plot, value=1, command=update_exc_plot).pack()
    tk.Radiobutton(frame1, text="rozdzielczość zgodna z symulacją odpowiedzi", variable=selected_exc_plot, value=2, command=update_exc_plot).pack()

    frame0.pack(fill="both",expand=True)
    frame1.pack(fill="both",expand=True)
    frame2.pack(fill="both",expand=True)
    frame3.pack(fill="both",expand=True)
    frame4.pack(fill="both",expand=True)
    notebook.add(frame0,text='parametry symulacji')
    notebook.add(frame1,text='pobudzenie')
    notebook.add(frame2,text='porównanie metod')
    notebook.add(frame3,text='metoda Eulera')
    notebook.add(frame4,text='metoda RK4')

    okno.protocol("WM_DELETE_WINDOW", lambda: on_closing(okno))
    okno.mainloop()
    sys.exit(0)

if __name__ == "__main__":
    main()
