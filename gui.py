import tkinter as tk
from tkinter import ttk, messagebox
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import numpy as np
from calcs import run_calculations, compute_poly_AI, rotate_around_point
from blade_profile import get_blade_profile

class BladeCalculatorGUI:
    def __init__(self, root):
        self.root = root
        self.root.title("Blade Calculator")
        
        # Default values
        self.default_values = {
            'E': 650000000,  # Young's modulus (N/m^3)
            'L': 0.05,       # Length (m)
            'rho': 1710,     # Density (kg/m^3)
            'twist': 0.1,    # Twist angle
            'layers': 25     # Number of layers for 3D plot
        }
        
        # Create main frame
        self.main_frame = ttk.Frame(root, padding="10")
        self.main_frame.grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))
        
        # Create input frame
        self.create_input_frame()
        
        # Create results frame
        self.create_results_frame()
        
        # Create plots frame
        self.create_plots_frame()
        
        # Generate button
        self.generate_btn = ttk.Button(self.main_frame, text="Generate", command=self.generate)
        self.generate_btn.grid(row=3, column=0, columnspan=2, pady=10)
        
        # Initialize with default values
        self.set_default_values()
        
        # Generate initial results and plots
        self.generate()
        
    def create_input_frame(self):
        input_frame = ttk.LabelFrame(self.main_frame, text="Input Parameters", padding="5")
        input_frame.grid(row=0, column=0, padx=5, pady=5, sticky=(tk.W, tk.E))
        
        # Create input fields
        self.input_vars = {}
        row = 0
        for param, default in self.default_values.items():
            ttk.Label(input_frame, text=f"{param}:").grid(row=row, column=0, sticky=tk.W)
            var = tk.StringVar(value=str(default))
            self.input_vars[param] = var
            entry = ttk.Entry(input_frame, textvariable=var)
            entry.grid(row=row, column=1, padx=5, pady=2)
            row += 1
            
    def create_results_frame(self):
        results_frame = ttk.LabelFrame(self.main_frame, text="Results", padding="5")
        results_frame.grid(row=0, column=1, padx=5, pady=5, sticky=(tk.W, tk.E))
        
        # Create results text widget with scrollbar
        self.results_text = tk.Text(results_frame, height=15, width=50)
        scrollbar = ttk.Scrollbar(results_frame, orient="vertical", command=self.results_text.yview)
        self.results_text.configure(yscrollcommand=scrollbar.set)
        
        self.results_text.grid(row=0, column=0, padx=5, pady=5, sticky=(tk.W, tk.E, tk.N, tk.S))
        scrollbar.grid(row=0, column=1, sticky=(tk.N, tk.S))
        
    def create_plots_frame(self):
        plots_frame = ttk.LabelFrame(self.main_frame, text="Visualizations", padding="5")
        plots_frame.grid(row=1, column=0, columnspan=2, padx=5, pady=5, sticky=(tk.W, tk.E))
        
        # Create matplotlib figure with subplots
        self.fig = plt.figure(figsize=(12, 8))
        
        # Create canvas
        self.canvas = FigureCanvasTkAgg(self.fig, master=plots_frame)
        self.canvas.get_tk_widget().grid(row=0, column=0, padx=5, pady=5)
        
    def set_default_values(self):
        for param, var in self.input_vars.items():
            var.set(str(self.default_values[param]))
            
    def generate(self):
        try:
            # Get input values
            E = float(self.input_vars['E'].get())
            L = float(self.input_vars['L'].get())
            rho = float(self.input_vars['rho'].get())
            twist = float(self.input_vars['twist'].get())
            layers = int(self.input_vars['layers'].get())
            
            # Run calculations
            results = run_calculations(E, L, rho)
            
            # Update results text
            self.results_text.delete(1.0, tk.END)
            self.results_text.insert(tk.END, results)
            
            # Generate plots
            self.fig.clear()
            
            # === Temporary simple plot ===
            ax = self.fig.add_subplot(111)
            ax.plot([1, 2, 3], [1, 4, 9], marker='o') # Added markers
            ax.set_title("Test Plot - Does this appear?")
            ax.grid(True)
            self.canvas.draw()
            # === End temporary plot ===

            # === Original plotting logic (commented out) ===
            # # Get blade profile points
            # X, Y = get_blade_profile()
            # # Create subplots
            # ax1 = self.fig.add_subplot(121, projection='3d')
            # ax2 = self.fig.add_subplot(122)
            # # Plot 3D blade
            # Z = np.linspace(0, 1, layers)
            # twist_per_loop = 2 * np.pi * twist / layers
            # X_orig, Y_orig = X.copy(), Y.copy() # Keep original profile for 2D plot
            # for z in Z:
            #     ax1.plot3D(X, Y, z, 'b')
            #     # Use compute_poly_AI and rotate_around_point from calcs
            #     A, Ix, Iy, Ixy, Cx, Cy = compute_poly_AI(X, Y)
            #     X, Y = rotate_around_point(X, Y, twist_per_loop * z, Cx, Cy)
            # ax1.set_aspect("equal")
            # ax1.set_title("3D Blade Visualization")
            # # Plot 2D profile
            # ax2.plot(X_orig, Y_orig, 'b-') # Use original coordinates
            # ax2.set_aspect("equal")
            # ax2.set_title("2D Blade Profile")
            # ax2.grid(True)
            # self.canvas.draw()
            # === End original logic ===
            
        except ValueError as e:
            messagebox.showerror("Error", f"Invalid input: {str(e)}")
        except Exception as e:
            # Catch other potential errors during calculation/plotting
            messagebox.showerror("Error", f"An unexpected error occurred: {str(e)}")
            print(f"Error during generate: {e}") # Print error to console
            
def main():
    root = tk.Tk()
    app = BladeCalculatorGUI(root)
    root.mainloop()

if __name__ == "__main__":
    main() 