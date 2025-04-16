import sys
import math
import numpy as np
import matplotlib.pyplot as plt
import naca # Import the naca module
from calcs import run_calculations, compute_poly_AI, rotate_around_point, Euler_Bournoulli_Clamped_Free_mode_shape_function

def main():
	# Define input parameters
	E = 650000000  # Young's modulus (N/m^2)
	L = 0.05       # Length (m)
	rho = 1710     # Density (kg/m^3)
	twist = 0.1    # Twist factor (used for visualization)
	layers = 25     # Number of layers for 3D plot visualization
	profile_number = '2412'
	n_points_per_side = 300 # Leads to 2*n+1 = 601 total points

	# --- Generate Profile --- 
	print(f"Generating NACA {profile_number} profile...")
	X_base, Y_base = naca.naca(profile_number, n_points_per_side)
	X_base = np.array(X_base)
	Y_base = np.array(Y_base)

	print("Running Blade Calculations...")

	# --- Perform Calculations & Print Results ---
	# Get both the results string and the frequency data dictionary
	results_string, freq_data = run_calculations(X_base, Y_base, E, L, rho)
	print(results_string)

	print("Generating Plots...")

	# --- Generate Plots ---
	fig, axes = plt.subplots(2, 2, figsize=(12, 10))
	fig.suptitle('Blade Analysis Plots', fontsize=16)

	# --- Plot 1: 2D Profile (Top-Left) ---
	ax1 = axes[0, 0]
	ax1.plot(X_base, Y_base, 'b-', label='Profile')
	
	# Calculate and plot the centroid
	_, _, _, _, Cx, Cy = compute_poly_AI(X_base, Y_base)
	ax1.plot(Cx, Cy, 'ro', markersize=8, label='Centroid (CoM)') # 'ro' for red circle
	
	ax1.set_title(f'2D Blade Profile (NACA {profile_number})') # Update title
	ax1.set_xlabel('Chord Fraction')
	ax1.set_ylabel('Thickness / Camber Fraction') # More accurate label
	ax1.set_aspect('equal', adjustable='box')
	ax1.grid(True)
	
	# Expand y-axis limits for legend space (increase expansion factor)
	current_ylim = ax1.get_ylim()
	y_center = (current_ylim[1] + current_ylim[0]) / 2
	y_range = current_ylim[1] - current_ylim[0]
	# Apply a 4x expansion relative to the original auto-range (2 * previous doubled range)
	ax1.set_ylim(y_center - y_range * 2, y_center + y_range * 2)

	# Position legend inside the plot area, top right
	ax1.legend(loc='upper right')

	# --- Plot 2: 3D Twisted Blade (Top-Right) ---
	ax2 = fig.add_subplot(2, 2, 2, projection='3d') # Use fig.add_subplot for 3D
	
	Z_coords = np.linspace(0, L, layers) # Use actual length L for z-coordinates
	# Calculate twist per layer based on total twist angle (e.g., twist * pi)
	total_twist_angle = twist * np.pi # Example: twist=0.1 -> 0.1*pi total twist
	twist_per_layer = total_twist_angle / (layers - 1) if layers > 1 else 0

	X_current, Y_current = X_base.copy(), Y_base.copy()
	# Compute initial centroid for rotation
	A_init, _, _, _, Cx_init, Cy_init = compute_poly_AI(X_current, Y_current)
	
	# Store initial coords for plot limits
	all_x = [X_current]
	all_y = [Y_current]
	
	for i, z_val in enumerate(Z_coords):
		ax2.plot(X_current, Y_current, zs=z_val, zdir='z', color='b')
		if i < layers - 1: # Don't rotate after plotting the last layer
			 # Rotate around the initial centroid for consistent axis
			X_current, Y_current = rotate_around_point(X_current, Y_current, twist_per_layer, Cx_init, Cy_init)
			all_x.append(X_current)
			all_y.append(Y_current)

	ax2.set_title('3D Twisted Blade Visualization')
	ax2.set_xlabel('X')
	ax2.set_ylabel('Y')
	ax2.set_zlabel('Z (Length)')
	
	# Set consistent limits for 3D plot based on rotated profiles
	min_x, max_x = np.min(np.concatenate(all_x)), np.max(np.concatenate(all_x))
	min_y, max_y = np.min(np.concatenate(all_y)), np.max(np.concatenate(all_y))
	range_x = max_x - min_x
	range_y = max_y - min_y
	mid_x, mid_y = (max_x + min_x)/2, (max_y + min_y)/2
	plot_range = max(range_x, range_y, L) # Make bounding box roughly cubic
	ax2.set_xlim(mid_x - plot_range/2, mid_x + plot_range/2)
	ax2.set_ylim(mid_y - plot_range/2, mid_y + plot_range/2)
	ax2.set_zlim(0, L)
	# ax2.set_aspect('equal') # Often tricky with 3D plots, using manual limits

	# --- Plot 3: Mode Shapes (Bottom-Left) ---
	ax3 = axes[1, 0]
	num_modes_to_plot = 4
	x_mode = np.linspace(0, L, 100) # Points along the length to plot shapes
	
	for mode_index in range(num_modes_to_plot):
		# Mode numbers typically start from 1, but function might use 0-based index?
		# Let's assume function expects 1-based index for mode_number.
		# Check calcs.py: Euler_Bournoulli_Clamped_Free_mode_shape_function uses Bnl[mode_number]
		# Bnl likely 0-indexed, but the function itself seems to take mode_number (1, 2, 3...)
		# Let's pass mode_number = mode_index + 1 
		mode_number = mode_index + 1
		# Now call the directly imported function
		y_mode = Euler_Bournoulli_Clamped_Free_mode_shape_function(L, mode_number, x_mode)
		ax3.plot(x_mode, y_mode, label=f'Mode {mode_number}')

	# ax3.text(0.5, 0.5, 'Plot Area 3 (Empty)', ha='center', va='center', fontsize=12, alpha=0.5)
	# ax3.set_xticks([])
	# ax3.set_yticks([])
	# ax3.spines['top'].set_visible(False)
	# ax3.spines['right'].set_visible(False)
	# ax3.spines['bottom'].set_visible(False)
	# ax3.spines['left'].set_visible(False)
	ax3.set_title('Mode Shapes')
	ax3.set_xlabel('Position along Length (m)')
	ax3.set_ylabel('Relative Amplitude')
	ax3.grid(True)
	ax3.legend()

	# --- Plot 4: Summary Text (Bottom-Right) ---
	ax4 = axes[1, 1]
	ax4.set_xticks([])
	ax4.set_yticks([])
	ax4.spines['top'].set_visible(False)
	ax4.spines['right'].set_visible(False)
	ax4.spines['bottom'].set_visible(False)
	ax4.spines['left'].set_visible(False)

	# Construct the summary text using f-strings and data
	# Using polynomial results as they are generally the most accurate for the profile
	summary_text = (
		f"Input Summary:\n"
		f"-----------------\n"
		f"Profile: NACA {profile_number}\n"
		f"E: {E:.2e} N/m^2\n"
		f"L: {L:.3f} m\n"
		f"rho: {rho:.1f} kg/m^3\n"
		f"Twist Factor: {twist}\n"
		f"\n"
		f"Modal Frequencies (Poly, Hz):\n"
		f"------------------------------\n"
		f"Lumped: {freq_data['lumped_poly']:.2f}\n"
		f"Dist. M1: {freq_data['dist_poly'][0]:.2f}\n"
		f"Dist. M2: {freq_data['dist_poly'][1]:.2f}\n"
		f"Dist. M3: {freq_data['dist_poly'][2]:.2f}\n"
		f"Dist. M4: {freq_data['dist_poly'][3]:.2f}"
	)

	# Add text to the axes
	ax4.text(0.05, 0.95, summary_text, 
			 ha='left', va='top', fontsize=10, 
			 wrap=True, # Allow text wrapping
			 family='monospace') # Use monospace font for alignment

	# --- Final Adjustments & Display ---
	plt.tight_layout(rect=[0, 0.03, 1, 0.95])
	plt.show()

	print("Done.")

if __name__ == '__main__':
	main()
	sys.exit()



