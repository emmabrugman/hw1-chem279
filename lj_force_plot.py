import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import linregress

# Data 
h_values = np.array([0.1, 0.01, 0.001, 0.0001])
log_h_values = np.log(h_values)
log_forward_error = np.array([-1.45477, -1.47962, -1.48247, -1.48276])
log_central_error = np.array([-23.0259, -23.0259, -23.0259, -23.0259])

# Calculate slopes
slope_forward, intercept_forward, _, _, _ = linregress(log_h_values, log_forward_error)
slope_central, intercept_central, _, _, _ = linregress(log_h_values, log_central_error)

# Plot data
plt.figure(figsize=(10, 6))
plt.plot(log_h_values, log_forward_error, 'o-', label=f'Forward Difference (Slope = {slope_forward:.2f})', markersize=8)
plt.plot(log_h_values, log_central_error, 'o-', label=f'Central Difference (Slope = {slope_central:.2f})', markersize=8)
plt.xlabel('Log(Step Size)')
plt.ylabel('Log(Error)')
plt.title('Log(Error) vs Log(Step Size) and Slopes')
plt.legend()
plt.grid(True)
plt.show()
