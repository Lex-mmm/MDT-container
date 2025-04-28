
import numpy as np
data = {'time': np.float64(1.0), 'values': {'heart_rate': np.float64(68.95), 'SaO2': np.float64(97.0), 'MAP': np.float64(10.28), 'RR': np.float64(11.99), 'etCO2': np.float64(37.3)}}
print(data['values']['SaO2'])