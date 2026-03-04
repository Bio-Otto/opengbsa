
import numpy as np
val = np.float64(1.0)
print(f"Value: {val}, Type: {type(val)}")
print(f"Is float? {isinstance(val, float)}")
print(f"Is int/float? {isinstance(val, (int, float))}")
print(f"Is np.number? {isinstance(val, np.number)}")
