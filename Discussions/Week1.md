# Generation of Skyplots

Given the time (and other data) and user position, calculate the position of the satellites in GNSS constellations and their azimuth and elevation from certain user position.

**Inputs**
1. Time: 2024-01-29:16:57:23+0300
2. User Position: 22.3085316N, 39.104593E

### 1. Convert Time to UTC then to GPS Time (Week)
how to convert:
converting from a certain timezone to UTC is just subtracting the timezone. Note that, if the result is negative we have to add 24 to get the time and subtract a day.
$$t_{utc} = t_{local} - t_{z}$$

```python
from datetime import timezone
dt = datetime(2015, 10, 19)
timestamp = dt.replace(tzinfo=timezone.utc).timestamp()
print(timestamp)
```
### 2. Obtain YUMA Files for the GPS Week and Read Yuma Files

### 3. Calculate the True anomaly from the mean anomaly
$$M = M0 + n(t - t_p) $$
Obtain $E$ from the following equation:
$$M = E - e \sin(E)$$ 
solve by Newton-Raphson
Use $E$ to calculate the true anomaly $\nu$
$$\nu = \tan^{-1}\left[ { \frac{\sqrt{1 - e^2} \sin{E}}{\cos{E} - e}} \right] $$ 


### 4. Calculate the position of the satellites (at least GPS)

### 5. Convert User positon to WGS-84 - XYZ format (meters)

### 6. Plot the User position (on sphere) and Satellite positions 

### 7. Calculate the pointing vector from the User to satellites and its projections

### 8. Create Skyplot