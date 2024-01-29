# Generation of Skyplots

Given the time (and other data) and user position, calculate the position of the satellites in GNSS constellations and their azimuth and elevation from certain user position.

**Inputs**
1. Time: 2024-01-29:16:57:23+0300
2. User Position: 22.3085316N, 39.104593E

### 1. Convert Time to UTC then to GPS Time (Week)

### 2. Obtain YUMA Files for the GPS Week and Read Yuma Files

### 3. Calculate the True anomaly from the mean anomaly

### 4. Calculate the position of the satellites (at least GPS)

### 5. Convert User positon to WGS-84 - XYZ format (meters)

### 6. Plot the User position (on sphere) and Satellite positions 

### 7. Calculate the pointing vector from the User to satellites and its projections

### 8. Create Skyplot