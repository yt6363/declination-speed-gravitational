import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from astropy.time import Time
from astropy.coordinates import solar_system_ephemeris, get_body_barycentric, get_sun, get_moon, SkyCoord, GCRS
from astropy.constants import G, M_sun, M_earth
import certifi
import os

def orbital_speed(distance, mu):
    return np.sqrt(mu / distance)

def gravitational_force(mass, distance):
    return (G.value * M_sun.value * mass) / (distance ** 2)

def get_declination(obj, time):
    if obj == 'sun':
        obj_position = get_sun(time)
    elif obj == 'moon':
        obj_position = get_moon(time)
    else:
        with solar_system_ephemeris.set('jpl'):
            obj_position = get_body_barycentric(obj, time)
    
    # Convert to SkyCoord for proper transformations
    sky_coord = SkyCoord(obj_position, frame='gcrs', obstime=time)
    return sky_coord.icrs.dec.degree

def get_object_position(obj, time):
    if obj == 'sun':
        return get_sun(time)
    elif obj == 'moon':
        return get_moon(time)
    else:
        with solar_system_ephemeris.set('jpl'):
            return get_body_barycentric(obj, time)

# Set SSL_CERT_FILE to use certifi's certificate bundle
os.environ['SSL_CERT_FILE'] = certifi.where()

# User input for objects, duration, and start date
objects = input("Enter the object names separated by commas (e.g., 'mars,jupiter,sun,moon'): ").strip().lower().split(',')
duration_unit = input("Enter the duration unit (years, months, weeks, days): ").strip().lower()
duration_value = int(input(f"Enter the number of {duration_unit}: "))
start_date = input("Enter the start date (YYYY-MM-DD): ").strip()

# Convert duration to days
if duration_unit == 'years':
    duration_days = duration_value * 365
elif duration_unit == 'months':
    duration_days = duration_value * 30
elif duration_unit == 'weeks':
    duration_days = duration_value * 7
elif duration_unit == 'days':
    duration_days = duration_value
else:
    raise ValueError("Invalid duration unit. Please enter 'years', 'months', 'weeks', or 'days'.")

# Set the observation start time
observation_start_time = Time(start_date)

# Create an array of observation times over the specified duration
observation_times = Time(np.linspace(observation_start_time.jd, observation_start_time.jd + duration_days, duration_days), format='jd')

# Masses of the planets and Moon (in kg)
object_masses = {
    'mercury': 3.3011e23,
    'venus': 4.8675e24,
    'earth': 5.97237e24,
    'mars': 6.4171e23,
    'jupiter': 1.8982e27,
    'saturn': 5.6834e26,
    'uranus': 8.6810e25,
    'neptune': 1.02413e26,
    'moon': M_earth.value * 0.0123,  # Approximate mass of the Moon
    'sun': M_sun.value
}

# Calculate and plot normalized gravitational forces and declination for each object
plt.figure(figsize=(15, 10))

# Gravitational force plot
plt.subplot(2, 1, 1)
for obj in objects:
    if obj not in object_masses:
        print(f"Unknown object: {obj}")
        continue
    
    mass = object_masses[obj]
    forces = []
    for t in observation_times:
        obj_position = get_object_position(obj, t)
        distance = obj_position.norm().to('m').value
        force = gravitational_force(mass, distance)
        forces.append(force)
    
    # Normalize forces
    forces = np.array(forces)
    normalized_forces = (forces - np.min(forces)) / (np.max(forces) - np.min(forces))
    
    # Plot the normalized forces with respect to the observation times
    plt.plot(observation_times.plot_date, normalized_forces, label=obj.capitalize())

# Customize the plot
plt.xlabel('Date')
plt.ylabel('Normalized Gravitational Force')
plt.title(f'Normalized Gravitational Forces from {start_date} over {duration_value} {duration_unit.capitalize()}')
plt.grid(True)

# Set date format on the x-axis to show year and month
date_format = mdates.DateFormatter('%Y-%m')
plt.gca().xaxis.set_major_formatter(date_format)

# Add a legend to differentiate the objects
plt.legend()

# Enable the display of the x, y coordinates on hover
def format_coord(x, y):
    date = mdates.num2date(x)
    return f'{date:%Y-%m-%d}, {y:.2f}'

plt.gca().format_coord = format_coord

# Declination plot
plt.subplot(2, 1, 2)
for obj in objects:
    declinations = []
    for t in observation_times:
        dec = get_declination(obj, t)
        declinations.append(dec)
    
    # Normalize declinations
    declinations = np.array(declinations)
    normalized_declinations = (declinations - np.min(declinations)) / (np.max(declinations) - np.min(declinations))
    
    # Plot the normalized declinations with respect to the observation times
    plt.plot(observation_times.plot_date, normalized_declinations, label=obj.capitalize())

# Customize the plot
plt.xlabel('Date')
plt.ylabel('Normalized Declination')
plt.title(f'Normalized Declinations from {start_date} over {duration_value} {duration_unit.capitalize()}')
plt.grid(True)

# Set date format on the x-axis to show year and month
plt.gca().xaxis.set_major_formatter(date_format)

# Add a legend to differentiate the objects
plt.legend()

# Enable the display of the x, y coordinates on hover
plt.gca().format_coord = format_coord

# Adjust layout and show the plot
plt.tight_layout()
plt.show()
