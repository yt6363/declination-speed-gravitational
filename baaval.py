import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import mplcursors
from astropy.time import Time
from astropy.coordinates import solar_system_ephemeris, get_body, get_sun, get_moon, SkyCoord, GCRS, get_body_barycentric_posvel
from astropy.constants import G, M_sun, M_earth
import tkinter as tk
from tkinter import ttk
import certifi
import os

def orbital_speed(distance, mu):
    return np.sqrt(mu / distance)

def gravitational_force(mass1, mass2, distance):
    return (G.value * mass1 * mass2) / (distance ** 2)

def get_declination(obj, time):
    if obj == 'sun':
        obj_position = get_sun(time)
    elif obj == 'moon':
        obj_position = get_moon(time)
    else:
        with solar_system_ephemeris.set('jpl'):
            obj_position = get_body(obj, time)
    
    sky_coord = SkyCoord(obj_position)
    gcrs_coord = sky_coord.transform_to(GCRS(obstime=time))
    return gcrs_coord.dec.degree

def get_object_position(obj, time):
    if obj == 'sun':
        return get_sun(time)
    elif obj == 'moon':
        return get_moon(time)
    else:
        with solar_system_ephemeris.set('jpl'):
            return get_body(obj, time)

def get_distance(body_name, time):
    earth = get_body_barycentric_posvel('earth', time)
    body = get_body_barycentric_posvel(body_name, time)
    distance = (body[0] - earth[0]).norm().to('m').value
    return distance

# Set SSL_CERT_FILE to use certifi's certificate bundle
os.environ['SSL_CERT_FILE'] = certifi.where()

def plot_data(objects, duration_unit, duration_value, start_date, plot_speed, plot_grav_force, plot_declination, plot_grav_on_earth):
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

    mu = G.value * M_sun.value  # Gravitational parameter for the Sun

    # Create the plot figure
    plt.figure(figsize=(15, 20))

    # Plot orbital speed
    if plot_speed:
        ax1 = plt.subplot(4, 1, 1)
        for obj in objects:
            if obj not in object_masses:
                print(f"Unknown object: {obj}")
                continue

            speeds = []
            for t in observation_times:
                obj_position = get_object_position(obj, t)
                distance = obj_position.distance.to('m').value
                speed = orbital_speed(distance, mu)
                speeds.append(speed)

            speeds = np.array(speeds)
            normalized_speeds = (speeds - np.min(speeds)) / (np.max(speeds) - np.min(speeds))
            plt.plot(observation_times.plot_date, normalized_speeds, label=obj.capitalize())

        plt.xlabel('Date')
        plt.ylabel('Normalized Orbital Speed')
        plt.title(f'Normalized Orbital Speeds from {start_date} over {duration_value} {duration_unit.capitalize()}')
        plt.grid(True)
        date_format = mdates.DateFormatter('%Y-%m-%d')
        plt.gca().xaxis.set_major_formatter(date_format)
        plt.legend()

        # Enable mplcursors
        cursor1 = mplcursors.cursor(ax1, hover=True)
        cursor1.connect("add", lambda sel: sel.annotation.set_text(f'{mdates.num2date(sel.target[0]):%Y-%m-%d}, {sel.target[1]:.2f}'))

    # Plot gravitational force
    if plot_grav_force:
        ax2 = plt.subplot(4, 1, 2)
        for obj in objects:
            if obj not in object_masses:
                print(f"Unknown object: {obj}")
                continue

            mass = object_masses[obj]
            forces = []
            for t in observation_times:
                obj_position = get_object_position(obj, t)
                distance = obj_position.distance.to('m').value
                force = gravitational_force(mass, M_sun.value, distance)
                forces.append(force)

            forces = np.array(forces)
            normalized_forces = (forces - np.min(forces)) / (np.max(forces) - np.min(forces))
            plt.plot(observation_times.plot_date, normalized_forces, label=obj.capitalize())

        plt.xlabel('Date')
        plt.ylabel('Normalized Gravitational Force')
        plt.title(f'Normalized Gravitational Forces from {start_date} over {duration_value} {duration_unit.capitalize()}')
        plt.grid(True)
        date_format = mdates.DateFormatter('%Y-%m-%d')
        plt.gca().xaxis.set_major_formatter(date_format)
        plt.legend()

        # Enable mplcursors
        cursor2 = mplcursors.cursor(ax2, hover=True)
        cursor2.connect("add", lambda sel: sel.annotation.set_text(f'{mdates.num2date(sel.target[0]):%Y-%m-%d}, {sel.target[1]:.2f}'))

    # Plot declination
    if plot_declination:
        ax3 = plt.subplot(4, 1, 3)
        for obj in objects:
            declinations = []
            for t in observation_times:
                dec = get_declination(obj, t)
                declinations.append(dec)

            declinations = np.array(declinations)
            normalized_declinations = (declinations - np.min(declinations)) / (np.max(declinations) - np.min(declinations))
            plt.plot(observation_times.plot_date, normalized_declinations, label=obj.capitalize())

        plt.xlabel('Date')
        plt.ylabel('Normalized Declination')
        plt.title(f'Normalized Declinations from {start_date} over {duration_value} {duration_unit.capitalize()}')
        plt.grid(True)
        date_format = mdates.DateFormatter('%Y-%m-%d')
        plt.gca().xaxis.set_major_formatter(date_format)
        plt.legend()

        # Enable mplcursors
        cursor3 = mplcursors.cursor(ax3, hover=True)
        cursor3.connect("add", lambda sel: sel.annotation.set_text(f'{mdates.num2date(sel.target[0]):%Y-%m-%d}, {sel.target[1]:.2f}'))

    # Plot gravitational force on Earth
    if plot_grav_on_earth:
        ax4 = plt.subplot(4, 1, 4)
        for planet in object_masses:
            if planet == 'earth':
                continue
            
            forces = []
            for t in observation_times:
                distance = get_distance(planet, t)
                force = gravitational_force(object_masses[planet], M_earth.value, distance)
                forces.append(force)
            
            forces = np.array(forces)
            normalized_forces = (forces - np.min(forces)) / (np.max(forces) - np.min(forces))
            plt.plot(observation_times.plot_date, normalized_forces, label=planet.capitalize())
        
        plt.xlabel('Date')
        plt.ylabel('Normalized Gravitational Force on Earth')
        plt.title(f'Normalized Gravitational Forces Acting on Earth from {start_date} over {duration_value} {duration_unit.capitalize()}')
        plt.grid(True)
        date_format = mdates.DateFormatter('%Y-%m-%d')
        plt.gca().xaxis.set_major_formatter(date_format)
        plt.legend()

        # Enable mplcursors
        cursor4 = mplcursors.cursor(ax4, hover=True)
        cursor4.connect("add", lambda sel: sel.annotation.set_text(f'{mdates.num2date(sel.target[0]):%Y-%m-%d}, {sel.target[1]:.2f}'))

    plt.tight_layout()
    plt.show()

# GUI for user input
class UserInputDialog(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("Orbital and Gravitational Plotter")

        self.selected_plots = {
            "speed": tk.BooleanVar(value=True),
            "grav_force": tk.BooleanVar(value=True),
            "declination": tk.BooleanVar(value=True),
            "grav_on_earth": tk.BooleanVar(value=True)
        }

        self.create_widgets()

    def create_widgets(self):
        # Objects
        ttk.Label(self, text="Enter the object names (e.g., 'mars,jupiter,sun,moon')").grid(row=0, column=0, padx=10, pady=5, sticky="w")
        self.objects_entry = ttk.Entry(self)
        self.objects_entry.grid(row=0, column=1, padx=10, pady=5)

        # Duration unit
        ttk.Label(self, text="Enter the duration unit (years, months, weeks, days)").grid(row=1, column=0, padx=10, pady=5, sticky="w")
        self.duration_unit_entry = ttk.Entry(self)
        self.duration_unit_entry.grid(row=1, column=1, padx=10, pady=5)

        # Duration value
        ttk.Label(self, text="Enter the number of duration units").grid(row=2, column=0, padx=10, pady=5, sticky="w")
        self.duration_value_entry = ttk.Entry(self)
        self.duration_value_entry.grid(row=2, column=1, padx=10, pady=5)

        # Start date
        ttk.Label(self, text="Enter the start date (YYYY-MM-DD)").grid(row=3, column=0, padx=10, pady=5, sticky="w")
        self.start_date_entry = ttk.Entry(self)
        self.start_date_entry.grid(row=3, column=1, padx=10, pady=5)

        # Plot selection
        ttk.Label(self, text="Select plots to display").grid(row=4, column=0, padx=10, pady=5, sticky="w")
        ttk.Checkbutton(self, text="Orbital Speed", variable=self.selected_plots["speed"]).grid(row=4, column=1, padx=10, pady=5, sticky="w")
        ttk.Checkbutton(self, text="Gravitational Force", variable=self.selected_plots["grav_force"]).grid(row=5, column=1, padx=10, pady=5, sticky="w")
        ttk.Checkbutton(self, text="Declination", variable=self.selected_plots["declination"]).grid(row=6, column=1, padx=10, pady=5, sticky="w")
        ttk.Checkbutton(self, text="Gravitational Force on Earth", variable=self.selected_plots["grav_on_earth"]).grid(row=7, column=1, padx=10, pady=5, sticky="w")

        # Submit button
        self.submit_button = ttk.Button(self, text="Submit", command=self.submit)
        self.submit_button.grid(row=8, column=0, columnspan=2, pady=10)

    def submit(self):
        objects = self.objects_entry.get().strip().lower().split(',')
        duration_unit = self.duration_unit_entry.get().strip().lower()
        duration_value = int(self.duration_value_entry.get().strip())
        start_date = self.start_date_entry.get().strip()

        plot_speed = self.selected_plots["speed"].get()
        plot_grav_force = self.selected_plots["grav_force"].get()
        plot_declination = self.selected_plots["declination"].get()
        plot_grav_on_earth = self.selected_plots["grav_on_earth"].get()

        self.destroy()

        plot_data(objects, duration_unit, duration_value, start_date, plot_speed, plot_grav_force, plot_declination, plot_grav_on_earth)

# Run the GUI
if __name__ == "__main__":
    app = UserInputDialog()
    app.mainloop()
