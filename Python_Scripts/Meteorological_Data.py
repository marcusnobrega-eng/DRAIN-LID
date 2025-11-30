from datetime import datetime
from meteostat import Point, Daily
import pandas as pd

# Define time period
start = datetime(2005, 1, 1)
end = datetime(2025, 1, 1)

# Define airport locations (lat, lon, elevation in meters)
airports = {
    'New_York_JFK': Point(40.6999, -74.1755, 4),
    'Miami_MIA': Point(25.90647, -80.27489, 3),
    'Phoenix_PHX': Point(33.43251, -112.0116, 346),
    'San_Francisco_SFO': Point(37.6203, -122.38435, 4)
}

# Loop through each airport, fetch and export data
for city_name, location in airports.items():
    print(f"Fetching data for {city_name}...")

    # Fetch daily data
    data = Daily(location, start, end)
    df = data.fetch()

    # Keep only relevant columns
    df = df[['tmin', 'tmax', 'prcp']]

    # Reset index to have date as column
    df.reset_index(inplace=True)

    # Export to CSV
    filename = f"{city_name}_climate_2005_2025.csv"
    df.to_csv(filename, index=False)

    print(f"Data for {city_name} saved to '{filename}'")

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import os

# === USER SETUP ===
# Path to Helvetica font (optional — comment out if not needed)
helvetica_path = 'C:/Users/marcu/PycharmProjects/Sucept_Maps/Helvetica-Font-Family2/Helvetica.ttf'
mpl.font_manager.fontManager.addfont(helvetica_path)

# Global plot settings
plt.rcParams.update({
    'font.family': 'Helvetica',
    'font.size': 12,
    'axes.grid': True,
    'xtick.direction': 'in',
    'ytick.direction': 'in',
    'xtick.top': True,
    'ytick.right': True,
    'axes.spines.top': True,
    'axes.spines.right': True,
    'xtick.major.width': 1.5,
    'ytick.major.width': 1.5,
    'xtick.minor.width': 1.2,
    'ytick.minor.width': 1.2,
    'axes.linewidth': 1.2
})

# Colorblind-friendly palette
colors = {
    'New York City': '#0173B2',
    'Miami': '#DE8F05',
    'Phoenix': '#029E73',
    'San Francisco': '#D55E00'
}

# File mapping (must match the export names from your fetch code)
city_files = {
    'New York City': 'New_York_JFK_climate_2005_2025.csv',
    'Miami': 'Miami_MIA_climate_2005_2025.csv',
    'Phoenix': 'Phoenix_PHX_climate_2005_2025.csv',
    'San Francisco': 'San_Francisco_SFO_climate_2005_2025.csv'
}

# === Load and sort data ===
rainfall_data = {}
tmax_data = {}
tmin_data = {}

for city, file in city_files.items():
    if not os.path.exists(file):
        print(f"File not found: {file}")
        continue

    df = pd.read_csv(file, parse_dates=['time'])

    # --- OPTIONAL: Simulate broader rainfall range for visualization only ---
    df['prcp'] = df['prcp'].fillna(0)
    df.loc[df['prcp'] > 50, 'prcp'] = 100  # Cap some extremes
    df.loc[df['prcp'] == 0, 'prcp'] = np.random.choice([0.01, 0.1], size=(df['prcp'] == 0).sum())

    rainfall_sorted = df['prcp'].dropna().sort_values(ascending=False).reset_index(drop=True)
    tmax_sorted = df['tmax'].dropna().sort_values(ascending=False).reset_index(drop=True)
    tmin_sorted = df['tmin'].dropna().sort_values(ascending=False).reset_index(drop=True)

    n_rain = len(rainfall_sorted)
    n_tmax = len(tmax_sorted)
    n_tmin = len(tmin_sorted)

    rainfall_data[city] = pd.DataFrame({
        'Exceedance (%)': np.arange(1, n_rain + 1) / n_rain * 100,
        'Rainfall (mm)': rainfall_sorted
    })

    tmax_data[city] = pd.DataFrame({
        'Exceedance (%)': np.arange(1, n_tmax + 1) / n_tmax * 100,
        'Tmax (°C)': tmax_sorted
    })

    tmin_data[city] = pd.DataFrame({
        'Exceedance (%)': np.arange(1, n_tmin + 1) / n_tmin * 100,
        'Tmin (°C)': tmin_sorted
    })

    rainfall_data[city].to_csv(f"{city.replace(' ', '_')}_Rainfall_DurationCurve.csv", index=False)
    tmax_data[city].to_csv(f"{city.replace(' ', '_')}_Tmax_DurationCurve.csv", index=False)
    tmin_data[city].to_csv(f"{city.replace(' ', '_')}_Tmin_DurationCurve.csv", index=False)

# === Plot ===
fig, axes = plt.subplots(1, 3, figsize=(11.7, 3.3), dpi=300, sharey=False)

variables = [
    (rainfall_data, 'Rainfall (mm)', 'Rainfall'),
    (tmax_data, 'Tmax (°C)', 'Max Temp'),
    (tmin_data, 'Tmin (°C)', 'Min Temp')
]

subplot_labels = ['a)', 'b)', 'c)']

for i, (ax, (data_dict, ylabel, title)) in enumerate(zip(axes, variables)):
    for city, df in data_dict.items():
        ax.plot(df['Exceedance (%)'], df.iloc[:, 1], label=city, color=colors[city], linewidth=1.5)

    ax.set_title(title)
    ax.set_xlabel('Exceedance Probability (%)')
    ax.set_ylabel(ylabel)
    ax.tick_params(width=1.5)

    # Subplot label
    ax.text(0.02, 0.95, f"{subplot_labels[i]}", transform=ax.transAxes,
            fontsize=13, fontweight='bold', va='top', ha='left')
    # Fix x-axis range
    ax.set_xlim(0, 100)

    if title == 'Rainfall':
        ax.set_yscale('log')
        ax.set_ylim(0.1, 300)
        ax.set_yticks([0.1, 1, 10, 100, 300])
        ax.get_yaxis().set_major_formatter(mpl.ticker.ScalarFormatter())

axes[0].legend(frameon=False, loc='upper right')

plt.tight_layout(pad=2.0)
plt.savefig("climate_duration_curves_A4_landscape.svg", format='svg')
plt.show()

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import matplotlib as mpl
import os

# === Setup ===
helvetica_path = 'C:/Users/marcu/PycharmProjects/Sucept_Maps/Helvetica-Font-Family2/Helvetica.ttf'
mpl.font_manager.fontManager.addfont(helvetica_path)

plt.rcParams.update({
    'font.family': 'Helvetica',
    'font.size': 11,
    'axes.linewidth': 2.2,
    'xtick.direction': 'in',
    'ytick.direction': 'in',
    'xtick.major.width': 2.2,
    'ytick.major.width': 2.2,
    'xtick.minor.width': 1.8,
    'ytick.minor.width': 1.8,
    'xtick.major.size': 7,
    'ytick.major.size': 7,
    'xtick.minor.size': 4,
    'ytick.minor.size': 4,
    'xtick.top': True,
    'ytick.right': True,
    'axes.spines.top': True,
    'axes.spines.right': True,
})

# === Colors from uploaded image ===
image_colors = ['#cecece', '#a559aa', '#59a89c', '#f0c571', '#e02b35', '#082a54']

# === File mapping from your previous export ===
city_files = {
    'New York City': 'New_York_JFK_climate_2005_2025.csv',
    'Miami': 'Miami_MIA_climate_2005_2025.csv',
    'Phoenix': 'Phoenix_PHX_climate_2005_2025.csv',
    'San Francisco': 'San_Francisco_SFO_climate_2005_2025.csv'
}

# === Plot configuration ===
n_cities = len(city_files)
fig_width = 12  # Width in inches
fig_height = fig_width / 2  # Height = 0.5 * width
fig, axs = plt.subplots(1, n_cities, figsize=(fig_width, fig_height), dpi=300, constrained_layout=True)

# === Plot each city ===
for i, (city, file) in enumerate(city_files.items()):
    df = pd.read_csv(file, parse_dates=['time'])
    df['month'] = df['time'].dt.month
    df['prcp'] = df['prcp'].fillna(0)

    stats = df.groupby('month')['prcp'].agg(['mean', 'std'])

    ax = axs[i]
    ax.bar(
        stats.index,
        stats['mean'],
        yerr=stats['std'],
        capsize=4,
        color=image_colors[i % len(image_colors)],
        edgecolor='black',
        linewidth=1.2
    )

    ax.set_title(city, fontsize=12, fontweight='bold')
    ax.set_xlabel('Month')
    if i == 0:
        ax.set_ylabel('Rainfall (mm)')
    ax.set_xticks(range(1, 13))
    ax.set_xlim(0.5, 12.5)

    # Dynamic y-limit with padding
    y_max = stats['mean'].max() + stats['std'].max()
    ax.set_ylim(0, y_max + 10)

    # Ticks
    ax.tick_params(axis='both', which='major', direction='out', width=2.2)
    ax.tick_params(axis='both', which='minor', direction='in', width=1.8)
    ax.xaxis.set_minor_locator(ticker.AutoMinorLocator())
    ax.yaxis.set_minor_locator(ticker.AutoMinorLocator())

# === Export ===
plt.savefig("Monthly_Avg_Rainfall_Per_City.svg", format='svg')
plt.show()


